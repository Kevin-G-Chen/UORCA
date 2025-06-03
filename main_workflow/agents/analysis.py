#############################
# SECTION: Logging and Utilities (updated)
#############################
from __future__ import annotations
import logging, os, re, subprocess, json, glob, argparse, asyncio, pathlib, datetime
from typing import List, Optional, Dict, Any, Union, Tuple, Literal

import pandas as pd
import numpy as np
from dataclasses import dataclass
from dotenv import load_dotenv
from pydantic import BaseModel, Field, ConfigDict
from pydantic_ai import Agent, RunContext
from shared import AnalysisContext, CheckpointStatus
from shared.workflow_logging import log_tool, log_agent_tool, log_tool_for_reflection
from agents.metadata import metadata_agent, MetadataContext
from unidecode import unidecode
import datetime
from openai import OpenAI
import matplotlib.pyplot as plt
import nest_asyncio
from concurrent.futures import ThreadPoolExecutor, as_completed

# ── configure python‑logging ────────────────────────────────────────────────

logger = logging.getLogger(__name__)

#############################
# SECTION: Custom Exceptions
#############################

class WorkflowError(Exception):
    """Base exception for workflow-related errors that agents can understand and recover from."""
    pass

class PrerequisiteError(WorkflowError):
    """Raised when a required prerequisite step has not been completed."""
    pass

class FileNotFoundError(WorkflowError):
    """Raised when a required file is missing."""
    pass

class ValidationError(WorkflowError):
    """Raised when data validation fails."""
    pass

#############################
# SECTION: Data classes and env
#############################

load_dotenv()
openai_api_key = os.getenv("OPENAI_API_KEY")
client = OpenAI()

#############################
# SECTION: Agent definitions
#############################
# ––– RNA‑seq analysis agent –––
analysis_prompt_path = "./main_workflow/prompts/analysis.txt"
try:
    system_prompt = pathlib.Path(analysis_prompt_path).read_text()
except Exception as e:
    logger.warning("Could not read analysis system prompt: %s – using fallback", e)
    system_prompt = """
    # RNA‑seq Data Analysis Expert (fallback)
    You are an expert RNA‑seq data analyst. Analyse RNA‑seq data using standard bioinformatics workflows while logging each major step.
    """

rnaseq_agent = Agent(
    "openai:o4-mini",
    deps_type=AnalysisContext,
    system_prompt=system_prompt)

# ----------------------------
# Define output schema for contrasts
# ----------------------------
class Contrast_format(BaseModel):
    """Schema for each contrast."""
    name: str
    expression: str
    description: Optional[str] = Field(
        description="Biological interpretation of the contrast")
    justification: Optional[str] = Field(
        description="Justification for the contrast design, in terms of value to the scientific community")
class Contrasts(BaseModel):
    """Schema for the output of the candidate contrasts."""
    contrasts: List[Contrast_format]
    model_config = ConfigDict(extra="allow")

# ----------------------------
# Create an RNAseq analysis agent
# ----------------------------

# Read system prompt from file
data_analysis_agent_system_prompt_path = "./main_workflow/prompts/analysis.txt"
metadata_prompt_path = "./main_workflow/prompts/metadata.txt"

try:
    with open(data_analysis_agent_system_prompt_path, 'r') as f:
        system_prompt = f.read()

    # Print the first and last few lines for verification
    lines = system_prompt.split('\n')
    print(f"\n--- SYSTEM PROMPT FILE LOADED ---")
    print("First 3 lines:")
    for line in lines[:3]:
        print(f"> {line}")

# If the file is not found, use a fallback system prompt

except Exception as e:
    print(f"Failed to read system prompt file: {str(e)}")
    system_prompt = """
    # RNA-seq Data Analysis Expert

    You are an expert RNAseq data analyst. Your task is to analyze RNAseq data using a series of bioinformatics tools.

    ## General Analysis Principles

    1. Work systematically through the RNA-seq analysis workflow
    2. Validate inputs at each step
    3. Provide clear explanations of what's happening
    4. Handle errors gracefully
    5. Generate appropriate visualizations when needed
    6. Be comprehensive in your analysis
    7. Track output files carefully
    """
    print("Using fallback system prompt instead.")

rnaseq_agent = Agent(
    'openai:gpt-4o',  # Change to more powerful model
    deps_type=AnalysisContext,
    system_prompt=system_prompt
)

# ----------------------------
# Utility Functions
# ----------------------------
@rnaseq_agent.tool
@log_tool
@log_tool_for_reflection
async def list_files(ctx: RunContext[AnalysisContext],
                     directory: Optional[str] = None,
                     pattern: Optional[str] = None,
                     recursive: bool = True) -> List[str]:
    """
    Search for and list files in a directory matching a specified pattern.
    Found files are automatically added to ctx.deps.files for future reference.
    In general, do not assume directory locations exist, unless they are explicitly specified in the context. Do not perform generic searches, such as in "/" or "~", or for the "*" pattern - these will be hugely problematic.

    Parameters:
      - directory (str, optional):
          The directory to search in. If not provided, uses the current directory.
      - pattern (str, optional):
          File pattern to match (e.g., "*.txt", "*.fastq.gz").
          Can be a glob pattern or file extension.
          If not provided, matches all files.
      - recursive (bool, default=True):
          Whether to search recursively through subdirectories.

    Returns:
      A list of file paths matching the criteria.
    """
    # Resolve search scope
    search_dir = directory or "."
    pat = pattern or "*"

    # SAFETY CHECKS: Validate the search parameters before proceeding

    # Check 1: Dangerous root directories
    dangerous_dirs = ['/', '/etc', '/usr', '/bin', '/home']
    if search_dir in dangerous_dirs:
        error_msg = f"Safety Error: Searching in '{search_dir}' is not allowed for system safety. Please specify a subdirectory within your working environment."
        logger.error("⚠️ %s", error_msg)
        return [error_msg]

    # Check 2: Home directory
    if search_dir.startswith('~'):
        error_msg = f"Safety Error: Searching in home directory '{search_dir}' is not allowed. Please specify a project-specific directory."
        logger.error("⚠️ %s", error_msg)
        return [error_msg]

    # Check 4: Too generic pattern with wide search
    if pat == "*" and (search_dir == "." or recursive):
        error_msg = f"Safety Error: Too generic search pattern '{pat}' in {search_dir} with recursive={recursive}. Please specify a more specific pattern."
        logger.error("⚠️ %s", error_msg)
        return [error_msg]

    logger.info("🔍 Searching for '%s' in %s (recursive=%s)",
                pat, os.path.abspath(search_dir), recursive)

    # Check if directory exists
    if not os.path.isdir(search_dir):
        msg = f"Error: Directory '{search_dir}' does not exist"
        logger.error("❌ %s", msg)
        return [msg]

    # Find files matching pattern
    try:
        # Use glob to find files
        glob_pattern = os.path.join(search_dir, "**", pat) if recursive else os.path.join(search_dir, pat)
        matched_files = sorted(glob.glob(glob_pattern, recursive=recursive))

        # Update context with found files
        if matched_files:
            # Initialize files list if it doesn't exist
            if ctx.deps.files is None:
                ctx.deps.files = []

            # Add new files to the existing list (avoiding duplicates)
            ctx.deps.files = list(set(ctx.deps.files).union(matched_files))

            # Log results
            logger.info("💾 Added %d files to context (total %d)",
                        len(matched_files), len(ctx.deps.files))

            # Show preview of found files
            preview = matched_files[:5] if len(matched_files) > 5 else matched_files
            logger.info("✅ Found %d files%s",
                        len(matched_files),
                        f": {preview}" if preview else "")
        else:
            logger.warning("⚠️ No files found for given pattern")

        return matched_files

    except Exception as e:
        error_msg = f"Error searching for files: {str(e)}"
        logger.error("❌ %s", error_msg)
        return [error_msg]



@rnaseq_agent.tool
@log_tool
@log_tool_for_reflection
async def run_kallisto_quantification(ctx: RunContext[AnalysisContext],
                                     kallisto_index: Optional[str] = None) -> str:
    """
    Run Kallisto quantification on paired-end FASTQ files and record the resulting abundance files.

    Uses the agent's intelligence to select an appropriate Kallisto index from previously
    discovered files (ctx.deps.files) based on the organism, or accepts a direct path.

    Inputs:
      - ctx: RunContext[AnalysisContext]
          Must include:
             • fastq_dir: The directory where FASTQ files are stored.
             • resource_dir: The directory containing resources.
             • output_dir: The directory to store Kallisto outputs.
      - kallisto_index: Optional[str]
          Direct path to the Kallisto index to use. If not provided, will intelligently
          select an appropriate index from ctx.deps.files based on the organism.
    """
    try:
        logger.info("🧬 run_kallisto_quantification started")
        
        # Update checkpoint: kallisto index selection started
        if hasattr(ctx.deps, 'checkpoints') and ctx.deps.checkpoints:
            ctx.deps.checkpoints.kallisto_index_selection.status = CheckpointStatus.IN_PROGRESS
            ctx.deps.checkpoints.kallisto_index_selection.timestamp = datetime.datetime.now().isoformat()
            logger.info("🔍 Checkpoint: Kallisto index selection started")

        # Validate fastq_dir is set and exists
        if not ctx.deps.fastq_dir:
            msg = "Error: No FASTQ directory specified in ctx.deps.fastq_dir"
            logger.error("❌ %s", msg)
            raise PrerequisiteError(msg)

        if not os.path.isdir(ctx.deps.fastq_dir):
            msg = f"Error: FASTQ directory {ctx.deps.fastq_dir} does not exist"
            logger.error("❌ %s", msg)
            raise FileNotFoundError(msg)

        # Always explicitly search in the fastq_dir for FASTQ files, ignoring any previously discovered files
        logger.info("🔍 Searching for FASTQ files in: %s", ctx.deps.fastq_dir)
        fastq_files = await list_files(ctx, directory=ctx.deps.fastq_dir, pattern="*.fastq.gz")

        if not fastq_files or len(fastq_files) == 0 or (isinstance(fastq_files[0], str) and fastq_files[0].startswith("Error")):
            msg = f"Error: No FASTQ files found in {ctx.deps.fastq_dir}"
            logger.error("❌ %s", msg)
            raise FileNotFoundError(msg)

        # Get Kallisto index - either provided directly or intelligently selected
        index_path = kallisto_index

        if not index_path or not os.path.exists(index_path):
            msg = (
                f"No Kallisto index provided. Please use list_files to find a suitable index file in {ctx.deps.resource_dir}. Note that the index file will end with the .idx suffix. Also note that you should use the index file that is appropriate for {ctx.deps.organism}."
            )
            logger.error("❌ %s", msg.replace('\n', ' | '))
            
            # Update checkpoint: kallisto index selection failed
            if hasattr(ctx.deps, 'checkpoints') and ctx.deps.checkpoints:
                ctx.deps.checkpoints.kallisto_index_selection.status = CheckpointStatus.FAILED
                ctx.deps.checkpoints.kallisto_index_selection.error_message = msg
                ctx.deps.checkpoints.kallisto_index_selection.timestamp = datetime.datetime.now().isoformat()
            
            raise FileNotFoundError(msg)

        # Update checkpoint: kallisto index selection completed
        if hasattr(ctx.deps, 'checkpoints') and ctx.deps.checkpoints:
            ctx.deps.checkpoints.kallisto_index_selection.status = CheckpointStatus.COMPLETED
            ctx.deps.checkpoints.kallisto_index_selection.details = f"Selected index: {index_path}"
            ctx.deps.checkpoints.kallisto_index_selection.timestamp = datetime.datetime.now().isoformat()
            logger.info("✅ Checkpoint: Kallisto index selection completed")
        
        # Update checkpoint: kallisto quantification started
        if hasattr(ctx.deps, 'checkpoints') and ctx.deps.checkpoints:
            ctx.deps.checkpoints.kallisto_quantification.status = CheckpointStatus.IN_PROGRESS
            ctx.deps.checkpoints.kallisto_quantification.timestamp = datetime.datetime.now().isoformat()
            logger.info("🧬 Checkpoint: Kallisto quantification started")

        # Rest of the Kallisto implementation remains unchanged
        output_dir = ctx.deps.output_dir
        os.makedirs(output_dir, exist_ok=True)

        # Identify paired files (R1/R2)
        r1_pat = re.compile(r'.*_(R1|1)\.fastq\.gz$')
        r2_pat = re.compile(r'.*_(R2|2)\.fastq\.gz$')
        pairs: Dict[str, Tuple[str,str]] = {}
        for f in fastq_files:
            base = os.path.basename(f)
            if r1_pat.match(base):
                mate = base.replace('_R1', '_R2').replace('_1.fastq', '_2.fastq')
                mate_path = os.path.join(os.path.dirname(f), mate)
                if mate_path in fastq_files:
                    sample = base.split('_R1')[0].split('_1.fastq')[0]
                    pairs[sample] = (f, mate_path)

        logger.info("🔍 Identified %d paired samples", len(pairs))
        if not pairs:
            msg = "Error: No paired FASTQ files found"
            logger.error("❌ %s", msg)
            raise ValidationError(msg)

        total_cpus = int(os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count() or 1))
        num_samples = len(pairs)
        desired_parallel = min(num_samples, total_cpus)
        threads_per_job = max(total_cpus // desired_parallel, 1)

        logger.info("⚙️ Running Kallisto using index: %s", index_path)
        logger.info("⚙️ Running %d Kallisto jobs in parallel, %d threads/job",
                    desired_parallel, threads_per_job)

        def run_kallisto_job(sample, r1, r2):
            sample_out = os.path.join(output_dir, "abundance", sample)
            os.makedirs(sample_out, exist_ok=True)
            cmd = [
                "kallisto", "quant", "--rf-stranded",
                "-i", index_path,
                "-o", sample_out,
                "-t", str(threads_per_job),
                "--plaintext", r1, r2
            ]
            result = subprocess.run(cmd, capture_output=True, text=True)
            return {
                "sample": sample,
                "returncode": result.returncode,
                "stdout": result.stdout,
                "stderr": result.stderr,
                "abundance_path": os.path.join(sample_out, "abundance.tsv") if result.returncode == 0 else None
            }

        # Use ThreadPoolExecutor to launch jobs in parallel
        results = []
        abundance_files = []
        with ThreadPoolExecutor(max_workers=desired_parallel) as pool:
            future_to_sample = {
                pool.submit(run_kallisto_job, sample, pair[0], pair[1]): sample
                for sample, pair in pairs.items()
            }
            for fut in as_completed(future_to_sample):
                res = fut.result()
                sample = res["sample"]
                if res["returncode"] == 0 and os.path.exists(res["abundance_path"]):
                    abundance_files.append(res["abundance_path"])
                    logger.info("✅ Kallisto completed for %s", sample)
                else:
                    logger.error("❌ Kallisto failed for %s. STDERR: %s", sample, res['stderr'])
                results.append(res)

        # Store abundance files
        ctx.deps.abundance_files = abundance_files

        # Add abundance files to context files list
        if not hasattr(ctx.deps, 'files') or ctx.deps.files is None:
            ctx.deps.files = []

        # Add abundance files to the files list, avoiding duplicates
        existing_files = set(ctx.deps.files)
        ctx.deps.files = list(existing_files.union(set(abundance_files)))
        logger.info("💾 Stored %d abundance files in context", len(abundance_files))

        # Update checkpoint: kallisto quantification completed successfully
        if hasattr(ctx.deps, 'checkpoints') and ctx.deps.checkpoints:
            ctx.deps.checkpoints.kallisto_quantification.status = CheckpointStatus.COMPLETED
            ctx.deps.checkpoints.kallisto_quantification.details = f"Generated {len(abundance_files)} abundance files"
            ctx.deps.checkpoints.kallisto_quantification.timestamp = datetime.datetime.now().isoformat()
            logger.info("✅ Checkpoint: Kallisto quantification completed successfully")

        out_lines = [f"Kallisto quant: {r['sample']} (code:{r['returncode']})" for r in results]
        return f"""
Parallel Kallisto quantification finished – {len(abundance_files)} abundance files ready.

Job info:
{chr(10).join(out_lines)}

Each job used {threads_per_job} threads; {desired_parallel} jobs run in parallel.

Files used:
- Kallisto index: {index_path}
"""
    except Exception as e:
        # Update checkpoint: kallisto quantification failed
        if hasattr(ctx.deps, 'checkpoints') and ctx.deps.checkpoints:
            ctx.deps.checkpoints.kallisto_quantification.status = CheckpointStatus.FAILED
            ctx.deps.checkpoints.kallisto_quantification.error_message = str(e)
            ctx.deps.checkpoints.kallisto_quantification.timestamp = datetime.datetime.now().isoformat()
        
        logger.exception("❌ run_kallisto_quantification crashed: %s", e)
        return f"Error running Kallisto quantification: {e}"

#############################
# Section: Differential Expression Analysis Tools
#############################

@rnaseq_agent.tool
@log_tool
@log_tool_for_reflection
async def prepare_edgeR_analysis(ctx: RunContext[AnalysisContext]) -> str:
    """
    Prepare a sample mapping table for downstream edgeR differential expression analysis
    by matching Kallisto abundance files with sample metadata.

    Inputs:
      - ctx: RunContext[AnalysisContext]
          Must include:
             • abundance_files: A list of abundance file paths from the Kallisto quantification step.
             • metadata_df: The loaded metadata DataFrame.
             • merged_column: The column name determined by previous analysis or merge steps.
             • output_dir: Where to save the prepared edgeR sample mapping CSV.

    Process:
      1. Validates that Kallisto quantification has already been run (abundance_files exist).
      2. Confirms that metadata has been loaded and a merged group column has been identified.
      3. Extracts sample names from abundance file paths.
      4. Matches these sample names to metadata rows using various strategies (exact matching,
         substring matching, reverse matching).
      5. Constructs a pandas DataFrame that maps each sample to its metadata.
      6. Saves this mapping DataFrame to a CSV file (named "edgeR_analysis_samples.csv") in the output directory.

    Output:
      Returns a detailed multiline string that includes:
         • The number of samples successfully mapped.
         • A preview of the sample-to-metadata mapping.
         • A summary of the group counts based on the merged analysis column.

    Purpose in pipeline:
      This tool bridges the quantification and differential expression steps by ensuring that each
      sample's expression data is accurately linked to its experimental metadata, a prerequisite for edgeR analysis.
      It assumes Kallisto quantification has already been completed in a previous step.
    """
    try:
        logger.info("📊 prepare_edgeR_analysis started")
        
        # Update checkpoint: edgeR/limma preparation started
        if hasattr(ctx.deps, 'checkpoints') and ctx.deps.checkpoints:
            ctx.deps.checkpoints.edger_limma_preparation.status = CheckpointStatus.IN_PROGRESS
            ctx.deps.checkpoints.edger_limma_preparation.timestamp = datetime.datetime.now().isoformat()
            logger.info("🔬 Checkpoint: edgeR/limma preparation started")

        # Check if we have abundance files
        if not ctx.deps.abundance_files:
            error_msg = "Error: No abundance files found. Please run Kallisto quantification first."
            logger.error("❌ %s", error_msg)
            raise PrerequisiteError(error_msg)

        # Check if we have metadata
        if ctx.deps.metadata_df is None:
            error_msg = "Error: Metadata not loaded. Please use the metadata agent."
            logger.error("❌ %s", error_msg)
            raise PrerequisiteError(error_msg)

        # Check if we have merged column
        if ctx.deps.merged_column is None:
            error_msg = "Error: Analysis column not identified. Please use the metadata agent."
            logger.error("❌ %s", error_msg)
            raise PrerequisiteError(error_msg)

        # Create output directory if it doesn't exist
        os.makedirs(ctx.deps.output_dir, exist_ok=True)

        # Get sample names from abundance file paths
        sample_names = [os.path.basename(os.path.dirname(f))
                        for f in ctx.deps.abundance_files]
        logger.info("🔍 Extracted %d sample names from abundance files", len(sample_names))

        # Create a mapping between abundance files and metadata
        # First, try to match based on exact sample names
        metadata_df = ctx.deps.metadata_df
        matched_samples = {}
        unmatched_samples = []

        for i, sample_name in enumerate(sample_names):
            # Try different ways to match samples
            # 1. Exact match in any column
            exact_matches = []
            for col in metadata_df.columns:
                matches = metadata_df[metadata_df[col].astype(
                    str) == sample_name].index.tolist()
                exact_matches.extend(matches)

            if exact_matches:
                matched_samples[sample_name] = {
                    'abundance_file': ctx.deps.abundance_files[i],
                    'metadata_row': exact_matches[0]
                }
            else:
                # 2. Try to find a substring match
                substring_matches = []
                for col in metadata_df.columns:
                    # Check if the sample name is contained in any column value
                    matches = metadata_df[metadata_df[col].astype(str).str.contains(
                        sample_name, case=False, na=False)].index.tolist()
                    substring_matches.extend(matches)

                if substring_matches:
                    matched_samples[sample_name] = {
                        'abundance_file': ctx.deps.abundance_files[i],
                        'metadata_row': substring_matches[0]
                    }
                else:
                    unmatched_samples.append(sample_name)

        if unmatched_samples:
            # Try reverse matching - look for metadata values in sample names
            for sample_name in unmatched_samples[:]:
                for col in metadata_df.columns:
                    for idx, value in metadata_df[col].items():
                        if str(value) in sample_name:
                            matched_samples[sample_name] = {
                                'abundance_file': ctx.deps.abundance_files[sample_names.index(sample_name)],
                                'metadata_row': idx
                            }
                            unmatched_samples.remove(sample_name)
                            break
                    if sample_name not in unmatched_samples:
                        break

        # If there are still unmatched samples, try a more aggressive approach
        if unmatched_samples and len(unmatched_samples) < len(sample_names) // 2:
            # If we have matched most samples, we can guess the rest based on order
            if len(metadata_df) == len(sample_names):
                sorted_samples = sorted(sample_names)
                sorted_metadata = metadata_df.sort_index()

                for i, sample_name in enumerate(sorted_samples):
                    if sample_name in unmatched_samples:
                        matched_samples[sample_name] = {
                            'abundance_file': ctx.deps.abundance_files[sample_names.index(sample_name)],
                            'metadata_row': sorted_metadata.index[i]
                        }
                unmatched_samples = []

        if unmatched_samples:
            error_msg = f"""
Error: Could not match {len(unmatched_samples)} of {len(sample_names)} samples to metadata.
Unmatched samples: {', '.join(unmatched_samples)}

Please check that sample names in the FASTQ files correspond to identifiers in the metadata.
            """
            logger.error("❌ %s", error_msg)
            raise ValidationError(error_msg)

        # Create a DataFrame for edgeR analysis
        analysis_df = pd.DataFrame(index=list(matched_samples.keys()))

        # Add the file paths
        analysis_df['abundance_file'] = [matched_samples[s]
                                         ['abundance_file'] for s in analysis_df.index]

        # Add the metadata
        for col in metadata_df.columns:
            analysis_df[col] = [metadata_df.loc[matched_samples[s]
                                                ['metadata_row'], col] for s in analysis_df.index]

        # Save the analysis dataframe for later use in metadata directory
        metadata_dir = os.path.join(ctx.deps.output_dir, "metadata")
        os.makedirs(metadata_dir, exist_ok=True)
        analysis_df_path = os.path.join(metadata_dir, "edger_analysis_samples.csv")
        analysis_df.to_csv(analysis_df_path)

        # Store in a registry for later reference:
        ctx.deps.file_registry = getattr(ctx.deps, 'file_registry', {})
        ctx.deps.file_registry['sample_mapping'] = analysis_df_path
        ctx.deps.sample_mapping = analysis_df_path
        logger.info("💾 Saved sample mapping to %s", analysis_df_path)

        # Update checkpoint: edgeR/limma preparation completed successfully
        if hasattr(ctx.deps, 'checkpoints') and ctx.deps.checkpoints:
            ctx.deps.checkpoints.edger_limma_preparation.status = CheckpointStatus.COMPLETED
            ctx.deps.checkpoints.edger_limma_preparation.details = f"Prepared {len(analysis_df)} samples for analysis"
            ctx.deps.checkpoints.edger_limma_preparation.timestamp = datetime.datetime.now().isoformat()
            logger.info("✅ Checkpoint: edgeR/limma preparation completed successfully")

        result = f"""
Successfully prepared data for DESeq2 analysis with {len(analysis_df)} samples. The sample mapping file can be found at {analysis_df_path}.

Sample mapping:
{analysis_df[['abundance_file', ctx.deps.merged_column]].head().to_string()}
{' ... ' if len(analysis_df) > 5 else ''}
{analysis_df[['abundance_file', ctx.deps.merged_column]].tail().to_string() if len(analysis_df) > 5 else ''}

Group counts:
{analysis_df[ctx.deps.merged_column].value_counts().to_string()}

Analysis is ready to proceed with the following groups: {', '.join(analysis_df[ctx.deps.merged_column].unique())}
        """
        logger.info("✅ Successfully prepared edgeR analysis with %d samples", len(analysis_df))
        return result
    except Exception as e:
        error_msg = f"Error preparing edgeR analysis: {str(e)}"
        logger.error("❌ %s", error_msg, exc_info=True)
        
        # Update checkpoint: edgeR/limma preparation failed
        if hasattr(ctx.deps, 'checkpoints') and ctx.deps.checkpoints:
            ctx.deps.checkpoints.edger_limma_preparation.status = CheckpointStatus.FAILED
            ctx.deps.checkpoints.edger_limma_preparation.error_message = error_msg
            ctx.deps.checkpoints.edger_limma_preparation.timestamp = datetime.datetime.now().isoformat()
        
        return error_msg


@rnaseq_agent.tool
@log_tool
@log_tool_for_reflection
async def run_edger_limma_analysis(ctx: RunContext[AnalysisContext],
    tx2gene_path: Optional[str] = None) -> str:
    """
    Run edgeR/limma analysis using the metadata agent's contrasts.

    This tool performs differential expression analysis using:
    1. The sample mapping created in a previous step
    2. The contrasts defined by the metadata agent
    3. An R script that handles the analysis

    Parameters:
    -----------
    ctx : RunContext[AnalysisContext]
        The context containing metadata, sample mapping, and contrasts information

    tx2gene_path : Optional[str]
        Path to a transcript-to-gene mapping file that maps transcript IDs to gene IDs.
        This file is essential for summarizing transcript-level Kallisto quantification
        to gene-level counts for differential expression analysis.
        Leave this option blank if you need help locating the file - once you are confident in where it is, you should specify it here.

    Finding the tx2gene file:
    ------------------------
    If tx2gene_path is not provided or invalid, you must find the appropriate file:

    1. DO NOT use the Kallisto index (.idx) file, or any other files that is NOT tx2gene file - this is incorrect and will cause errors
    2. Use the list_files tool to search in the appropriate directory with pattern "*t2g.txt"
    3. If you do not provide a tx2gene file, the error message will help you identify where specifically to look and what tool you should use.

    Analysis workflow:
    ----------------
    1. Loads quantification data from Kallisto via tximport
    2. Creates and normalizes a DGEList
    3. Builds design matrix using the grouping column
    4. Performs voom transformation and fits limma model
    5. Applies contrasts defined by the metadata agent
    6. Saves results and generates various plots

    Returns:
    -------
    str: Summary of the analysis process including contrasts used and output locations

    Note:
    -----
    If this function runs without tx2gene_path and returns an error message, examine
    the error carefully as it often contains useful information about where to find
    appropriate tx2gene files for your organism.
    """
    try:
        logger.info("📊 run_edger_limma_analysis started")
        
        # Update checkpoint: RNAseq analysis started
        if hasattr(ctx.deps, 'checkpoints') and ctx.deps.checkpoints:
            ctx.deps.checkpoints.rnaseq_analysis.status = CheckpointStatus.IN_PROGRESS
            ctx.deps.checkpoints.rnaseq_analysis.timestamp = datetime.datetime.now().isoformat()
            logger.info("🧪 Checkpoint: RNAseq analysis started")

        # Ensure the output directory exists
        os.makedirs(ctx.deps.output_dir, exist_ok=True)

        # Check if we have contrasts available
        if getattr(ctx.deps, 'contrast_matrix_df', None) is None and not getattr(ctx.deps, 'contrasts', None):
            error_msg = "No contrasts defined. Please run process_metadata_with_agent first."
            logger.error("❌ %s", error_msg)
            raise PrerequisiteError(error_msg)

        # Check if sample_mapping exists and is valid
        sample_mapping_path = getattr(ctx.deps, 'sample_mapping', None)
        if not sample_mapping_path:
            error_msg = "Error: Sample mapping not defined. Please run prepare_edgeR_analysis first to create the sample metadata mapping."
            logger.error("❌ %s", error_msg)
            raise PrerequisiteError(error_msg)

        if not os.path.exists(sample_mapping_path):
            error_msg = f"Error: Sample mapping file not found at {sample_mapping_path}. Please run prepare_edgeR_analysis to generate this file."
            logger.error("❌ %s", error_msg)
            raise FileNotFoundError(error_msg)

        # If we have contrasts but they're not in DataFrame format, convert them
        contrasts = getattr(ctx.deps, 'contrasts', None)
        contrast_matrix_df = getattr(ctx.deps, 'contrast_matrix_df', None)
        
        if contrast_matrix_df is None and contrasts:
            logger.info("🔄 Converting contrasts from agent output to DataFrame format")
            contrast_data = []
            for contrast in contrasts.data.contrasts:
                contrast_data.append({
                    'name': contrast.name,
                    'expression': contrast.expression,
                    'description': contrast.description if hasattr(contrast, 'description') else "",
                    'justification': contrast.justification if hasattr(contrast, 'justification') else ""
                })
            if contrast_data:
                setattr(ctx.deps, 'contrast_matrix_df', pd.DataFrame(contrast_data))

                # Save contrasts to a CSV file if not already done
                contrast_path = getattr(ctx.deps, 'contrast_path', None)
                if contrast_path is None:
                    metadata_dir = os.path.join(ctx.deps.output_dir, "metadata")
                    os.makedirs(metadata_dir, exist_ok=True)
                    contrast_path = os.path.join(metadata_dir, "contrasts.csv")
                    pd.DataFrame(contrast_data).to_csv(contrast_path, index=False)
                    setattr(ctx.deps, 'contrast_path', contrast_path)
                    logger.info("💾 Saved contrasts to %s", contrast_path)

        # Get the sample mapping path (redundant check, but kept for safety)
        sample_mapping_path = getattr(ctx.deps, 'sample_mapping', None)
        if not sample_mapping_path or not os.path.exists(sample_mapping_path):
            error_msg = "Error: Sample mapping file not found. Please run prepare_edgeR_analysis first."
            logger.error("❌ %s", error_msg)
            return error_msg

        # Determine the tx2gene file argument; if not provided, raise an error
        tx2gene_arg = tx2gene_path
        if not tx2gene_arg or not os.path.exists(tx2gene_arg):
            error_msg = f"Error: tx2gene file not found. Please use the list_files function, and specify the t2g.txt pattern in the {ctx.deps.resource_dir} directory to identify the possible tx2gene files. Please note the file is called EXACTLY t2g.txt, though the directory it is in will vary. Doing this step will provide you with a list of candidate files, which are applicable for different species. Ensure you select the correct file for the identified species: {ctx.deps.organism}."
            logger.error("❌ %s", error_msg)
            raise FileNotFoundError(error_msg)

        # Use the main script
        main_r_script_path = "./main_workflow/additional_scripts/RNAseq.R"

     # Execute the R script with the necessary arguments
        cmd = ['Rscript', main_r_script_path, sample_mapping_path,
               ctx.deps.merged_column, ctx.deps.output_dir, tx2gene_arg]

        # Add the contrasts CSV file as the fifth argument if it exists
        contrast_path = getattr(ctx.deps, 'contrast_path', None)
        if contrast_path and os.path.exists(contrast_path):
            cmd.append(contrast_path)

        # Check for None values in the command
        for i, arg in enumerate(cmd):
            if arg is None:
                logger.warning("⚠️ Command argument %d is None!", i)

        logger.info("⚙️ Running R script: %s", ' '.join(cmd))
        process = subprocess.run(cmd, capture_output=True, text=True)

        # Log the R script output at appropriate levels
        stdout_output = process.stdout or ""
        stderr_output = process.stderr or ""

        # Store in context for later reference by the agent
        ctx.deps.r_script_stdout = stdout_output
        ctx.deps.r_script_stderr = stderr_output

        if stdout_output:
            logger.info("📜 R script stdout:\n%s", stdout_output)

        if stderr_output:
            if process.returncode != 0:
                logger.error("❌ R script stderr:\n%s", stderr_output)
            else:
                # Some R packages output warnings to stderr even on success
                logger.warning("⚠️ R script stderr:\n%s", stderr_output)

        # Check the return code
        if process.returncode != 0:
            error_msg = f"edgeR/limma analysis failed with return code: {process.returncode}."
            detailed_error = f"{error_msg} See error output above.\n\nSTDEERR:\n{stderr_output[:1000]}...(truncated)" if len(stderr_output) > 1000 else f"{error_msg}\n\nSTDEERR:\n{stderr_output}"
            logger.error("❌ %s", error_msg)
            raise WorkflowError(detailed_error)


        # Find and record DEG result files in the new contrast-specific directories
        analysis_dir = os.path.join(ctx.deps.output_dir, "RNAseqAnalysis")
        if os.path.exists(analysis_dir):
            # Look for DEG files in each contrast subdirectory
            deg_files = []
            for subdir in os.listdir(analysis_dir):
                subdir_path = os.path.join(analysis_dir, subdir)
                if os.path.isdir(subdir_path) and subdir != "logs":
                    # Check for DEG.csv in this contrast directory
                    deg_file = os.path.join(subdir_path, "DEG.csv")
                    if os.path.exists(deg_file):
                        deg_files.append(deg_file)
            
            if deg_files:
                setattr(ctx.deps, 'deg_results_path', deg_files[0] if len(deg_files) == 1 else deg_files)
                logger.info("💾 Found %d differential expression result files in contrast directories", len(deg_files))

        # Update checkpoint: RNAseq analysis completed successfully
        if hasattr(ctx.deps, 'checkpoints') and ctx.deps.checkpoints:
            ctx.deps.checkpoints.rnaseq_analysis.status = CheckpointStatus.COMPLETED
            ctx.deps.checkpoints.rnaseq_analysis.details = "DEG analysis completed, files generated"
            ctx.deps.checkpoints.rnaseq_analysis.timestamp = datetime.datetime.now().isoformat()
            logger.info("✅ Checkpoint: RNAseq analysis completed successfully")

        logger.info("✅ edgeR/limma analysis completed successfully")

        return f"""
edgeR/limma analysis completed successfully.

Analysis performed using:
- {len(ctx.deps.contrast_matrix_df)} contrasts defined by the metadata agent
- Grouping column: {ctx.deps.merged_column}
- {len(ctx.deps.unique_groups or [])} unique groups

Results saved to {ctx.deps.output_dir}
Analysis directory: {os.path.join(ctx.deps.output_dir, "RNAseqAnalysis")}
Dataset-wide plots (MDS, CPM, etc.) are in the main analysis directory
Each contrast has its own subdirectory with DEG.csv and specific plots

==== R Script Output (truncated) ====
{stdout_output[:1500]}
{"..." if len(stdout_output) > 1500 else ""}

==== R Script Warnings/Messages (if any) ====
{stderr_output[:1000]}
{"..." if len(stderr_output) > 1000 else ""}
"""

    except Exception as e:
        error_msg = f"Error in run_edger_limma_analysis: {str(e)}"
        logger.error("❌ %s", error_msg, exc_info=True)
        
        # Update checkpoint: RNAseq analysis failed
        if hasattr(ctx.deps, 'checkpoints') and ctx.deps.checkpoints:
            ctx.deps.checkpoints.rnaseq_analysis.status = CheckpointStatus.FAILED
            ctx.deps.checkpoints.rnaseq_analysis.error_message = error_msg
            ctx.deps.checkpoints.rnaseq_analysis.timestamp = datetime.datetime.now().isoformat()
        
        return error_msg

@rnaseq_agent.tool
@log_agent_tool
@log_tool_for_reflection
async def process_metadata_with_agent(ctx: RunContext[AnalysisContext]) -> str:
    """
    Process metadata using a specialized metadata agent that has access to analysis context.
    The metadata agent will receive information about any previous analysis attempts and reflections.
    """
    try:
        logger.info("📋 process_metadata_with_agent started")
        
        # Update checkpoint: metadata analysis started
        if hasattr(ctx.deps, 'checkpoints') and ctx.deps.checkpoints:
            ctx.deps.checkpoints.metadata_analysis.status = CheckpointStatus.IN_PROGRESS
            ctx.deps.checkpoints.metadata_analysis.timestamp = datetime.datetime.now().isoformat()
            logger.info("📊 Checkpoint: Metadata analysis started")
        
        metadata_path = getattr(ctx.deps, 'metadata_path', None)
        logger.info("🔍 Processing metadata at: %s", metadata_path)

        # Check if metadata file exists
        if not metadata_path or not os.path.exists(metadata_path):
            error_msg = f"Error: Metadata file not found at {metadata_path}"
            logger.error("❌ %s", error_msg)
            raise FileNotFoundError(error_msg)

        # Import the metadata agent module
        try:
            from agents import metadata
            logger.info("✅ Successfully imported metadata module")
        except ImportError as e:
            error_msg = f"Error importing metadata_agent: {str(e)}"
            logger.error("❌ %s", error_msg, exc_info=True)
            raise WorkflowError(error_msg)

        # Create a MetadataContext instance specifically for the metadata agent
        metadata_deps = metadata.MetadataContext(metadata_path=metadata_path)

        # Generate the metadata prompt, incorporating analysis history if available
        base_prompt = """
        Please analyze the RNAseq metadata file and perform the following tasks:
        1. Process and clean the metadata
        2. Identify biologically relevant columns for analysis
        3. Create a final grouping variable (merging columns if needed)
        4. Extract the unique values found in the analysis column
        5. Design appropriate contrasts for differential expression analysis

        You should handle any errors or special cases in the data, and make appropriate decisions
        about which steps to take based on the data characteristics.
        """

        # Include dataset information context if available
        dataset_context = ""
        dataset_information = getattr(ctx.deps, 'dataset_information', None)
        if dataset_information:
            dataset_context = f"""

            IMPORTANT DATASET CONTEXT: Below is information about the dataset you are analyzing:

            {dataset_information}

            Please consider this information when selecting relevant columns and designing contrasts.
            The biological context above may help you understand which experimental factors are most important.
            """
            
        # Include reflection information if available
        reflection_context = ""
        reflections = getattr(ctx.deps, 'reflections', [])
        if reflections:
            reflection_list = "\n".join([f"Reflection {i+1}: {r}" for i, r in enumerate(reflections)])
            reflection_context = f"""

            IMPORTANT REFLECTIONS FROM PREVIOUS ATTEMPTS: 
            The analysis agent has provided the following reflections about previous attempts.
            Please pay close attention to these reflections as they identify specific issues
            that need to be addressed:

            {reflection_list}

            These reflections highlight issues that were encountered in previous analysis attempts.
            Please consider them carefully when selecting columns and designing contrasts.
            """

        # Include analysis history context if available
        analysis_context = ""
        analysis_history = getattr(ctx.deps, 'analysis_history', [])
        if analysis_history:
            # Find the most recent entry with a non-None output
            latest = None
            for entry in reversed(analysis_history):
                if entry.get("output") is not None:
                    latest = entry
                    break

            # Add analysis context to prompt if we found valid output
            if latest and latest.get("output"):
                output_text = latest.get("output", "")
                truncated_output = output_text[:2000] if output_text else ""

                analysis_context = f"""

                IMPORTANT CONTEXT: I am providing you with information about the ongoing analysis process.
                The analysis agent is currently on iteration {latest.get('iteration', '?')} and has provided the following output:

                ----- ANALYSIS AGENT OUTPUT -----
                {truncated_output}
                {"..." if output_text and len(output_text) > 2000 else ""}
                ----- END ANALYSIS AGENT OUTPUT -----

                If you see any errors or issues related to metadata processing in the analysis output,
                please keep them in mind as you analyze the metadata. Your goal is to select the most
                appropriate columns for analysis to prevent downstream errors.
                """

        # Combine all contexts into the final prompt
        metadata_prompt = base_prompt + dataset_context + reflection_context + analysis_context

        # Run the metadata agent with the Contrasts result type
        logger.info("🤖 Running metadata agent...")
        metadata_result = await metadata.run_agent_async(
            metadata_prompt,
            deps=metadata_deps,
            output_type=Contrasts
        )

        if not metadata_result or not hasattr(metadata_result, 'output'):
            error_msg = "Metadata agent returned no valid results"
            logger.error("❌ %s", error_msg)
            raise WorkflowError(error_msg)

        logger.info("✅ Metadata agent completed successfully")

        # Transfer the key information from the metadata agent back to the main agent context
        setattr(ctx.deps, 'metadata_df', metadata_deps.metadata_df)
        setattr(ctx.deps, 'merged_column', metadata_deps.merged_column)
        setattr(ctx.deps, 'unique_groups', metadata_deps.unique_groups)
        setattr(ctx.deps, 'contrasts', metadata_result)

        # Create a contrast DataFrame from the agent's output
        contrast_data = []
        for contrast in metadata_result.output.contrasts:
            contrast_data.append({
                'name': contrast.name,
                'expression': contrast.expression,
                'description': contrast.description if hasattr(contrast, 'description') else "",
                'justification': contrast.justification if hasattr(contrast, 'justification') else ""
            })

        # Convert to DataFrame and store in the context
        if contrast_data:
            setattr(ctx.deps, 'contrast_matrix_df', pd.DataFrame(contrast_data))
            logger.info("📊 Created contrast matrix with %d contrasts", len(contrast_data))

            # Save contrasts to a CSV file in metadata directory for later use
            output_dir = getattr(ctx.deps, 'output_dir', '.')
            metadata_dir = os.path.join(output_dir, "metadata")
            os.makedirs(metadata_dir, exist_ok=True)
            contrast_path = os.path.join(metadata_dir, "contrasts.csv")
            pd.DataFrame(contrast_data).to_csv(contrast_path, index=False)
            setattr(ctx.deps, 'contrast_path', contrast_path)
            logger.info("💾 Saved contrasts to %s", contrast_path)
        else:
            logger.warning("⚠️ No contrasts were generated by the metadata agent")

        # Generate a summary for the main agent
        merged_column = getattr(ctx.deps, 'merged_column', 'Not defined')
        unique_groups = getattr(ctx.deps, 'unique_groups', 'Not identified')
        
        summary = f"""
Metadata processing completed successfully.

Selected analysis column: {merged_column}
Unique groups identified: {unique_groups}

Designed contrasts:
"""
        for contrast in metadata_result.output.contrasts:
            summary += f"- {contrast.name}: {contrast.expression}\n"

        if hasattr(metadata_result.output, 'summary'):
            summary += f"\nSummary:\n{metadata_result.output.summary}\n"

        contrast_path = getattr(ctx.deps, 'contrast_path', None)
        if contrast_path:
            summary += f"\nContrasts saved to: {contrast_path}\n"

        # Update checkpoint: metadata analysis completed successfully
        if hasattr(ctx.deps, 'checkpoints') and ctx.deps.checkpoints:
            ctx.deps.checkpoints.metadata_analysis.status = CheckpointStatus.COMPLETED
            ctx.deps.checkpoints.metadata_analysis.details = f"Selected analysis column: {merged_column}"
            ctx.deps.checkpoints.metadata_analysis.timestamp = datetime.datetime.now().isoformat()
            logger.info("✅ Checkpoint: Metadata analysis completed successfully")

        logger.info("✅ Metadata processing pipeline completed")
        return summary

    except Exception as e:
        error_msg = f"Error processing metadata with agent: {str(e)}"
        logger.error("❌ %s", error_msg, exc_info=True)
        
        # Update checkpoint: metadata analysis failed
        if hasattr(ctx.deps, 'checkpoints') and ctx.deps.checkpoints:
            ctx.deps.checkpoints.metadata_analysis.status = CheckpointStatus.FAILED
            ctx.deps.checkpoints.metadata_analysis.error_message = error_msg
            ctx.deps.checkpoints.metadata_analysis.timestamp = datetime.datetime.now().isoformat()
        
        return error_msg

@log_agent_tool
async def run_agent_async(prompt: str, deps: AnalysisContext, usage=None):
    logger.info("🚀 Analysis agent invoked with prompt: %s", prompt)
    
    # Ensure all required attributes exist
    for attr_name, default_value in [
        ('tool_logs', []),
        ('reflections', []),
        ('analysis_history', []),
        ('kallisto_index_used', None),
        ('tx2gene_file_used', None)
    ]:
        if not hasattr(deps, attr_name):
            setattr(deps, attr_name, default_value)
    
    # Save the starting length of tool logs for later reporting
    tool_logs_start_len = len(deps.tool_logs)
    
    # Run the agent
    result = await rnaseq_agent.run(prompt, deps=deps, usage=usage)

    # Log the agent's output
    logger.info("📄 Analysis agent output: %s", result.output)
    
    # Report on tool logs that were added during this run
    new_logs_count = len(deps.tool_logs) - tool_logs_start_len
    if new_logs_count > 0:
        logger.info("🔧 Collected %d tool logs during analysis run", new_logs_count)

    # Log usage statistics if available
    if hasattr(result, 'usage') and result.usage:
        try:
            usage_stats = result.usage()
            logger.info("📊 Analysis agent usage: %s", usage_stats)
        except Exception as e:
            logger.debug("Could not get usage stats: %s", e)

    # Write tool logs to a file if we have an output directory
    if hasattr(deps, 'output_dir') and deps.output_dir:
        try:
            log_dir = pathlib.Path(deps.output_dir) / "logs"
            log_dir.mkdir(parents=True, exist_ok=True)
            
            # Create a timestamped log file to preserve history
            timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            tool_log_path = log_dir / f"analysis_tool_logs_{timestamp}.json"
            
            # Also update the main log file
            main_log_path = log_dir / "analysis_tool_logs_latest.json"
            
            # Save to timestamped file
            with open(tool_log_path, 'w') as f:
                json.dump(deps.tool_logs, f, indent=2)
            
            # Update the latest log file
            with open(main_log_path, 'w') as f:
                json.dump(deps.tool_logs, f, indent=2)
                
            logger.info("💾 Saved %d tool logs to %s and %s", 
                       len(deps.tool_logs), tool_log_path.name, main_log_path.name)
        except Exception as e:
            logger.warning("⚠️ Failed to save tool logs: %s", str(e))

    return result

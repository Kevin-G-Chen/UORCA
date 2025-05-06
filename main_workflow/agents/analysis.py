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
from rich.console import Console
from pydantic import BaseModel, Field, ConfigDict
from pydantic_ai import Agent, RunContext
from shared import AnalysisContext
from shared.workflow_logging import log_tool
from agents.metadata import metadata_agent, MetadataContext
import gseapy
from unidecode import unidecode
from openai import OpenAI
import matplotlib.pyplot as plt
import nest_asyncio
from concurrent.futures import ThreadPoolExecutor, as_completed


# ── configure python‑logging ────────────────────────────────────────────────
logging.basicConfig(
    format="%(asctime)s  %(levelname)-8s  %(name)s ▶  %(message)s",
    level=logging.INFO,
    datefmt="%H:%M:%S"
)
logger = logging.getLogger(__name__)


#############################
# SECTION: Data classes and env
#############################

load_dotenv()
openai_api_key = os.getenv("OPENAI_API_KEY")
client = OpenAI()

#############################
# SECTION: Pydantic schemas
#############################
class Contrast_format(BaseModel):
    name: str
    expression: str
    description: Optional[str] = Field(description="Biological interpretation of the contrast")
    justification: Optional[str] = Field(description="Justification for the contrast design, in terms of value to the scientific community")
class Contrasts(BaseModel):
    contrasts: List[Contrast_format]
    model_config = ConfigDict(extra="allow")

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
async def list_fastq_files(ctx: RunContext[AnalysisContext]) -> str:
    """
    List all FASTQ files in the fastq_dir directory from the context.
    This tool automatically gets the fastq_dir from the context and searches for fastq.gz files.
    """
    """
    List all FASTQ files in the fastq_dir directory from the context.
    This tool automatically gets the fastq_dir from the context and searches for fastq.gz files.
    """
    logger.info("🔍 list_fastq_files started - searching in %s", ctx.deps.fastq_dir)
    fastq_dir = ctx.deps.fastq_dir

    # Check if directory exists
    if not os.path.exists(fastq_dir):
        error_msg = f"Error: Directory '{fastq_dir}' does not exist"
        logger.error("❌ %s", error_msg)
        return error_msg

    # Use find_files to find fastq.gz files
    fastq_files = await find_files(ctx, fastq_dir, 'fastq.gz')

    if not fastq_files:
        error_msg = f"No fastq.gz files found in {fastq_dir}. Directory contents: {os.listdir(fastq_dir) if os.path.isdir(fastq_dir) else 'Not a directory'}"
        logger.warning("⚠️ %s", error_msg)
        return error_msg

    result = f"Found {len(fastq_files)} fastq.gz files in {fastq_dir}"
    if len(fastq_files) <= 5:
        result += f": {', '.join(fastq_files)}"
    else:
        result += f" (first 5): {', '.join(fastq_files[:5])}..."
    logger.info("✅ %s", result)
    return result


@rnaseq_agent.tool
@log_tool
async def find_files(ctx: RunContext[AnalysisContext], directory: str, suffix: Union[str, List[str]]) -> List[str]:
    """
    Recursively search for and return a sorted list of files within the specified directory that have the given suffix.

    Inputs:
      - ctx: RunContext[AnalysisContext]
          Contains the dependency context (AnalysisContext) that holds directory information and other runtime parameters.
      - directory (str):
          The root directory path in which to search for files. This can be either an absolute or relative path.
      - suffix (Union[str, List[str]]):
          The file suffix (for example "fastq.gz" for FASTQ files) or a list of suffixes that will be used to filter matching files.

    Process:
      1. Logs the current context details (only once per run, thanks to a flag in ctx.deps).
      2. Uses os.walk to recursively traverse the directory structure and check every file name.
      3. For each file that ends with the specified suffix, concatenates its full path and adds it to a list.
      4. Returns the sorted list of matching file paths.
      5. Reports progress by logging the number of files found.

    Output:
      A sorted list of absolute file path strings that match the file suffix provided. In case of errors
      (e.g. directory not found), an error message string is returned inside a list.

    Purpose in pipeline:
      This tool locates critical input FASTQ files using the fastq.gz suffix (or other types based on suffix) from the file system,
      enabling subsequent steps (such as quantification with Kallisto) to process the correct data.
    """
    try:
        logger.info("🔍 find_files started - searching for %s in %s", suffix, directory)

        # Execute the actual file search
        matched_files = []
        for root, _, files in os.walk(directory):
            for f in files:
                if isinstance(suffix, str):
                    condition = f.endswith(suffix)
                else:
                    condition = any(f.endswith(s) for s in suffix)
                if condition:
                    matched_files.append(os.path.join(root, f))

        # Sort the files and prepare result
        matched_files = sorted(matched_files)

        # Log the result
        if matched_files:
            if len(matched_files) <= 5:
                logger.info("✅ Found %d files matching suffix '%s'. Files: %s",
                           len(matched_files), suffix, matched_files)
            else:
                logger.info("✅ Found %d files matching suffix '%s'. First 5: %s...",
                           len(matched_files), suffix, matched_files[:5])
        else:
            logger.warning("⚠️ No files matching suffix '%s' found in %s", suffix, directory)

        return matched_files

    except FileNotFoundError:
        error_msg = f"Error: Directory '{directory}' not found."
        logger.error("❌ %s", error_msg)
        return [error_msg]

    except Exception as e:
        error_msg = f"Error in find_files: {str(e)}"
        logger.error("❌ %s", error_msg, exc_info=True)
        return [error_msg]


#############################
# Section: Kallisto Quantification Tools
#############################


@rnaseq_agent.tool
@log_tool
async def find_kallisto_index(ctx: RunContext[AnalysisContext]) -> str:
    """
    Search for and return the file path of an appropriate Kallisto index based on the organism specified in AnalysisContext.

    Inputs:
      - ctx: RunContext[AnalysisContext]
          Contains:
            • kallisto_index_dir: The directory where Kallisto index (.idx) files are stored.
            • organism: The organism name (e.g., "human") which the tool uses to filter the available index files.

    Process:
      1. Logs the organism being used for the search.
      2. Calls the find_files utility to locate files with the suffix ".idx" in the specified kallisto_index_dir.
      3. Filters the list of index files to find those whose parent directory name matches the organism (case-insensitive).
      4. If a matching index is found, selects it and further checks for a transcript-to-gene mapping file (.txt) in the same directory;
         if found, sets ctx.deps.tx2gene_path accordingly.
      5. If no organism-specific index is found, defaults to the first available index.

    Output:
      Returns a string message indicating the path to the found Kallisto index and, if applicable, that the transcript-to-gene mapping file was also set.

    Purpose in pipeline:
      Finding the correct Kallisto index is essential for running quantification of RNAseq data, ensuring that the
      quantification step uses the appropriate reference for the organism under study.
    """
    try:
        organism = ctx.deps.organism.lower()
        index_dir = ctx.deps.kallisto_index_dir

        logger.info("🔍 Searching for %s Kallisto index in %s", organism, index_dir)

        # Look for index files in the specified directory
        index_files = await find_files(ctx, ctx.deps.kallisto_index_dir, '.idx')

        if not index_files:
            error_msg = f"Error: No Kallisto index files found in {index_dir}"
            logger.error("❌ %s", error_msg)
            return error_msg

        # Try to find an index matching the organism
        matching_indices = [idx for idx in index_files if organism in os.path.basename(
            os.path.dirname(idx)).lower()]

        if matching_indices:
            index_path = matching_indices[0]
            logger.info("✅ Found organism-specific Kallisto index: %s", index_path)

            # Also find the transcript-to-gene mapping file if available
            tx2gene_files = await find_files(ctx, os.path.dirname(index_path), '.txt')
            if tx2gene_files:
                t2g_files = [f for f in tx2gene_files if any(
                    x in os.path.basename(f).lower() for x in ['t2g', 'tx2gene'])]
                if t2g_files:
                    ctx.deps.tx2gene_path = t2g_files[0]
                    logger.info("📄 Found transcript-to-gene mapping file: %s", ctx.deps.tx2gene_path)

            return f"Found Kallisto index for {organism}: {index_path}"
        else:
            # If no organism-specific index found, return the first one
            fallback_index = index_files[0]
            logger.warning("⚠️ No organism-specific index for %s found, using fallback: %s", organism, fallback_index)
            return f"No index specific to {organism} found. Using the first available index: {fallback_index}"

    except Exception as e:
        error_msg = f"Error finding Kallisto index: {str(e)}"
        logger.error("❌ %s", error_msg, exc_info=True)
        return error_msg


@rnaseq_agent.tool
@log_tool
async def run_kallisto_quantification(ctx: RunContext[AnalysisContext]) -> str:
    """
    Run Kallisto quantification on paired-end FASTQ files and record the resulting abundance files.

    Inputs:
      - ctx: RunContext[AnalysisContext]
          Must include:
             • fastq_dir: The directory where FASTQ files are stored.
             • kallisto_index_dir: The directory containing the Kallisto index files.
             • output_dir: The directory to store Kallisto outputs.
             • organism: Organism information to select the correct index.

    Process:
      1. Logs current FASTQ directory and full context details.
      2. Uses find_files to search for FASTQ files with the suffix "fastq.gz" in the designated fastq_dir.
      3. If no FASTQ files are found, returns an error message.
      4. Calls find_kallisto_index to determine the correct index (and transcript-to-gene mapping if available).
      5. Extracts the index path from the returned result.
      6. Creates the output directory if it does not exist.
      7. Searches for paired FASTQ files using regular expression patterns for R1 and R2.
      8. Matches R1 with R2 files to form sample pairs.
      9. For each paired sample, builds the Kallisto command and executes it via subprocess.
      10. Logs progress before and after running quantification for each sample.
      11. Collects the paths to generated abundance files (e.g., abundance.tsv) and stores them in ctx.deps.abundance_files.

    Output:
      A string summary reporting:
         • The number of sample pairs processed.
         • The result (success or error) for each sample.
         • The total number of abundance files found for downstream analysis.

    Purpose in pipeline:
      This tool integrates the quantification step of the pipeline by using Kallisto to convert raw FASTQ reads into
      transcript abundance estimates, which are later used for differential expression analysis.
    """
    try:
        logger.info("🧬 run_kallisto_quantification started – fastq_dir=%s", ctx.deps.fastq_dir)

        fastq_files = await find_files(ctx, ctx.deps.fastq_dir, 'fastq.gz')
        if not fastq_files:
            msg = f"Error: No FASTQ files found in {ctx.deps.fastq_dir}"
            logger.error("❌ %s", msg)
            return msg

        index_result = await find_kallisto_index(ctx)
        if "Error" in index_result:
            logger.error("❌ %s", index_result)
            return index_result

        # Extract .idx path from index_result
        index_path = None
        for tok in index_result.split():
            if tok.endswith('.idx'):
                index_path = tok
                break
        if not index_path or not os.path.exists(index_path):
            msg = "Error: Could not determine Kallisto index path"
            logger.error("❌ %s", msg)
            return msg

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
            return msg

        total_cpus = int(os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count() or 1))
        num_samples = len(pairs)
        desired_parallel = min(num_samples, total_cpus)
        threads_per_job = max(total_cpus // desired_parallel, 1)

        logger.info("⚙️ Running Kallisto for %d samples using %d jobs in parallel, %d threads/job",
                    num_samples, desired_parallel, threads_per_job)

        def run_kallisto_job(sample, r1, r2):
            sample_out = os.path.join(output_dir, sample)
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

        ctx.deps.abundance_files = abundance_files
        logger.info("💾 Stored %d abundance files in ctx.deps", len(abundance_files))
        out_lines = [f"Kallisto quant: {r['sample']} (code:{r['returncode']})" for r in results]
        return f"""
Parallel Kallisto quantification finished – {len(abundance_files)} abundance files ready.

Job info:
{chr(10).join(out_lines)}

Each job used {threads_per_job} threads; {desired_parallel} jobs run in parallel.
"""
    except Exception as e:
        logger.exception("❌ run_kallisto_quantification crashed: %s", e)
        return f"Error running Kallisto quantification: {e}"

#############################
# Section: Differential Expression Analysis Tools
#############################


@rnaseq_agent.tool
@log_tool
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

        # Check if we have abundance files
        if not ctx.deps.abundance_files:
            error_msg = "Error: No abundance files found. Please run Kallisto quantification first."
            logger.error("❌ %s", error_msg)
            return error_msg

        # Check if we have metadata
        if ctx.deps.metadata_df is None:
            error_msg = "Error: Metadata not loaded. Please run load_metadata first."
            logger.error("❌ %s", error_msg)
            return error_msg

        # Check if we have merged column
        if ctx.deps.merged_column is None:
            error_msg = "Error: Analysis column not identified. Please run identify_analysis_columns first."
            logger.error("❌ %s", error_msg)
            return error_msg

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
            warning_msg = f"""
Warning: Could not match {len(unmatched_samples)} of {len(sample_names)} samples to metadata.
Unmatched samples: {', '.join(unmatched_samples)}

Please check that sample names in the FASTQ files correspond to identifiers in the metadata.
            """
            logger.warning("⚠️ %s", warning_msg)
            return warning_msg

        # Create a DataFrame for edgeR analysis
        analysis_df = pd.DataFrame(index=list(matched_samples.keys()))

        # Add the file paths
        analysis_df['abundance_file'] = [matched_samples[s]
                                         ['abundance_file'] for s in analysis_df.index]

        # Add the metadata
        for col in metadata_df.columns:
            analysis_df[col] = [metadata_df.loc[matched_samples[s]
                                                ['metadata_row'], col] for s in analysis_df.index]

        # Save the analysis dataframe for later use
        analysis_df_path = os.path.join(
            ctx.deps.output_dir, "edger_analysis_samples.csv")
        analysis_df.to_csv(analysis_df_path)

        # Store in a registry for later reference:
        ctx.deps.file_registry = getattr(ctx.deps, 'file_registry', {})
        ctx.deps.file_registry['sample_mapping'] = analysis_df_path
        ctx.deps.sample_mapping = analysis_df_path
        logger.info("💾 Saved sample mapping to %s", analysis_df_path)

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
        return error_msg


@rnaseq_agent.tool
@log_tool
async def run_edger_limma_analysis(ctx: RunContext[AnalysisContext]) -> str:
    """
    Run edgeR/limma analysis using the metadata agent's contrasts.

    This tool performs differential expression analysis using:
    1. The sample mapping created in a previous step
    2. The contrasts defined by the metadata agent
    3. An R script that handles the analysis

    The analysis workflow includes:
    - Loading quantification data via tximport
    - Creating and normalizing a DGEList
    - Building a design matrix using the grouping column
    - Performing voom transformation and fitting a limma model
    - Applying the contrasts defined by the metadata agent
    - Saving results and generating plots

    Results are saved in the output directory specified in the context.
    """
    try:
        logger.info("📊 run_edger_limma_analysis started")

        # Ensure the output directory exists
        os.makedirs(ctx.deps.output_dir, exist_ok=True)

        # Check if we have contrasts available
        if ctx.deps.contrast_matrix_df is None and not ctx.deps.contrasts:
            error_msg = "No contrasts defined. Please run process_metadata_with_agent first."
            logger.warning("⚠️ %s", error_msg)
            return error_msg

        # If we have contrasts but they're not in DataFrame format, convert them
        if ctx.deps.contrast_matrix_df is None and ctx.deps.contrasts:
            logger.info("🔄 Converting contrasts from agent output to DataFrame format")
            contrast_data = []
            for contrast in ctx.deps.contrasts.data.contrasts:
                contrast_data.append({
                    'name': contrast.name,
                    'expression': contrast.expression,
                    'description': contrast.description if hasattr(contrast, 'description') else "",
                    'justification': contrast.justification if hasattr(contrast, 'justification') else ""
                })
            if contrast_data:
                ctx.deps.contrast_matrix_df = pd.DataFrame(contrast_data)

                # Save contrasts to a CSV file if not already done
                if ctx.deps.contrast_path is None:
                    contrast_path = os.path.join(ctx.deps.output_dir, "contrasts.csv")
                    pd.DataFrame(contrast_data).to_csv(contrast_path, index=False)
                    ctx.deps.contrast_path = contrast_path
                    logger.info("💾 Saved contrasts to %s", contrast_path)

        # Get the sample mapping path
        sample_mapping_path = ctx.deps.sample_mapping

        # Determine the tx2gene file argument; if not provided, pass "NA"
        tx2gene_arg = ctx.deps.tx2gene_path if ctx.deps.tx2gene_path and os.path.exists(
            ctx.deps.tx2gene_path) else "NA"

        # Use the main script
        main_r_script_path = "./main_workflow/additional_scripts/RNAseq.R"

        # Check if the R script exists
        if not os.path.exists(main_r_script_path):
            # Try an alternative path if the first one doesn't exist
            main_r_script_path = "aim1/UORCA/script_development/experiments/sample_RNAseq.R"
            if not os.path.exists(main_r_script_path):
                error_msg = "Error: R script not found at expected paths. Please verify the R script location."
                logger.error("❌ %s", error_msg)
                return error_msg

        # Execute the R script with the necessary arguments
        cmd = ['Rscript', main_r_script_path, sample_mapping_path,
               ctx.deps.merged_column, ctx.deps.output_dir, tx2gene_arg]

        # Add the contrasts CSV file as the fifth argument if it exists
        if ctx.deps.contrast_path and os.path.exists(ctx.deps.contrast_path):
            cmd.append(ctx.deps.contrast_path)

        # Check for None values in the command
        for i, arg in enumerate(cmd):
            if arg is None:
                logger.warning("⚠️ Command argument %d is None!", i)

        logger.info("⚙️ Running R script: %s", ' '.join(cmd))
        process = subprocess.run(cmd, capture_output=True, text=True)

        # Log the R script output at appropriate levels
        if process.stdout:
            logger.info("📜 R script stdout:\n%s", process.stdout)

        if process.stderr:
            if process.returncode != 0:
                logger.error("❌ R script stderr:\n%s", process.stderr)
            else:
                # Some R packages output warnings to stderr even on success
                logger.warning("⚠️ R script stderr:\n%s", process.stderr)

        # Check the return code
        if process.returncode != 0:
            error_msg = f"edgeR/limma analysis failed with return code: {process.returncode}."
            logger.error("❌ %s", error_msg)
            return f"{error_msg} See error output above."

        # Find and record DEG result files
        deg_dir = os.path.join(ctx.deps.output_dir, "DEG")
        if os.path.exists(deg_dir):
            deg_files = [os.path.join(deg_dir, f) for f in os.listdir(deg_dir) if f.endswith('.csv')]
            if deg_files:
                ctx.deps.deg_results_path = deg_files[0] if len(deg_files) == 1 else deg_files
                logger.info("💾 Found %d differential expression result files", len(deg_files))

        logger.info("✅ edgeR/limma analysis completed successfully")

        return f"""
edgeR/limma analysis completed successfully.

Analysis performed using:
- {len(ctx.deps.contrast_matrix_df)} contrasts defined by the metadata agent
- Grouping column: {ctx.deps.merged_column}
- {len(ctx.deps.unique_groups or [])} unique groups

Results saved to {ctx.deps.output_dir}
Differential expression results saved to {deg_dir}
"""

    except Exception as e:
        error_msg = f"Error in run_edger_limma_analysis: {str(e)}"
        logger.error("❌ %s", error_msg, exc_info=True)
        return error_msg

@rnaseq_agent.tool
@log_tool
async def process_metadata_with_agent(ctx: RunContext[AnalysisContext]) -> str:
    """
    Process metadata using a specialized metadata agent.

    This tool creates a dedicated MetadataAgent to analyze the metadata file,
    identify biologically relevant columns, merge columns if needed, extract unique values,
    and generate appropriate contrasts for differential expression analysis.

    The results are stored back in the main RNAseq context for downstream analysis.
    The contrasts are also saved to a CSV file and converted to a DataFrame for use in edgeR/limma analysis.
    """
    try:
        logger.info("📋 process_metadata_with_agent started")
        logger.info("🔍 Processing metadata at: %s", ctx.deps.metadata_path)

        # Check if metadata file exists
        if not ctx.deps.metadata_path or not os.path.exists(ctx.deps.metadata_path):
            error_msg = f"Error: Metadata file not found at {ctx.deps.metadata_path}"
            logger.error("❌ %s", error_msg)
            return error_msg

        # Import the metadata agent module
        try:
            from agents.metadata import metadata_agent, MetadataContext
            logger.info("✅ Successfully imported metadata_agent module")
        except ImportError as e:
            error_msg = f"Error importing metadata_agent: {str(e)}"
            logger.error("❌ %s", error_msg, exc_info=True)
            return error_msg

        # Create a MetadataContext instance specifically for the metadata agent
        metadata_deps = MetadataContext(metadata_path=ctx.deps.metadata_path)

        # Prompt for the metadata agent
        metadata_prompt = """
        Please analyze the RNAseq metadata file and perform the following tasks:
        1. Process and clean the metadata
        2. Identify biologically relevant columns for analysis
        3. Create a final grouping variable (merging columns if needed)
        4. Extract the unique values found in the analysis column
        5. Design appropriate contrasts for differential expression analysis based on these unique values

        You should handle any errors or special cases in the data, and make appropriate decisions
        about which steps to take based on the data characteristics.
        """

        # Run the metadata agent with the Contrasts result type
        logger.info("🤖 Running metadata agent...")
        from agents.metadata import Contrasts
        metadata_result = await metadata_agent.run(
            metadata_prompt,
            deps=metadata_deps,
            output_type=Contrasts
        )

        if not metadata_result or not hasattr(metadata_result, 'output'):
            error_msg = "Metadata agent returned no valid results"
            logger.error("❌ %s", error_msg)
            return error_msg

        logger.info("✅ Metadata agent completed successfully")

        # Transfer the key information from the metadata agent back to the main agent context
        ctx.deps.metadata_df = metadata_deps.metadata_df
        ctx.deps.merged_column = metadata_deps.merged_column
        ctx.deps.unique_groups = metadata_deps.unique_groups
        ctx.deps.contrasts = metadata_result

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
            ctx.deps.contrast_matrix_df = pd.DataFrame(contrast_data)
            logger.info("📊 Created contrast matrix with %d contrasts", len(contrast_data))

            # Save contrasts to a CSV file for later use
            os.makedirs(ctx.deps.output_dir, exist_ok=True)
            contrast_path = os.path.join(ctx.deps.output_dir, "contrasts.csv")
            pd.DataFrame(contrast_data).to_csv(contrast_path, index=False)
            ctx.deps.contrast_path = contrast_path
            logger.info("💾 Saved contrasts to %s", contrast_path)
        else:
            logger.warning("⚠️ No contrasts were generated by the metadata agent")

        # Generate a summary for the main agent
        summary = f"""
Metadata processing completed successfully.

Selected analysis column: {ctx.deps.merged_column}
Unique groups identified: {ctx.deps.unique_groups}

Designed contrasts:
"""
        for contrast in metadata_result.output.contrasts:
            summary += f"- {contrast.name}: {contrast.expression}\n"

        if hasattr(metadata_result.output, 'summary'):
            summary += f"\nSummary:\n{metadata_result.output.summary}\n"

        if ctx.deps.contrast_path:
            summary += f"\nContrasts saved to: {ctx.deps.contrast_path}\n"

        logger.info("✅ Metadata processing pipeline completed")
        return summary

    except Exception as e:
        error_msg = f"Error processing metadata with agent: {str(e)}"
        logger.error("❌ %s", error_msg, exc_info=True)
        return error_msg

@log_tool
async def run_agent_async(prompt: str, deps: AnalysisContext, usage=None):
    logger.info("🚀 Analysis agent invoked by Master Agent – prompt: %s", prompt)
    return await rnaseq_agent.run(prompt, deps=deps, usage=usage)

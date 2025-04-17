#############################
# SECTION: Imports and Setup
#############################
from typing import List, Optional, Dict, Any, Union
from dataclasses import dataclass, field
from pydantic import BaseModel, Field, ConfigDict
from pydantic_ai import Agent, RunContext
import argparse, os, re, json, subprocess
import pandas as pd, numpy as np
from dotenv import load_dotenv
from rich.console import Console
import nest_asyncio
from unidecode import unidecode
from openai import OpenAI

nest_asyncio.apply()
console = Console()
load_dotenv()
openai_api_key = os.getenv("OPENAI_API_KEY")
client = OpenAI()

# Simplified logging
class LogLevel:
    MINIMAL, NORMAL, VERBOSE, DEBUG = range(4)
CURRENT_LOG_LEVEL = LogLevel.NORMAL
def log(msg, level=LogLevel.NORMAL, style=""):
    if CURRENT_LOG_LEVEL >= level:
        console.print(msg, style=style or None)

def log_tool_header(name, params=None):
    log(f"TOOL: {name} {params or ''}", LogLevel.NORMAL, "bold blue")

def log_tool_result(res):
    log(res, LogLevel.NORMAL)

#############################
# Section: Utility Functions
#############################

class Contrast_format(BaseModel):
    """Schema for each contrast."""
    name: str
    expression: str
    description: Optional[str] = Field(
        default="",
        description="Biological interpretation of the contrast"
    )
    justification: Optional[str] = Field(
        default="",
        description="Justification for the contrast design"
    )

class Contrasts(BaseModel):
    """Schema for candidate contrasts."""
    contrasts: List[Contrast_format]
    summary: Optional[str] = Field(
        default="",
        description="Summary of the designed contrasts and reflection log"
    )
    model_config = ConfigDict(extra="allow")

# And add a matching ReflectionResult so the reflection_agent can type‑check:
class ReflectionResult(BaseModel):
    contrasts: List[Contrast_format]
    summary: Optional[str] = Field(default="", description="Reflection process summary")
    model_config = ConfigDict(extra="allow")


@dataclass
class RNAseqData:
    """Container for RNAseq analysis data and paths."""
    # Input data
    fastq_dir: str
    metadata_path: str
    kallisto_index_dir: str
    organism: str = "human"  # default to human
    output_dir: str = "output"
    kallisto_index_path: Optional[str] = None
    tx2gene_path: Optional[str] = None

    # Runtime data that gets populated during analysis
    metadata_df: Optional[pd.DataFrame] = None
    abundance_files: List[str] = None
    merged_column: Optional[str] = None
    contrast_groups: Dict[str, Dict[str, str]] = None
    sample_mapping: Optional[pd.DataFrame] = None
    deg_results_path: Optional[str] = None
    contrasts: Optional[Contrasts] = None

    # Additional fields for contrast handling
    contrast_path: Optional[str] = None  # Path to the JSON file containing contrasts
    contrast_matrix_df: Optional[pd.DataFrame] = None  # DataFrame representation of contrasts

    # Additional fields for extended functionality
    unique_groups: Optional[List[str]] = None
    file_registry: Dict[str, str] = Field(default_factory=dict)

# Class to prepare contrasts is defined above - no need to redefine


# Load environment variables from .env file
load_dotenv()

# Get the API key
openai_api_key = os.getenv("OPENAI_API_KEY")
# Define Pydantic class for the structured output
client = OpenAI()

# ----------------------------
# Create an RNAseq analysis agent
# ----------------------------

# Adding temporary rough information for the agent until I develop a proper way of implementing documentation

# Read system prompt from file
data_analysis_agent_system_prompt_path = "../script_development/prompt_development/data_analysis_agent_systemprompt.txt"
metadata_prompt_path = "../script_development/prompt_development/metadata_processing_systemprompt.txt"

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
    'openai:o4-mini',  # Change to more powerful model
    deps_type=RNAseqData,
    system_prompt=system_prompt
)



# Immediately after you create metadata_agent, add:

# Load combined reflection prompt (reuse or fallback as you had)
def load_reflection_prompt():
    metadata_prompt = system_prompt  # Use the already loaded metadata agent prompt
    reflection_prompt_path = "../script_development/prompt_development/reflection_agent_system_prompt.txt"

    try:
        with open(reflection_prompt_path, 'r') as f:
            reflection_prompt = f.read()

        # Print the first and last few lines for verification
        lines = reflection_prompt.split('\n')
        print("\n--- REFLECTION PROMPT FILE LOADED ---")
        print("First 3 lines:")
        for line in lines[:3]:
            print(f"> {line}")

        print("...")

        print("Last 3 lines:")
        for line in lines[-3:]:
            print(f"> {line}")
        print(f"--- END OF PREVIEW ({len(lines)} total lines) ---\n")

        # Combine the prompts
        combined_prompt = f"""
        #### Combined Metadata Processing and Reflection Agent System Prompt

        # Part 1: Metadata Processing Guidelines
        {metadata_prompt}

        # Part 2: Reflection-Specific Guidelines
        {reflection_prompt}
        """
        return combined_prompt

    except Exception:
        fallback_prompt = """
        You are a reflection agent designed to analyze and suggest additional contrasts for RNAseq metadata analysis. You will receive a list of existing contrasts and the unique values in the analysis column. Your task is to reflect on these contrasts and suggest additional biologically meaningful contrasts that might have been missed. Consider the pairwise comparisons that could be made and any potentially meaningful group combinations. If you believe the current set of contrasts is comprehensive, please explain why no additional contrasts are needed. Do not go out of your way to try and "force" new contrasts - only suggest those that are truly meaningful and relevant to the analysis. Your output should include the new contrasts and a summary of your reflection process.
        """
        print("Using fallback reflection prompt instead.")
        return fallback_prompt

# Load the combined prompt
reflection_system_prompt = load_reflection_prompt()


reflection_agent = Agent(
    'openai:o4-mini',
    deps_type=RNAseqData,
    system_prompt=reflection_system_prompt
)
# ----------------------------
# Utility Functions
# ----------------------------


@rnaseq_agent.tool
async def list_fastq_files(ctx: RunContext[RNAseqData]) -> str:
    """
    List all FASTQ files in the fastq_dir directory from the context.
    This tool automatically gets the fastq_dir from the context and searches for fastq.gz files.
    """
    """
    List all FASTQ files in the fastq_dir directory from the context.
    This tool automatically gets the fastq_dir from the context and searches for fastq.gz files.
    """
    log_tool_header("list_fastq_files")
    fastq_dir = ctx.deps.fastq_dir

    # Check if directory exists
    if not os.path.exists(fastq_dir):
        error_msg = f"Error: Directory '{fastq_dir}' does not exist"
        log_tool_result(error_msg)
        return error_msg

    # Use find_files to find fastq.gz files
    fastq_files = await find_files(ctx, fastq_dir, 'fastq.gz')

    if not fastq_files:
        error_msg = f"No fastq.gz files found in {fastq_dir}. Directory contents: {os.listdir(fastq_dir) if os.path.isdir(fastq_dir) else 'Not a directory'}"
        log_tool_result(error_msg)
        return error_msg

    result = f"Found {len(fastq_files)} fastq.gz files in {fastq_dir}"
    if CURRENT_LOG_LEVEL >= LogLevel.VERBOSE:
        result += f": {', '.join(fastq_files)}"
    log_tool_result(result)
    return result


@rnaseq_agent.tool
async def find_files(ctx: RunContext[RNAseqData], directory: str, suffix: Union[str, List[str]]) -> List[str]:
    """
    Recursively search for and return a sorted list of files within the specified directory that have the given suffix.

    Inputs:
      - ctx: RunContext[RNAseqData]
          Contains the dependency context (RNAseqData) that holds directory information and other runtime parameters.
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
        # Log the initial context only once per run
        if not hasattr(ctx.deps, '_logged_context'):
            log(f"Initial Context.deps details:\n{vars(ctx.deps)}",
                level=LogLevel.VERBOSE, style="bold blue")
            setattr(ctx.deps, '_logged_context', True)

        # Log tool call with parameters
        log_tool_header(
            "find_files", {"directory": directory, "suffix": suffix})

        # Log context details only in verbose mode
        log(f"Context.deps details:\n{vars(ctx.deps)}",
            level=LogLevel.VERBOSE, style="bold blue")
        if hasattr(ctx, "message_history"):
            log(f"Message History: {ctx.message_history}",
                level=LogLevel.DEBUG, style="bold magenta")

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
        result_msg = f"Found {len(matched_files)} files matching suffix '{suffix}' in directory: {directory}"
        if matched_files and CURRENT_LOG_LEVEL >= LogLevel.VERBOSE:
            result_msg += f"\nFirst few files: {matched_files[:3]}"
            if len(matched_files) > 3:
                result_msg += f"\n... and {len(matched_files) - 3} more"

        log_tool_result(result_msg)
        return matched_files

    except FileNotFoundError:
        error_msg = f"Error: Directory '{directory}' not found."
        log(error_msg, style="bold red")
        log_tool_result(error_msg)
        return [error_msg]

    except Exception as e:
        error_msg = f"Error: {str(e)}"
        log(error_msg, style="bold red")
        log_tool_result(error_msg)
        return [error_msg]


@rnaseq_agent.tool
async def print_dependency_paths(ctx: RunContext[RNAseqData]) -> str:
    """Print out all paths and settings in the dependency object."""
    log_tool_header("print_dependency_paths")
    result = f"""
    Dependency paths:
    - fastq_dir: {ctx.deps.fastq_dir}
    - metadata_path: {ctx.deps.metadata_path}
    - kallisto_index_dir: {ctx.deps.kallisto_index_dir}
    - organism: {ctx.deps.organism}
    - output_dir: {ctx.deps.output_dir}
    - tx2gene_path: {ctx.deps.tx2gene_path}
    """
    log_tool_result(result)
    return result

#############################
# Section: Kallisto Quantification Tools
#############################


@rnaseq_agent.tool
async def find_kallisto_index(ctx: RunContext[RNAseqData]) -> str:
    """
    Search for and return the file path of an appropriate Kallisto index based on the organism specified in RNAseqData.

    Inputs:
      - ctx: RunContext[RNAseqData]
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
        console.log(
            f"[bold blue]Finding Kallisto index for organism:[/] {ctx.deps.organism}")
        organism = ctx.deps.organism.lower()
        index_dir = ctx.deps.kallisto_index_dir

        # Look for index files in the specified directory
        index_files = await find_files(ctx, ctx.deps.kallisto_index_dir, '.idx')

        if not index_files:
            return f"Error: No Kallisto index files found in {index_dir}"

        # Try to find an index matching the organism
        matching_indices = [idx for idx in index_files if organism in os.path.basename(
            os.path.dirname(idx)).lower()]

        if matching_indices:
            index_path = matching_indices[0]
            # Also find the transcript-to-gene mapping file if available
            tx2gene_files = await find_files(ctx, os.path.dirname(index_path), '.txt')
            if tx2gene_files:
                t2g_files = [f for f in tx2gene_files if any(
                    x in os.path.basename(f).lower() for x in ['t2g', 'tx2gene'])]
                if t2g_files:
                    ctx.deps.tx2gene_path = t2g_files[0]

            return f"Found Kallisto index for {organism}: {index_path}"
        else:
            # If no organism-specific index found, return the first one
            return f"No index specific to {organism} found. Using the first available index: {index_files[0]}"

    except Exception as e:
        return f"Error finding Kallisto index: {str(e)}"


@rnaseq_agent.tool
async def run_kallisto_quantification(ctx: RunContext[RNAseqData], parallel: bool = True, max_workers: Optional[int] = None, threads_per_job: Optional[int] = None) -> str:
    """
    Run Kallisto quantification on paired-end FASTQ files with adaptive parallelization and record the resulting abundance files.

    Inputs:
      - ctx: RunContext[RNAseqData]
          Must include:
             • fastq_dir: The directory where FASTQ files are stored.
             • kallisto_index_dir: The directory containing the Kallisto index files.
             • output_dir: The directory to store Kallisto outputs.
             • organism: Organism information to select the correct index.
      - parallel: Whether to enable parallel processing (default: True)
      - max_workers: Maximum number of concurrent Kallisto processes (default: auto-determined)
      - threads_per_job: Number of threads per Kallisto instance (default: auto-determined)

    Process:
      1. Logs current FASTQ directory and full context details.
      2. Uses find_files to search for FASTQ files with the suffix "fastq.gz" in the designated fastq_dir.
      3. If no FASTQ files are found, returns an error message.
      4. Calls find_kallisto_index to determine the correct index (and transcript-to-gene mapping if available).
      5. Extracts the index path from the returned result.
      6. Creates the output directory if it does not exist.
      7. Searches for paired FASTQ files using regular expression patterns for R1 and R2.
      8. Matches R1 with R2 files to form sample pairs.
      9. Determines optimal parallelization strategy based on available system resources.
      10. For each paired sample, builds the Kallisto command and executes it via subprocess (in parallel if enabled).
      11. Logs progress before and after running quantification for each sample.
      12. Collects the paths to generated abundance files (e.g., abundance.tsv) and stores them in ctx.deps.abundance_files.

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
        import concurrent.futures
        import multiprocessing

        console.log(
            f"[bold cyan]Fastq Directory from context:[/] {ctx.deps.fastq_dir}")
        console.log(
            f"[bold blue]Full context.deps details:[/]\n{vars(ctx.deps)}")

        # Find paired FASTQ files using the fastq_dir from the dependency
        fastq_files = await find_files(ctx, ctx.deps.fastq_dir, 'fastq.gz')
        if not fastq_files:
            return f"Error: No FASTQ files found in {ctx.deps.fastq_dir}"

        # Find the Kallisto index
        index_result = await find_kallisto_index(ctx)
        if "Error" in index_result:
            return index_result

        # Extract the index path from the result
        index_path = None
        for line in index_result.splitlines():
            if '.idx' in line:
                # Extract the path, which should be after the colon
                if ':' in line:
                    index_path = line.split(':', 1)[1].strip()
                else:
                    # If no colon, look for a path with .idx
                    words = line.split()
                    for word in words:
                        if '.idx' in word:
                            index_path = word
                            break

        if not index_path:
            return "Error: Could not determine Kallisto index path"

        # Create output directory
        output_dir = ctx.deps.output_dir
        os.makedirs(output_dir, exist_ok=True)

        # Organize FASTQ files into pairs based on naming conventions
        paired_files = {}

        # Identify pairs using common naming patterns
        r1_pattern = re.compile(r'.*_(R1|1)\.fastq\.gz$')
        r2_pattern = re.compile(r'.*_(R2|2)\.fastq\.gz$')

        r1_files = [f for f in fastq_files if r1_pattern.match(f)]
        r2_files = [f for f in fastq_files if r2_pattern.match(f)]

        # Match R1 with R2 files
        for r1_file in r1_files:
            # Convert R1 to R2 in the filename
            expected_r2 = r1_file.replace(
                '_R1', '_R2').replace('_1.fastq', '_2.fastq')
            if expected_r2 in r2_files:
                # Extract sample name from filename
                sample_name = os.path.basename(r1_file).split('_R1')[
                    0].split('_1.fastq')[0]
                paired_files[sample_name] = (r1_file, expected_r2)

        console.log(
            f"[bold yellow]Progress:[/] Identified {len(paired_files)} paired FASTQ file groups.")
        if not paired_files:
            # If no pairs found, check if files are single-end
            single_end = all(not r1_pattern.match(
                f) and not r2_pattern.match(f) for f in fastq_files)
            if single_end:
                return "Error: Single-end reads detected. Kallisto requires paired-end reads or additional parameters for single-end analysis."
            else:
                return "Error: Could not identify paired FASTQ files"

        # Determine parallelization strategy based on available resources
        total_cpus = multiprocessing.cpu_count()
        total_samples = len(paired_files)

        # Default parallelization strategy if not specified
        if max_workers is None:
            # Use at most half of CPUs for concurrent processes, but at least 1
            max_workers = min(total_samples, max(1, total_cpus // 2))

        if threads_per_job is None:
            # Distribute available CPUs among workers, but at least 1 thread per job
            threads_per_job = max(1, total_cpus // max_workers)

        console.log(f"[bold cyan]System resources:[/] {total_cpus} CPUs available")
        console.log(f"[bold cyan]Parallelization strategy:[/] Running up to {max_workers} samples simultaneously with {threads_per_job} threads each")

        # Function to process one sample
        def process_sample(sample_info):
            sample_name, (r1, r2) = sample_info
            sample_output_dir = os.path.join(output_dir, sample_name)
            os.makedirs(sample_output_dir, exist_ok=True)

            # Build Kallisto command
            cmd = [
                "kallisto", "quant",
                "--rf-stranded",  # Temporary hardcode until I implement determination of strandedness
                "-i", index_path,
                "-o", sample_output_dir,
                "-t", str(threads_per_job),  # Use calculated threads per sample
                "--plaintext",  # Output plaintext instead of HDF5
                r1, r2
            ]

            console.log(
                f"[bold yellow]Progress:[/] Running Kallisto quantification for sample: {sample_name}")
            process = subprocess.run(cmd, capture_output=True, text=True)

            console.log(
                f"[bold yellow]Progress:[/] Completed Kallisto run for sample: {sample_name} (return code: {process.returncode})")

            result = {
                "sample": sample_name,
                "success": process.returncode == 0,
                "message": f"Successfully processed {sample_name}" if process.returncode == 0 else f"Error processing {sample_name}: {process.stderr}",
                "abundance_file": os.path.join(sample_output_dir, "abundance.tsv") if process.returncode == 0 else None
            }

            return result

        results = []
        abundance_files = []

        # Use parallel processing if enabled and we have multiple samples
        if parallel and total_samples > 1:
            with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
                future_to_sample = {executor.submit(process_sample, item): item for item in paired_files.items()}
                for future in concurrent.futures.as_completed(future_to_sample):
                    try:
                        result = future.result()
                        results.append(result["message"])
                        if result["success"] and result["abundance_file"] and os.path.exists(result["abundance_file"]):
                            abundance_files.append(result["abundance_file"])
                    except Exception as exc:
                        sample_name = future_to_sample[future][0]
                        results.append(f"Error processing {sample_name}: {str(exc)}")
        else:
            # Sequential processing as fallback
            for sample_name, (r1, r2) in paired_files.items():
                result = process_sample((sample_name, (r1, r2)))
                results.append(result["message"])
                if result["success"] and result["abundance_file"] and os.path.exists(result["abundance_file"]):
                    abundance_files.append(result["abundance_file"])

        # Store the abundance file paths
        ctx.deps.abundance_files = abundance_files

        return f"""
Kallisto quantification completed for {len(results)} sample pairs.
Parallelization: {max_workers} concurrent samples with {threads_per_job} threads each.

Results:
{chr(10).join(results)}

Found {len(abundance_files)} abundance files for downstream analysis.
        """
    except Exception as e:
        return f"Error running Kallisto quantification: {str(e)}"

#############################
# Section: Differential Expression Analysis Tools
#############################


@rnaseq_agent.tool
async def prepare_edgeR_analysis(ctx: RunContext[RNAseqData]) -> str:
    """
    Prepare a sample mapping table for downstream edgeR differential expression analysis
    by matching Kallisto abundance files with sample metadata.

    Inputs:
      - ctx: RunContext[RNAseqData]
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
        log_tool_header("prepare_edgeR_analysis")

        # Check if we have abundance files
        if not ctx.deps.abundance_files:
            error_msg = "Error: No abundance files found. Please run Kallisto quantification first."
            log_tool_result(error_msg)
            return error_msg

        # Check if we have metadata
        if ctx.deps.metadata_df is None:
            error_msg = "Error: Metadata not loaded. Please run load_metadata first."
            log_tool_result(error_msg)
            return error_msg

        # Check if we have merged column
        if ctx.deps.merged_column is None:
            error_msg = "Error: Analysis column not identified. Please run identify_analysis_columns first."
            log_tool_result(error_msg)
            return error_msg

        # Create output directory if it doesn't exist
        os.makedirs(ctx.deps.output_dir, exist_ok=True)

        # Get sample names from abundance file paths
        sample_names = [os.path.basename(os.path.dirname(f))
                        for f in ctx.deps.abundance_files]
        log(f"Extracted {len(sample_names)} sample names from abundance files",
            level=LogLevel.VERBOSE)

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
            log_tool_result(warning_msg)
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
        log(f"Saved sample mapping to {analysis_df_path}",
            level=LogLevel.NORMAL)
        # Optionally, store in a registry for later reference:
        ctx.deps.file_registry = getattr(ctx.deps, 'file_registry', {})
        ctx.deps.file_registry['sample_mapping'] = analysis_df_path
        ctx.deps.sample_mapping = analysis_df_path
        log(f"Saved sample mapping to {analysis_df_path}",
            level=LogLevel.NORMAL)

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
        log_tool_result(result)
        return result
    except Exception as e:
        return f"Error running Kallisto quantification: {str(e)}"


@rnaseq_agent.tool
async def run_edger_limma_analysis(ctx: RunContext[RNAseqData]) -> str:
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
        log_tool_header("run_edger_limma_analysis")

        # Ensure the output directory exists
        os.makedirs(ctx.deps.output_dir, exist_ok=True)

        # Check if we have contrasts available
        if ctx.deps.contrast_matrix_df is None and not ctx.deps.contrasts:
            log("No contrasts defined. Please run process_metadata_with_agent first.",
                style="bold yellow")
            return "Error: No contrasts defined for differential expression analysis."

        # If we have contrasts but they're not in DataFrame format, convert them
        if ctx.deps.contrast_matrix_df is None and ctx.deps.contrasts:
            log("Converting contrasts from agent output to DataFrame format",
                level=LogLevel.VERBOSE)
            contrast_data = []
            for contrast in ctx.deps.contrasts.contrasts:
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
                    log(f"Saved contrasts to {contrast_path}", level=LogLevel.NORMAL)

        # Get the sample mapping path
        sample_mapping_path = ctx.deps.sample_mapping

        # Determine the tx2gene file argument; if not provided, pass "NA"
        tx2gene_arg = ctx.deps.tx2gene_path if ctx.deps.tx2gene_path and os.path.exists(
            ctx.deps.tx2gene_path) else "NA"

        # Use the main script
        main_r_script_path = "../script_development/data_analysis_agent/RNAseq.R"

        # Check if the R script exists
        if not os.path.exists(main_r_script_path):
            # Try an alternative path if the first one doesn't exist
            main_r_script_path = "aim1/UORCA/script_development/experiments/sample_RNAseq.R"
            if not os.path.exists(main_r_script_path):
                return "Error: R script not found at expected paths. Please verify the R script location."

        # Execute the R script with the necessary arguments
        cmd = ['Rscript', main_r_script_path, sample_mapping_path,
               ctx.deps.merged_column, ctx.deps.output_dir, tx2gene_arg]

        # Add the contrasts CSV file as the fifth argument if it exists
        if ctx.deps.contrast_path and os.path.exists(ctx.deps.contrast_path):
            cmd.append(ctx.deps.contrast_path)

        log(f"Running R script: {' '.join(cmd)}", level=LogLevel.NORMAL)
        process = subprocess.run(cmd, capture_output=True, text=True)
        stdout = process.stdout
        stderr = process.stderr

        # Log the captured outputs for traceability
        log_tool_result(f"STDOUT:\n{stdout}")
        if stderr:
            log_tool_result(f"STDERR:\n{stderr}")

        # Check the return code
        if process.returncode != 0:
            return f"edgeR/limma analysis failed with return code: {process.returncode}. See error output above."

        # Find and record DEG result files
        deg_dir = os.path.join(ctx.deps.output_dir, "DEG")
        if os.path.exists(deg_dir):
            deg_files = [os.path.join(deg_dir, f) for f in os.listdir(deg_dir) if f.endswith('.csv')]
            if deg_files:
                ctx.deps.deg_results_path = deg_files[0] if len(deg_files) == 1 else deg_files
                log(f"Found {len(deg_files)} differential expression result files", level=LogLevel.NORMAL)

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
        log(error_msg, style="bold red")
        return error_msg


#############################
# SECTION: Metadata Parsing Agent tools
#############################

def clean_string(ctx: RunContext[RNAseqData], s: str) -> str:
    """
    (Copied from MetadataAgent)
    Normalize and clean an input string by removing non-ASCII characters, extra whitespace, and unwanted symbols.
    """
    if pd.isna(s):
        return "NA"
    s = str(s).strip()
    s = unidecode(s)
    s = s.replace(" ", "_")
    s = re.sub(r'[^\w]', '', s)
    return s


# ----------------------------
# Imports
# ----------------------------

# ----------------------------
# Using logging functions already defined above
# ----------------------------
# Note: LogLevel, log, log_tool_header, and log_tool_result are already defined earlier in the file

# ----------------------------
# Dependency Class
# ----------------------------


@dataclass
class MetadataContext:
    """Container for metadata analysis data with enhanced context tracking."""
    metadata_path: str
    metadata_df: Optional[pd.DataFrame] = None
    # This will hold the final analysis column name (or merged version)
    merged_column: Optional[str] = None
    # Store the unique groups found in the analysis column
    unique_groups: Optional[List[str]] = None
    # To store designed contrasts if needed
    contrast_details: Optional[Dict[str, Any]] = None
    # Track all iterations of contrasts for reflection
    all_contrasts_iterations: List[Dict[str, Any]] = Field(default_factory=list)


# ----------------------------
# Load environment variables and initialize OpenAI client for metadata agent
# ----------------------------
# Note: We're using the same OpenAI client that was initialized earlier


# ----------------------------
# Create RNAseq metadata analysis agent
# ----------------------------
# Try reading your system prompt, otherwise use the fallback prompt
system_prompt_path = "../script_development/prompt_development/metadata_processing_systemprompt.txt"
try:
    with open(system_prompt_path, 'r') as f:
        system_prompt = f.read()

    # Print the first and last few lines for verification
    lines = system_prompt.split('\n')
    print(f"\n--- SYSTEM PROMPT FILE LOADED ---")
    print("First 3 lines:")
    for line in lines[:3]:
        print(f"> {line}")

    print("...")

    print("Last 3 lines:")
    for line in lines[-3:]:
        print(f"> {line}")
    print(f"--- END OF PREVIEW ({len(lines)} total lines) ---\n")
except Exception as e:
    system_prompt = """
    #### Integrated Prompt for Metadata Processing and Grouping Variable Selection

    You are provided with RNAseq metadata from different experiments. Your task is to identify the column(s) that contain biologically relevant information for differential expression analysis and to merge them into a single grouping variable if necessary. This grouping variable will be used in a single-factor analysis with edgeR/limma. In doing so, you must also evaluate each column to decide which ones provide informative biological variation and which ones should be excluded.

    General Guidelines:
    1. Focus on Biologically Relevant Information:
    • Include columns that capture sample-specific biological attributes such as tissue/disease type, genotype, or treatment conditions.
    • Exclude technical columns (e.g., sample IDs, run/experiment numbers) and those with no variation (all values identical) or with unique values that do not group samples.
    2. Merging Columns:
    • If more than one column is informative (e.g., one column for tissue type and one for treatment), merge these into a single composite grouping variable (for example, merged_analysis_group).
    • Ensure that the final grouping factor includes only information that is biologically significant for differential expression analysis.
    3. Output:
    • Return the name(s) of the final grouping column(s) and a brief explanation of your selection process and why the other columns were excluded.
    """
    print("Using fallback system prompt instead.")

metadata_agent = Agent(
    'openai:o4-mini',         # Use a powerful model
    deps_type=MetadataContext,
    system_prompt=system_prompt
)

# ----------------------------
# Utility function: Clean a string
# ----------------------------


@metadata_agent.tool
def clean_string(ctx: RunContext[MetadataContext], s: str) -> str:
    """
    Normalize and clean an input string by removing non-ASCII characters, extra whitespace, and unwanted symbols.
    """
    if pd.isna(s):
        return "NA"
    s = str(s).strip()
    s = unidecode(s)
    s = s.replace(" ", "_")
    s = re.sub(r'[^\w]', '', s)
    return s

# ----------------------------
# Tool 1: Process metadata
# ----------------------------


@metadata_agent.tool
async def process_metadata(ctx: RunContext[MetadataContext]) -> dict:
    """
    Load metadata from ctx.deps.metadata_path, clean all column names and
    cell values using clean_string, and store the cleaned DataFrame in the context.

    Basic preprocessing applied:
    - Removes Run/Experiment columns to avoid technical identifiers interfering with biological grouping
    - Removes columns with identical values (all values are the same)
    - Removes columns with all unique values (every row has a different value)
    - Cleans column names and values for consistency

    Returns information about the processed data.
    """
    try:
        log_tool_header("process_metadata", {"metadata_path": ctx.deps.metadata_path})

        # Load metadata based on file extension
        if ctx.deps.metadata_path.endswith('.csv'):
            df = pd.read_csv(ctx.deps.metadata_path)
        elif ctx.deps.metadata_path.endswith(('.tsv', '.txt')):
            df = pd.read_csv(ctx.deps.metadata_path, sep='\t')
        else:
            df = pd.read_csv(ctx.deps.metadata_path, sep=None, engine='python')

        # Remove Run and Experiment columns - these are technical columns that interfere with biological grouping
        df = df.loc[:, ~df.columns.str.contains('Run|Experiment', case=False)]

        # Remove duplicate rows
        df = df.drop_duplicates()

        # Remove columns where all values are identical
        df = df.loc[:, df.nunique() > 1]

        # Remove columns where all values are different (i.e. all values are unique)
        df = df.loc[:, df.nunique() < df.shape[0]]

        # Clean column names
        new_columns = {col: clean_string(ctx, col) for col in df.columns}
        df.rename(columns=new_columns, inplace=True)

        # Clean each cell value
        for col in df.columns:
            df[col] = df[col].apply(
                lambda x: clean_string(ctx, x) if pd.notna(x) else x)

        # Store cleaned metadata in context
        ctx.deps.metadata_df = df

        # Count unique values for each column
        column_stats = {}
        for col in df.columns:
            unique_vals = df[col].unique()
            column_stats[col] = {
                "unique_count": len(unique_vals),
                "values": list(unique_vals) if len(unique_vals) < 20 else list(unique_vals[:20])
            }

        nrows, ncols = df.shape
        summary = f"Metadata processed: {nrows} rows and {ncols} columns.\nColumns: {', '.join(df.columns)}"
        log_tool_result(summary)

        # Return information about the processed data
        return {
            "message": summary,
            "columns": list(df.columns),
            "column_stats": column_stats,
            "shape": (nrows, ncols),
            "error": False
        }
    except Exception as e:
        error_msg = f"Error processing metadata: {str(e)}"
        log(error_msg, style="bold red")
        log_tool_result(error_msg)
        return {"message": error_msg, "error": True, "columns": [], "column_stats": {}, "shape": (0, 0)}

# ----------------------------
# Tool 2: Merge Analysis Columns
# ----------------------------


@metadata_agent.tool
async def merge_analysis_columns(ctx: RunContext[MetadataContext], columns_input: Union[str, List[str], dict] = None) -> dict:
    """
    Merge specified columns into a single analysis column.

    This tool can accept input in different formats:
    1. JSON string with {"columns": [...]} structure
    2. List of column names directly
    3. A single column name as string

    If only one column is provided, it's set as the analysis column.
    If multiple columns are provided, they are merged by joining values with an underscore.

    The result is stored in ctx.deps.merged_column.
    """
    try:
        # Handle different input formats flexibly
        candidate_cols = []

        # Case 1: JSON string input
        if isinstance(columns_input, str):
            try:
                # Try to parse as JSON
                parsed = json.loads(columns_input)
                if isinstance(parsed, dict) and "columns" in parsed:
                    candidate_cols = parsed["columns"]
                elif isinstance(parsed, list):
                    candidate_cols = parsed
                else:
                    # Treat as a single column name
                    candidate_cols = [columns_input]
            except json.JSONDecodeError:
                # Not JSON, treat as a single column name
                candidate_cols = [columns_input]

        # Case 2: List input
        elif isinstance(columns_input, list):
            candidate_cols = columns_input

        # Case 3: Dict input with "columns" key
        elif isinstance(columns_input, dict) and "columns" in columns_input:
            candidate_cols = columns_input["columns"]

        # Validate we have something to work with
        if not candidate_cols:
            msg = "No columns provided for merging."
            log_tool_result(msg)
            return {"message": msg, "success": False, "merged_column": None}

        df = ctx.deps.metadata_df.copy()

        # Filter to ensure columns actually exist in the dataframe
        valid_cols = [col for col in candidate_cols if col in df.columns]
        if not valid_cols:
            msg = f"None of the specified columns {candidate_cols} exist in the metadata."
            log_tool_result(msg)
            return {"message": msg, "success": False, "merged_column": None}

        if len(valid_cols) == 1:
            # Single column; use it directly
            ctx.deps.merged_column = valid_cols[0]
            msg = f"Single column '{valid_cols[0]}' selected as the analysis column."
        else:
            # Multiple columns; merge them
            merged_col = "merged_analysis_group"
            df[merged_col] = df[valid_cols].astype(str).apply(
                lambda row: "_".join(row.values), axis=1)
            ctx.deps.metadata_df = df  # update the DataFrame
            ctx.deps.merged_column = merged_col
            msg = f"Multiple columns {', '.join(valid_cols)} merged into column '{merged_col}'."

            # Show a preview of the merged values
            unique_merged = df[merged_col].unique().tolist()
            preview = unique_merged[:5] if len(
                unique_merged) > 5 else unique_merged
            msg += f"\nMerged values (preview): {preview}"

        log_tool_result(msg)
        return {
            "message": msg,
            "success": True,
            "merged_column": ctx.deps.merged_column,
            "input_columns": valid_cols
        }
    except Exception as e:
        error_msg = f"Error in merge_analysis_columns: {str(e)}"
        log(error_msg, style="bold red")
        log_tool_result(error_msg)
        return {"message": error_msg, "success": False, "merged_column": None}

# ----------------------------
# Tool 3: Extract Unique Values
# ----------------------------


@metadata_agent.tool
async def extract_unique_values(ctx: RunContext[MetadataContext]) -> dict:
    """
    Simply extracts and returns the unique values from the selected analysis column.

    This tool only:
    1. Identifies unique values in the selected analysis column
    2. Stores these values in ctx.deps.unique_groups
    3. Returns a dictionary with the unique values

    No additional analysis or contrast generation is performed.
    """
    try:
        if ctx.deps.metadata_df is None:
            return {"success": False, "message": "Error: Metadata has not been processed.", "unique_values": []}

        if not ctx.deps.merged_column:
            return {"success": False, "message": "Error: Analysis column not defined. Please run merge_analysis_columns first.", "unique_values": []}

        df = ctx.deps.metadata_df
        analysis_col = ctx.deps.merged_column

        # Extract unique values from the analysis column
        unique_values = sorted(df[analysis_col].dropna().unique().tolist())
        ctx.deps.unique_groups = unique_values

        # Log the unique values found
        log_tool_result(
            f"Extracted {len(unique_values)} unique values from column '{analysis_col}'")

        # Return only the unique values with basic metadata
        return {
            "success": True,
            "column": analysis_col,
            "unique_values": unique_values,
            "count": len(unique_values)
        }

    except Exception as e:
        error_msg = f"Error extracting unique values: {str(e)}"
        log(error_msg, style="bold red")
        return {"success": False, "message": error_msg, "unique_values": []}


# ----------------------------
# Reflection Process Function
# ----------------------------
def run_reflection_process(metadata_deps: RNAseqData, initial_contrasts: Contrasts) -> Contrasts:
    log("Starting reflection process...", style="bold green")

    # Copy initial contrasts and track unique expressions
    all_contrasts = initial_contrasts.contrasts.copy()
    unique_exprs = {c.expression.strip().lower() for c in all_contrasts}

    metadata_deps.all_contrasts_iterations = []
    report_sections = []

    # Log the initial set
    metadata_deps.all_contrasts_iterations.append({
        "iteration": 0,
        "type": "initial",
        "contrasts": [{"name": c.name, "expression": c.expression} for c in all_contrasts]
    })
    report_sections.append("## INITIAL CONTRASTS ##\n" +
                           "\n".join(f"- {c.name}: {c.expression}" for c in all_contrasts))

    iteration = 1
    while True:
        reflection_prompt = f"""
I have the following {len(all_contrasts)} contrasts:

{chr(10).join(f"- {c.name}: {c.expression}" for c in all_contrasts)}

Unique groups: {metadata_deps.unique_groups}

Please suggest any additional biologically meaningful contrasts that might be missing.
If none are needed, explain why.
"""
        # Run reflection
        result: ReflectionResult = reflection_agent.run_sync(
            reflection_prompt,
            deps=MetadataContext,
            result_type=ReflectionResult
        )

        # Collect brand‑new contrasts
        new_list = []
        for c in result.contrasts:
            expr = c.expression.strip().lower()
            if expr not in unique_exprs:
                unique_exprs.add(expr)
                new_list.append(c)
                all_contrasts.append(c)

        # Record iteration
        metadata_deps.all_contrasts_iterations.append({
            "iteration": iteration,
            "type": "reflection",
            "new_contrasts": [{"name": c.name, "expression": c.expression} for c in new_list],
            "summary": result.summary or ""
        })

        # Update report
        if new_list:
            report_sections.append(
                f"## REFLECTION {iteration} ##\nAdded {len(new_list)}:\n" +
                "\n".join(f"- {c.name}: {c.expression}" for c in new_list)
            )
        else:
            report_sections.append(f"## REFLECTION {iteration} ##\nNo new contrasts added.")
            break

        iteration += 1
        if iteration > 3:
            report_sections.append("Reached max iterations.")
            break

    # Build final summary
    summary = "\n\n".join(report_sections)
    summary += (
        f"\n\nTotal initial: {len(initial_contrasts.contrasts)}\n"
        f"Total added: {len(all_contrasts) - len(initial_contrasts.contrasts)}\n"
        f"Final total: {len(all_contrasts)}\n"
        f"Iterations: {iteration - 1}"
    )

    return Contrasts(contrasts=all_contrasts, summary=summary)

# ----------------------------
# Updated metadata processing tool in Data Analysis Agent
# ----------------------------
@rnaseq_agent.tool
async def process_metadata_with_agent(ctx: RunContext[RNAseqData]) -> str:
    """
    1. Load full metadata (for sample mapping).
    2. Create a separate ContrastContext for running metadata_agent + reflection_agent.
    3. Generate initial contrasts, then reflect to augment them.
    4. Write out full set of contrasts to CSV and update ctx.deps.
    """
    log_tool_header("process_metadata_with_agent")

    # 1) Load full metadata into ctx.deps.metadata_df
    full = pd.read_csv(ctx.deps.metadata_path)
    # apply same cleaning everywhere
    for col in full.columns:
        full[col] = full[col].apply(lambda x: await clean_string(ctx, x) if pd.notna(x) else x)
    ctx.deps.metadata_df = full

    # 2) Prepare ContrastContext (drop Run/Experiment)
    df2 = full.loc[:, ~full.columns.str.contains('Run|Experiment', case=False)].drop_duplicates()
    contrast_ctx = ContrastContext(metadata_path=ctx.deps.metadata_path)
    contrast_ctx.metadata_df = df2

    # 3) Run metadata_agent to get initial contrasts
    prompt = "<your metadata prompt>"
    initial: Contrasts = metadata_agent.run_sync(prompt, deps=contrast_ctx, result_type=Contrasts)
    contrast_ctx.contrasts = initial.contrasts
    contrast_ctx.summary = initial.summary

    # 4) Extract unique groups
    _ = metadata_agent.run_sync("", deps=contrast_ctx, tool="extract_unique_values")

    # 5) Reflect to add new contrasts
    all_cons = initial.contrasts.copy()
    exprs = {c.expression.strip().lower() for c in all_cons}
    iteration = 1
    while True:
        # build reflection prompt
        refl_prompt = f"Reflect on contrasts: { [c.expression for c in all_cons] }"
        refl: ReflectionResult = reflection_agent.run_sync(
            refl_prompt, deps=contrast_ctx, result_type=ReflectionResult)
        new = [c for c in refl.contrasts if c.expression.strip().lower() not in exprs]
        if not new or iteration>3: break
        for c in new:
            exprs.add(c.expression.strip().lower()); all_cons.append(c)
        iteration += 1
    # Final Contrasts object
    enhanced = Contrasts(contrasts=all_cons, summary="<combined summary>")

    # 6) Save back to main context
    ctx.deps.contrasts = enhanced
    data = [c.dict() for c in all_cons]
    dfc = pd.DataFrame(data)
    os.makedirs(ctx.deps.output_dir, exist_ok=True)
    path = os.path.join(ctx.deps.output_dir, 'contrasts.csv')
    dfc.to_csv(path, index=False)
    ctx.deps.contrast_matrix_df = dfc
    ctx.deps.contrast_path = path

    return f"Metadata+contrast processing complete. {len(all_cons)} contrasts saved to {path}."


# ----------------------------
# Main Execution
# ----------------------------
if __name__ == "__main__":
    # Parse command-line arguments for the RNAseq analysis data
    parser = argparse.ArgumentParser(
        description="RNAseq analysis pipeline parameters for testing")
    parser.add_argument("--fastq_dir", type=str, default="./TestRNAseqData_SETBP1/GSE262710/fastq",
                        help="Directory where FASTQ files are stored")
    parser.add_argument("--metadata_path", type=str, default="./TestRNAseqData_SETBP1/GSE262710/GSE262710_metadata.csv",
                        help="Path to the metadata CSV file")
    parser.add_argument("--kallisto_index_dir", type=str, default="../../../data/kallisto_indices",
                        help="Directory where Kallisto index files are stored")
    parser.add_argument("--organism", type=str,
                        default="human", help="Organism name")
    parser.add_argument("--output_dir", type=str, default="./analysis_output/GSE262710",
                        help="Directory to save analysis output")
    args = parser.parse_args()

    # Ensure the output directory exists
    os.makedirs(args.output_dir, exist_ok=True)
    # Create a dedicated logs folder inside the output directory
    logs_dir = os.path.join(args.output_dir, "logs")
    os.makedirs(logs_dir, exist_ok=True)
    # Generate a timestamp string formatted as YYYYMMDD_HHMMSS
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    # Create the log file name with the timestamp embedded
    log_file_path = os.path.join(logs_dir, f"log_{timestamp}.txt")
    log_file = open(log_file_path, "a")

    # Define a Tee class to duplicate output to both STDOUT and the log file
    class Tee:
        def __init__(self, *streams):
            self.streams = streams

        def write(self, data):
            for s in self.streams:
                s.write(data)

        def flush(self):
            for s in self.streams:
                s.flush()

    import sys
    import datetime
    sys.stdout = Tee(sys.stdout, log_file)

    # Create an RNAseqData instance (renamed to analysis_data)
    analysis_data = RNAseqData(
        fastq_dir=args.fastq_dir,
        metadata_path=args.metadata_path,
        kallisto_index_dir=args.kallisto_index_dir,
        organism=args.organism,
        output_dir=args.output_dir

    )

    # Initialize conversation with analysis steps (using your testing prompt)
    initial_prompt = """
    Use the provided tools to perform an RNAseq analysis. This should encompass:
        1. Kallisto quantification, after identifying appropriate files and indices. Note that the index files, FASTQ files, and metadata are already provided, and you should not need to perform additional tasks to generate these - instead, locate them using the provided tools, using your judgement to determine if it is appropriate.
        2. Processing of the metadata - you should use the dedicated metadata agent to do this.
        3. Preparing the edgeR analysis
        4. Running edgeR analysis for differential expression, including contrasts and results.
        5. For the moment, do NOT perform a GSEA.
    """

    # Run the agent
    try:
        # Run the agent synchronously using the new dependency instance
        result = rnaseq_agent.run_sync(
            initial_prompt, deps=analysis_data)
        console.print(
            Panel("Analysis completed successfully!", style="bold green"))
        console.print("\n[bold yellow]Agent response:[/bold yellow]")
        console.print(result.data)
    except Exception as e:
        console.print(
            Panel(f"Error during analysis: {str(e)}", style="bold red"))

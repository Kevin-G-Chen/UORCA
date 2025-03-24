# Create a global console instance
# console = Console()
# Imports
# ----------------------------

import os
import glob
import re
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import date
from rich.console import Console
from rich.panel import Panel
from dataclasses import dataclass
from typing import List, Dict, Optional, Union, Tuple, Any, Literal
from unidecode import unidecode
from pydantic import BaseModel, Field
import matplotlib.pyplot as plt
import seaborn as sns
import gseapy as gp
import nest_asyncio
nest_asyncio.apply()
console = Console()
from pydantic_ai import Agent, RunContext
import logfire
logfire.configure(token=os.environ.get("LOGFIRE_KEY"))
logfire.instrument_openai()


# Configure logging levels
class LogLevel:
    MINIMAL = 0   # Only critical information
    NORMAL = 1    # Default level - tool calls, parameters, and results
    VERBOSE = 2   # Detailed information including context for each call
    DEBUG = 3     # Maximum information for debugging

# Set the current log level (change this to adjust verbosity)
CURRENT_LOG_LEVEL = LogLevel.NORMAL

def log(message, level=LogLevel.NORMAL, style=""):
    """Log a message if the current log level is equal to or greater than the specified level"""
    if CURRENT_LOG_LEVEL >= level:
        if style:
            console.print(message, style=style)
        else:
            console.print(message)

def log_tool_header(tool_name, params=None):
    """Print a clear header when a tool is called"""
    if CURRENT_LOG_LEVEL >= LogLevel.NORMAL:
        console.print("╭" + "─" * 100 + "╮")
        header = f"│ TOOL: {tool_name}"
        header += " " * (100 - len(header) - 1) + "│"
        console.print(header, style="bold blue")

        if params and CURRENT_LOG_LEVEL >= LogLevel.NORMAL:
            param_str = "│ Parameters: " + str(params)
            param_str += " " * (100 - len(param_str) - 1) + "│"
            console.print(param_str)

        console.print("├" + "─" * 100 + "┤")

def log_tool_result(result):
    """Print a clearly formatted tool result"""
    if CURRENT_LOG_LEVEL >= LogLevel.NORMAL:
        result_lines = str(result).strip().split('\n')
        for line in result_lines:
            line_str = "│ " + line
            line_str += " " * (100 - len(line_str) - 1) + "│"
            console.print(line_str)
        console.print("╰" + "─" * 100 + "╯")
        console.print("")  # Add extra line for separation

# ----------------------------
# Define the dependency type
# ----------------------------
@dataclass
class RNAseqData:
    """Container for RNAseq analysis data and paths."""
    # Input data
    fastq_dir: str
    metadata_path: str
    kallisto_index_dir: str
    organism: str = "human"  # default to human
    output_dir: str = "output"
    tx2gene_path: Optional[str] = None

    # Runtime data that gets populated during analysis
    metadata_df: Optional[pd.DataFrame] = None
    abundance_files: List[str] = None
    merged_column: Optional[str] = None
    contrast_groups: Dict[str, Dict[str, str]] = None
    sample_mapping: Optional[pd.DataFrame] = None

# ----------------------------
# Create an RNAseq analysis agent
# ----------------------------

# Adding temporary rough information for the agent until I develop a proper way of implementing documentation

rnaseq_agent = Agent(
    'openai:o3-mini', # Change to more powerful model
    deps_type=RNAseqData,
    system_prompt="""
    You are an expert RNAseq data analyst. Your task is to analyze RNAseq data using a series of bioinformatics tools.

    Follow these general principles throughout your analysis:
    1. Work systematically through the RNA-seq analysis workflow
    2. Validate inputs at each step
    3. Provide clear explanations of what's happening
    4. Handle errors gracefully
    5. Generate appropriate visualizations when needed
    6. Be comprehensive, both in the analysis steps but also more routine steps. For example, if you cannot find a file, ensure you check other common file extensions.
    7. After completing each step, take careful note of any output files. Specifically, make note of the location and names of saved files, and ensure these are added to context.
    """
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
      2. Uses os.walk to recursively traverse the directory structure and check every file's name.
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
            log(f"Initial Context.deps details:\n{vars(ctx.deps)}", level=LogLevel.VERBOSE, style="bold blue")
            setattr(ctx.deps, '_logged_context', True)

        # Log tool call with parameters
        log_tool_header("find_files", {"directory": directory, "suffix": suffix})

        # Log context details only in verbose mode
        log(f"Context.deps details:\n{vars(ctx.deps)}", level=LogLevel.VERBOSE, style="bold blue")
        if hasattr(ctx, "message_history"):
            log(f"Message History: {ctx.message_history}", level=LogLevel.DEBUG, style="bold magenta")

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
async def run_kallisto_quantification(ctx: RunContext[RNAseqData]) -> str:
    """
    Run Kallisto quantification on paired-end FASTQ files and record the resulting abundance files.

    Inputs:
      - ctx: RunContext[RNAseqData]
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
        console.log(f"[bold cyan]Fastq Directory from context:[/] {ctx.deps.fastq_dir}")
        console.log(f"[bold blue]Full context.deps details:[/]\n{vars(ctx.deps)}")

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
            expected_r2 = r1_file.replace('_R1', '_R2').replace('_1.fastq', '_2.fastq')
            if expected_r2 in r2_files:
                # Extract sample name from filename
                sample_name = os.path.basename(r1_file).split('_R1')[0].split('_1.fastq')[0]
                paired_files[sample_name] = (r1_file, expected_r2)

        console.log(f"[bold yellow]Progress:[/] Identified {len(paired_files)} paired FASTQ file groups.")
        if not paired_files:
            # If no pairs found, check if files are single-end
            single_end = all(not r1_pattern.match(f) and not r2_pattern.match(f) for f in fastq_files)
            if single_end:
                return "Error: Single-end reads detected. Kallisto requires paired-end reads or additional parameters for single-end analysis."
            else:
                return "Error: Could not identify paired FASTQ files"

        # Run Kallisto for each pair
        results = []
        for sample_name, (r1, r2) in paired_files.items():
            sample_output_dir = os.path.join(output_dir, sample_name)
            os.makedirs(sample_output_dir, exist_ok=True)

            # Build Kallisto command
            cmd = [
                "kallisto", "quant",
                "-i", index_path,
                "-o", sample_output_dir,
                "-t", "4",  # Use 4 threads
                "--plaintext",  # Output plaintext instead of HDF5
                "--bootstrap-samples=10",  # Number of bootstrap samples, low at the moment to speed up testing
                r1, r2
            ]

            console.log(f"[bold yellow]Progress:[/] Running Kallisto quantification for sample: {sample_name}")
            process = subprocess.run(cmd, capture_output=True, text=True)

            console.log(f"[bold yellow]Progress:[/] Completed Kallisto run for sample: {sample_name} (return code: {process.returncode})")
            if process.returncode == 0:
                results.append(f"Successfully processed {sample_name}")
            else:
                results.append(f"Error processing {sample_name}: {process.stderr}")

        # Collect paths to abundance files
        abundance_files = []
        for sample_name in paired_files.keys():
            abundance_file = os.path.join(output_dir, sample_name, "abundance.tsv")
            if os.path.exists(abundance_file):
                abundance_files.append(abundance_file)

        # Store the abundance file paths
        ctx.deps.abundance_files = abundance_files

        return f"""
Kallisto quantification completed for {len(results)} sample pairs.

Results:
{chr(10).join(results)}

Found {len(abundance_files)} abundance files for downstream analysis.
        """
    except Exception as e:
        return f"Error running Kallisto quantification: {str(e)}"


# ----------------------------
# Main Execution
# ----------------------------
if __name__ == "__main__":
    # Create data instance for GSE262710
    test_data = RNAseqData(
        fastq_dir="./notebooks/Testing/PydanticAI/TestRNAseqData_SETBP1/GSE262710/fastq",
        metadata_path="./notebooks/Testing/PydanticAI/TestRNAseqData_SETBP1/GSE262710/GSE262710_metadata.csv",
        kallisto_index_dir="./data/kallisto_indices",
        organism="human",
        output_dir="./notebooks/Testing/PydanticAI/analysis_output/GSE262710"
    )
    test_data2 = RNAseqData( # Temporary data instance for testing
        fastq_dir="./TestRNAseqData_SETBP1/GSE262710/fastq",
        metadata_path="./TestRNAseqData_SETBP1/GSE262710/GSE262710_metadata.csv",
        kallisto_index_dir="../../../data/kallisto_indices",
        organism="human",
        output_dir="./analysis_output/GSE262710"
    )

    # Initialize conversation with analysis steps
    initial_prompt = """
    Please perform a Kallisto quantification with the information provided. Note that this is an initial test case, and so you will not be expected to perform the remaining analysis steps.
    """

    # Run the agent
    try:
        # Run the agent synchronously
        result = rnaseq_agent.run_sync(initial_prompt, deps=test_data2)
        console.print(Panel("Analysis completed successfully!", style="bold green"))
        console.print("\n[bold yellow]Agent response:[/bold yellow]")
        console.print(result.data)
    except Exception as e:
        console.print(Panel(f"Error during analysis: {str(e)}", style="bold red"))

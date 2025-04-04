
# Create a global console instance
# console = Console()
# Imports
# ----------------------------

import logfire
from pydantic_ai import Agent, RunContext
import argparse
import os
import glob
import re
import subprocess
import pandas as pd
import numpy as np
import json
from dotenv import load_dotenv
import matplotlib.pyplot as plt
from datetime import date
from rich.console import Console
from rich.panel import Panel
from dataclasses import dataclass
from typing import List, Dict, Optional, Union, Tuple, Any, Literal
from unidecode import unidecode
from pydantic import BaseModel, Field, ConfigDict
import matplotlib.pyplot as plt
import nest_asyncio
from openai import OpenAI
nest_asyncio.apply()
console = Console()
# logfire.configure(token=os.environ.get("LOGFIRE_KEY"))
# logfire.instrument_openai()


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
    if CURRENT_LOG_LEVEL >= LogLevel.NORMAL:
        console.print(f"TOOL: {tool_name}", style="bold blue")
        if params:
            console.print(f"Parameters: {params}")


def log_tool_result(result):
    if CURRENT_LOG_LEVEL >= LogLevel.NORMAL:
        console.print(result)

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
    kallisto_index_path: Optional[str] = None
    tx2gene_path: Optional[str] = None

    # Runtime data that gets populated during analysis
    metadata_df: Optional[pd.DataFrame] = None
    abundance_files: List[str] = None
    merged_column: Optional[str] = None
    contrast_groups: Dict[str, Dict[str, str]] = None
    sample_mapping: Optional[pd.DataFrame] = None


# Load environment variables from .env file
load_dotenv()

# Get the API key
openai_api_key = os.getenv("OPENAI_API_KEY")
# Define Pydantic class for the structured output
client = OpenAI()

# Define Pydantic class for the structured output


class RelevantColumns(BaseModel):
    model_config = ConfigDict(extra="forbid")
    columns: List[str] = Field(
        description="Column, or columns, that are best suited for grouping variables in the metadata table.")


schema = RelevantColumns.model_json_schema()
# ----------------------------
# Create an RNAseq analysis agent
# ----------------------------

# Adding temporary rough information for the agent until I develop a proper way of implementing documentation

# Read system prompt from file
system_prompt_path = "../script_development/prompt_development/data_analysis_agent_systemprompt.txt"
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
    deps_type=RNAseqData,
    system_prompt=system_prompt
)


@rnaseq_agent.tool
async def run_gsea_analysis(ctx: RunContext[RNAseqData], contrast_name: str) -> str:
    """
    Run Gene Set Enrichment Analysis (GSEA) on the normalized expression data for a specified contrast.

    Inputs:
      - ctx: RunContext[RNAseqData]
          The execution context containing an RNAseqData dependency instance. This object carries paths and runtime
          attributes such as:
            • fastq_dir: Directory for FASTQ files.
            • metadata_path: Path to the metadata CSV.
            • kallisto_index_dir: Directory where Kallisto index files are stored.
            • output_dir: Directory for saving output files.
            • Other optional or runtime-populated attributes (e.g. metadata_df, abundance_files, merged_column, contrast_groups).
      - contrast_name (str):
          A string identifier for the contrast, which is used to look up the specific experimental conditions
          (numerator and denominator) within the contrast_groups attribute.

    Process:
      1. Logs key information from the dependency context.
      2. Checks for the existence of a normalized counts CSV file (generated from edgeR analysis).
      3. Loads expression data from the file.
      4. Retrieves the contrast details (numerator and denominator groups) from ctx.deps.contrast_groups.
      5. Constructs a class vector for samples by comparing sample names in the expression data with entries in the
         metadata (using the merged_column).
      6. Calls the gseapy.gsea function with the expression data, contrast, and additional parameters.
      7. Generates plots (e.g. a top‐5 pathways plot) and saves them to the designated output directory.
      8. Constructs and returns a detailed summary string describing the performed analysis, including the number of
         pathways analyzed and the significant pathways identified.
    Output:
      Returns a multiline string summary that includes:
         • The contrast name.
         • Total number of pathways analyzed.
         • Number of significantly enriched pathways (e.g. FDR q-val < 0.25).
         • A brief display (via head) of the top enriched pathways.
         • A list of output directories and file names generated by the analysis.

    Purpose in pipeline:
      This tool is used in the later stages of the RNAseq analysis workflow to perform pathway enrichment analysis
      after differential expression has been quantified.
    """
    try:
        console.log(
            f"[bold blue]Tool Called:[/] run_gsea_analysis with contrast_name: {contrast_name}")
        if hasattr(ctx, "message_history"):
            console.log(
                f"[bold magenta]Message History:[/] {ctx.message_history}")
        console.log(f"[bold blue]Context.deps details:[/]\n{vars(ctx.deps)}")
        console.log(f"[bold cyan]Fastq Directory:[/] {ctx.deps.fastq_dir}")
        if hasattr(ctx, "message_history"):
            console.log(
                f"[bold magenta]Message History:[/] {ctx.message_history}")

        # Check if we have DEG results
        deg_file = os.path.join(ctx.deps.output_dir, f"deg_{contrast_name}.csv")
        if not os.path.exists(deg_file):
            msg = f"Error: DEG results file '{deg_file}' not found. Please run differential expression analysis first."
            console.log(f"[bold red]Tool Error:[/] {msg}")
            return msg

        deg_df = pd.read_csv(deg_file)
        console.log(f"[bold yellow]Progress:[/] Loaded DEG data with shape: {deg_df.shape}")

        # Build the rank list from the DEG CSV. The DEG file is assumed to have columns "Gene" and "logFC".
        rnk = deg_df[['Gene', 'logFC']].copy()
        rnk['Gene'] = rnk['Gene'].str.upper()  # ensure gene symbols are uppercase
        rnk = rnk.dropna().sort_values("logFC", ascending=False)
        console.log(f"[bold yellow]Progress:[/] Constructed rank list with {rnk.shape[0]} genes")

        # Define an output directory that clearly indicates the contrast and that this is a preranked GSEA run.
        this_gsea_out_dir = os.path.join(ctx.deps.output_dir, f"GSEA_{contrast_name}_prerank")
        os.makedirs(this_gsea_out_dir, exist_ok=True)
        console.log(f"[bold yellow]Progress:[/] Created GSEA output directory: {this_gsea_out_dir}")

        # Run preranked GSEA using the fixed GMT file
        gmt_path = "/data/tki_agpdev/kevin/phd/aim1/UORCA/scratch/msigdb/c2.all.v2024.1.Hs.symbols.gmt"
        pre_res = gp.prerank(
            rnk=rnk,
            gene_sets=gmt_path,
            permutation_num=1000,  # adjust if needed
            outdir=this_gsea_out_dir,
            seed=42,
            verbose=True
        )

        # Save the complete GSEA results to CSV.
        all_out = os.path.join(this_gsea_out_dir, f"{contrast_name}_gsea_results_all.csv")
        pre_res.res2d.to_csv(all_out)
        console.log(f"[bold green]Saved complete GSEA results to:{all_out}")

        # If the expected column is present, filter for significance.
        if "FDR q-val" in pre_res.res2d.columns:
            sig = pre_res.res2d[pre_res.res2d["FDR q-val"] < 0.05]
            sig_out = os.path.join(this_gsea_out_dir, f"{contrast_name}_gsea_results_sig.csv")
            sig.to_csv(sig_out)
            console.log(f"[bold green]Saved significant GSEA results to:{sig_out}")
            sig_msg = f"{sig.shape[0]} significant gene sets found"
        else:
            sig_msg = "No FDR q-val column found; significant results not extracted"

        msg = f"""GSEA preranked analysis completed for contrast: {contrast_name}
Total gene sets tested: {pre_res.res2d.shape[0]}
{sig_msg}
Complete results saved to: {all_out}
"""
        console.log(f"[bold green]Tool Completed:[/] run_gsea_analysis for contrast: {contrast_name}")
        return msg

        # Generate plots
        terms = gs_res.res2d.Term
        gs_res.plot(
            terms=terms[:5],
            show_ranking=True,
            ofname=os.path.join(ctx.deps.output_dir,
                                f"gsea_{contrast_name}_top5.png")
        )

        console.log(
            f"[bold yellow]Progress:[/] Generated GSEA plots in directory: {os.path.join(ctx.deps.output_dir, f'gsea_{contrast_name}')}")
        # Summarize results
        sig_pathways = gs_res.res2d[gs_res.res2d['FDR q-val'] < 0.25]
        msg = f"""
 GSEA Analysis completed for contrast: {contrast_name}

 Summary:
 - Total pathways analyzed: {len(gs_res.res2d)}
 - Significant pathways (FDR < 0.25): {len(sig_pathways)}
 - Top enriched pathways:
 {sig_pathways[['Term', 'NES', 'FDR q-val']].head().to_string()}

 Generated files:
 - GSEA results: gsea_{contrast_name}/
 - Top pathways plot: gsea_{contrast_name}_top5.png
 """
        console.log(
            f"[bold green]Tool Completed:[/] run_gsea_analysis for contrast: {contrast_name}")
        return msg

    except Exception as e:
        error_msg = f"Error running GSEA analysis: {str(e)}"
        console.log(f"[bold red]Tool Exception:[/] {error_msg}")
        return error_msg

    finally:
        # This finally block will always run, ensuring that any open figures are closed.
        plt.close('all')

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
def clean_string(ctx: RunContext[RNAseqData], s: str) -> str:
    """
    Normalize and clean an input string by removing non-ASCII characters, redundant white space, and unwanted symbols.

    Inputs:
      - ctx: RunContext[RNAseqData]
          Provides the current context (although not used for processing in this tool, it conforms to tool signature).
      - s (str):
          The string that is to be cleaned. If s is missing or contains special characters, these are normalized.

    Process:
      1. Checks if the input is missing (NaN) and returns "NA" if so.
      2. Converts the input to a string and trims leading/trailing whitespace.
      3. Uses unidecode to normalize special characters to standard ASCII.
      4. Replaces spaces with underscores.
      5. Uses a regular expression to remove any remaining non-word characters.

    Output:
      Returns the cleaned and normalized string.

    Purpose in pipeline:
      This tool is used to ensure consistency in sample naming across metadata and quantification outputs,
      which is vital for matching and merging data from different sources.
    """
    if pd.isna(s):
        return "NA"  # Handle missing values
    # Convert to string and remove leading/trailing whitespaces
    s = str(s).strip()
    s = unidecode(s)  # Normalize special characters to ASCII
    s = s.replace(" ", "_")  # Replace spaces with underscores
    s = re.sub(r'[^\w]', '', s)  # Remove non-word characters
    return s


@rnaseq_agent.tool
async def process_metadata(ctx: RunContext[RNAseqData]) -> str:
    """
    Load metadata from the file specified in ctx.deps.metadata_path, remove columns where all values are identical,
    clean column names using clean_string, and identify candidate grouping columns for differential expression analysis.

    If multiple candidate columns are found (based on common biological keywords), they are merged into a new column
    named "merged_analysis_group". If only one candidate is identified, it is used directly.

    The processed metadata and the selected grouping column are stored in ctx.deps.metadata_df and ctx.deps.merged_column,
    respectively.
    """
    try:
        log_tool_header("process_metadata", {
                        "metadata_path": ctx.deps.metadata_path})

        # Load metadata based on file extension
        if ctx.deps.metadata_path.endswith('.csv'):
            df = pd.read_csv(ctx.deps.metadata_path)
        elif ctx.deps.metadata_path.endswith('.tsv') or ctx.deps.metadata_path.endswith('.txt'):
            df = pd.read_csv(ctx.deps.metadata_path, sep='\t')
        else:
            df = pd.read_csv(ctx.deps.metadata_path, sep=None, engine='python')

        # Remove columns where all values are the same
        df = df.loc[:, df.nunique() > 1]

        # Clean column names using the clean_string tool function
        new_columns = {col: clean_string(ctx, col) for col in df.columns}
        df.rename(columns=new_columns, inplace=True)

        # Clean all cell values in the dataframe using clean_string
        for col in df.columns:
            df[col] = df[col].apply(
                lambda x: clean_string(ctx, x) if pd.notna(x) else x)

        # Store cleaned metadata in the context
        ctx.deps.metadata_df = df

        # Identify candidate grouping columns passing metadata via LLM

        # Prepare inputs for LLM prompt

        preview = ctx.deps.metadata_df.to_csv(index=False)
        columns = list(ctx.deps.metadata_df.columns)

        # Prepare LLM prompt

        prompt = f"""
        You are an expert in RNAseq metadata analysis. Given the following metadata table, determine the column, or columns, that are most suitable for grouping variables in the context of differential expression analysis.

        Please note that you should focus only on groups that contain biologically interested information. For example, columns like 'sample_id' or 'replicate' are not typically used as grouping variables. Furthermore, "differentiation_experiment," while good for quality control, does not contribute any biological findings, and should not be used as a grouping variable. In contrast, columns such as "genotype" "celltype" and the like are often used as grouping variables.

        Do not select redundant columns. For example, if one column contains values such as "WT" and another contains "Wild Type," you should only select one of them, as these columns are functionally equivalent.

        Metadata columns: {columns}
        Metadata preview: {preview}
        """

        response = client.responses.create(
            model="gpt-4o",
            input=prompt,
            text={
                "format": {
                    "type": "json_schema",
                    "name": "analysis_columns",
                    "schema": schema,
                    "strict": True
                }
            }
        )

        candidate_cols = json.loads(response.output_text)['columns']

        # Merge candidates if more than one is found; otherwise, use the single candidate
        if len(candidate_cols) > 1:
            merged_col = "merged_analysis_group"
            df[merged_col] = df[candidate_cols].apply(
                lambda row: '_'.join([str(val) for val in row.values]), axis=1)
            ctx.deps.merged_column = merged_col
            result_msg = f"Merged candidate columns ({', '.join(candidate_cols)}) into '{merged_col}'."
        else:
            ctx.deps.merged_column = candidate_cols[0]
            result_msg = f"Selected candidate column '{candidate_cols[0]}' as the grouping column."

        # Update the metadata in the context
        ctx.deps.metadata_df = df

        summary = f"Metadata processed: {df.shape[0]} samples, {df.shape[1]} columns.\n{result_msg}"
        log_tool_result(summary)
        return summary
    except Exception as e:
        error_msg = f"Error processing metadata: {str(e)}"
        log(error_msg, style="bold red")
        log_tool_result(error_msg)
        return error_msg


@rnaseq_agent.tool
async def design_contrasts(ctx: RunContext[RNAseqData]) -> str:
    """
    Design experimental contrasts for differential expression analysis based on the merged group information from metadata.

    Inputs:
      - ctx: RunContext[RNAseqData]
          Must include:
            • A metadata DataFrame (ctx.deps.metadata_df) with merged grouping specified in ctx.deps.merged_column.

    Process:
      1. Validates that the metadata and merged analysis column are present.
      2. Extracts the unique group identifiers from the merged column.
      3. Based on the number of groups, creates either a simple contrast (for two groups) or, for multiple groups, identifies
         a reference (control) group and generates pairwise or reference-based contrasts.
      4. Constructs a dictionary (contrast_details) that maps each contrast name to its details (numerator, denominator, expression).
      5. Stores the contrast details in ctx.deps.contrast_groups.
      6. Returns a detailed summary of the designed contrasts, including the groups used and a table-like representation of
         the contrast details.

    Output:
      A multiline string report detailing:
         • The unique groups identified in the metadata.
         • The contrast names with corresponding details.
         • Recommendations for using these contrasts in downstream differential expression analysis.

    Purpose in pipeline:
      Designing contrasts is a critical step that defines comparisons between experimental conditions for subsequent
      differential expression analysis (e.g., by edgeR).
    """
    try:
        if ctx.deps.metadata_df is None:
            return "Error: Metadata not loaded. Please run load_metadata first."

        if ctx.deps.merged_column is None:
            return "Error: Analysis column not identified. Please run identify_analysis_columns first."

        metadata_df = ctx.deps.metadata_df
        group_col = ctx.deps.merged_column

        # Get unique groups
        groups = metadata_df[group_col].unique().tolist()

        # Design contrasts based on the groups
        contrasts = []
        contrast_details = {}

        if len(groups) == 2:
            # Simple case with only two groups - one contrast
            contrasts.append(f"{groups[1]}_vs_{groups[0]}")
            contrast_details[f"{groups[1]}_vs_{groups[0]}"] = {
                'name': f"{groups[1]} vs {groups[0]}",
                'numerator': groups[1],
                'denominator': groups[0],
                'expression': f"{groups[1]} - {groups[0]}"
            }
        else:
            # For multiple groups, compare each to a reference
            # Try to identify a control/reference group
            potential_controls = [g for g in groups if any(kw in g.lower() for kw in
                                                           ['control', 'ctrl', 'reference', 'ref', 'normal', 'wt', 'wild', 'mock', 'dmso', 'pbs', 'untreated'])]

            if potential_controls:
                # Use the identified control
                reference = potential_controls[0]
                other_groups = [g for g in groups if g != reference]

                for group in other_groups:
                    contrast_name = f"{group}_vs_{reference}"
                    contrasts.append(contrast_name)
                    contrast_details[contrast_name] = {
                        'name': f"{group} vs {reference}",
                        'numerator': group,
                        'denominator': reference,
                        'expression': f"{group} - {reference}"
                    }
            else:
                # If no obvious control, pick the first group as reference
                reference = groups[0]
                other_groups = groups[1:]

                for group in other_groups:
                    contrast_name = f"{group}_vs_{reference}"
                    contrasts.append(contrast_name)
                    contrast_details[contrast_name] = {
                        'name': f"{group} vs {reference}",
                        'numerator': group,
                        'denominator': reference,
                        'expression': f"{group} - {reference}"
                    }

                # Also create pairwise comparisons between non-reference groups if there aren't too many
                if len(other_groups) <= 4:  # Limit to avoid too many comparisons
                    for i, group1 in enumerate(other_groups):
                        for group2 in other_groups[i+1:]:
                            contrast_name = f"{group1}_vs_{group2}"
                            contrasts.append(contrast_name)
                            contrast_details[contrast_name] = {
                                'name': f"{group1} vs {group2}",
                                'numerator': group1,
                                'denominator': group2,
                                'expression': f"{group1} - {group2}"
                            }

        # Store the contrasts
        ctx.deps.contrast_groups = contrast_details

        return f"""
Designed contrasts based on {group_col} column with values: {', '.join(groups)}

Contrasts:
{pd.DataFrame([contrast_details[c] for c in contrasts]).to_string()}
        """
    except Exception as e:
        return f"Error designing contrasts: {str(e)}"


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

# ----------------------------
# Kallisto Quantification Tools
# ----------------------------


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

        # Run Kallisto for each pair
        results = []
        for sample_name, (r1, r2) in paired_files.items():
            sample_output_dir = os.path.join(output_dir, sample_name)
            os.makedirs(sample_output_dir, exist_ok=True)

            # Build Kallisto command
            cmd = [
                "kallisto", "quant",
                "--rf-stranded",  # Temporary hardcode until I implement determination of strandedness
                "-i", index_path,
                "-o", sample_output_dir,
                "-t", "24",  # Use 4 threads
                "--plaintext",  # Output plaintext instead of HDF5
                r1, r2
            ]

            console.log(
                f"[bold yellow]Progress:[/] Running Kallisto quantification for sample: {sample_name} via command: {cmd}")
            process = subprocess.run(cmd, capture_output=True, text=True)

            console.log(
                f"[bold yellow]Progress:[/] Completed Kallisto run for sample: {sample_name} (return code: {process.returncode})")
            if process.returncode == 0:
                results.append(f"Successfully processed {sample_name}")
            else:
                results.append(
                    f"Error processing {sample_name}: {process.stderr}")

        # Collect paths to abundance files
        abundance_files = []
        for sample_name in paired_files.keys():
            abundance_file = os.path.join(
                output_dir, sample_name, "abundance.tsv")
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
# Differential Expression Analysis Tools
# ----------------------------


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
    Minimal tool to run an edgeR/limma analysis using a single R script.

    This tool performs the following steps:
      1. Loads required libraries (using pacman) and reads the metadata CSV from ctx.deps.metadata_path.
      2. Expects the metadata CSV to include an 'abundance_file' column with Kallisto output paths.
      3. Imports quantification data via tximport. If a tx2gene mapping file is provided in ctx.deps.tx2gene_path,
         it is used; otherwise, the script runs without it.
      4. Creates a DGEList and attaches the metadata.
      5. Normalizes the DGEList using calcNormFactors.
      6. Constructs a design matrix using the grouping column specified in ctx.deps.merged_column.
      7. Performs a voom transformation and fits a linear model with limma.
      8. If exactly two groups are present, a contrast is computed (group2 - group1); if more than two groups,
         top tables for each coefficient are generated.
      9. Saves differential expression results (CSV files) and the normalized DGEList (RDS file).

    Logging messages are printed to STDOUT/STDERR during the execution of the R script.

    Note: This is a baseline analysis template and does not include additional filtering or plotting steps.
    """
    try:
        log_tool_header("run_edger_limma_analysis")

        # Ensure the output directory exists
        os.makedirs(ctx.deps.output_dir, exist_ok=True)
        r_script_path = os.path.join(
            ctx.deps.output_dir, "edger_limma_analysis.R")

        # Determine the tx2gene file argument; if not provided, pass "NA"
        tx2gene_arg = ctx.deps.tx2gene_path if ctx.deps.tx2gene_path and os.path.exists(
            ctx.deps.tx2gene_path) else "NA"

        # Write the R script
        with open(r_script_path, "w") as f:
            f.write('''
user_lib <- Sys.getenv("R_LIBS_USER", unset="~/R/library")
if (!dir.exists(user_lib)) {
    dir.create(user_lib, recursive = TRUE)
}
.libPaths(c(user_lib, .libPaths()))

library(pacman)
# Load only the required packages (removing tidyverse)
p_load(edgeR, limma, tximport)

cat("=== R Script: edgeR/limma Analysis Start ===\n")

args <- commandArgs(trailingOnly = TRUE)
metadata_file <- args[1]
merged_group <- args[2]
output_dir <- args[3]
tx2gene_file <- args[4]

cat("Loading metadata from:", metadata_file, "\n")
metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)
str(metadata)
if(nrow(metadata) == 0) {
  stop("Metadata is empty!")
}

if(!"abundance_file" %in% colnames(metadata)) {
  stop("Metadata must contain an 'abundance_file' column!")
}

# Load tx2gene mapping if provided
if(tx2gene_file != "NA" && file.exists(tx2gene_file)){
  cat("Loading tx2gene mapping from:", tx2gene_file, "\n")
  tx2gene <- read.delim(tx2gene_file, header = FALSE, stringsAsFactors = FALSE)
  # Select columns 1 and 3 (like dplyr::select(1,3))
  tx2gene <- tx2gene[, c(1, 3)]
  # Rename columns (like setNames)
  names(tx2gene) <- c("TXNAME", "GENEID")
  # Filter rows where GENEID is not NA and not empty (like dplyr::filter)
  tx2gene <- tx2gene[!is.na(tx2gene$GENEID) & tx2gene$GENEID != "", ]
  str(tx2gene)
  use_tx2gene <- TRUE
} else {
  cat("No valid tx2gene file provided. Proceeding without tx2gene mapping.\n")
  use_tx2gene <- FALSE
}

cat("Importing Kallisto quantification data...\n")
if(use_tx2gene){
  kallisto <- tximport(metadata$abundance_file, type = "kallisto",
                       tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM",
                       ignoreAfterBar = TRUE)
} else {
  kallisto <- tximport(metadata$abundance_file, type = "kallisto",
                       txOut = TRUE, ignoreAfterBar = TRUE)
}

cat("Creating DGEList...\n")
DGE <- DGEList(counts = kallisto$counts)
# Combine DGE$samples with metadata (bind_cols replaced with cbind)
DGE$samples <- cbind(DGE$samples, metadata)

cat("Normalizing DGEList...\n")
DGE.norm <- calcNormFactors(DGE)

cat("Creating design matrix using grouping column:", merged_group, "\n")
design <- model.matrix(as.formula(paste("~0 +", merged_group)), data = DGE.norm$samples)
colnames(design) <- sub(merged_group, "", colnames(design))
cat("Design matrix:\n")
print(design)

cat("Performing voom transformation...\n")
v <- voom(DGE.norm, design, plot = FALSE)

cat("Fitting linear model...\n")
fit <- lmFit(v, design)
fit <- eBayes(fit)


if(ncol(design) == 2){
  cat("Exactly two groups detected. Calculating contrast (group2 - group1)...\n")
  contrast_name <- paste(colnames(design)[2], "-", colnames(design)[1])
  contrast <- makeContrasts(diff = contrast_name, levels = design)
  cat("Contrast matrix:")
  print(contrast)
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2)
  cat("Top differential expression results for contrast:\n")
  top_results <- topTable(fit2, number = Inf)
  print(head(top_results))
  top_results$Gene <- rownames(top_results)
  top_results <- top_results[, c("Gene", setdiff(names(top_results), "Gene"))]
  print(head(top_results))
  write.csv(top_results, file = file.path(output_dir, "DEG_results.csv"), row.names = FALSE)

} else {
  cat("Multiple groups detected. Generating top results for each coefficient...\n")
  for(i in 1:ncol(design)){
    coef_name <- colnames(design)[i]
    cat("Top results for", coef_name, ":\n")
    top_results <- topTable(fit, coef = i, number = Inf)

    top_results$Gene <- rownames(top_results)
    top_results <- top_results[, c("Gene", setdiff(names(top_results), "Gene"))]
    print(head(top_results))
    write.csv(top_results, file = file.path(output_dir, paste0("DEG_results_", coef_name, ".csv")), row.names = FALSE)
  }
}

cat("Saving normalized DGEList object...\n")
saveRDS(DGE.norm, file = file.path(output_dir, "DGE_norm.RDS"))

cat("=== R Script: edgeR/limma Analysis Completed ===\n")
            ''')
        os.chmod(r_script_path, 0o755)
        log(f"Created R script at {r_script_path}", level=LogLevel.NORMAL)

        # Execute the R script. Pass the metadata path, merged column, output directory, and tx2gene argument.
        cmd = ['Rscript', r_script_path, ctx.deps.sample_mapping,
               ctx.deps.merged_column, ctx.deps.output_dir, tx2gene_arg]
        print("Running R script for edgeR/limma analysis:", ' '.join(cmd))
        process = subprocess.run(cmd, capture_output=True, text=True)
        stdout = process.stdout
        stderr = process.stderr

        # Log the captured outputs for traceability
        log_tool_result(f"STDOUT:\n{stdout}")
        log_tool_result(f"STDERR:\n{stderr}")

        return f"edgeR/limma analysis completed with return code: {process.returncode}"

    except Exception as e:
        error_msg = f"Error in run_edger_limma_analysis: {str(e)}"
        log(error_msg, style="bold red")
        return error_msg

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
    log_file_path = os.path.join(args.output_dir, "log.txt")
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
        2. Preparation for the edgeR differential expression analysis, including sample mapping and metadata analysis
        3. Running edgeR analysis for differential expression, including contrasts and results.
        4. After performing the edgeR analysis, use these results to perform a GSEA.
    """

    # Run the agent
    try:
        # Run the agent synchronously using the new dependency instance
        result = rnaseq_agent.run_sync(initial_prompt, deps=analysis_data)
        console.print(
            Panel("Analysis completed successfully!", style="bold green"))
        console.print("\n[bold yellow]Agent response:[/bold yellow]")
        console.print(result.data)
    except Exception as e:
        console.print(
            Panel(f"Error during analysis: {str(e)}", style="bold red"))

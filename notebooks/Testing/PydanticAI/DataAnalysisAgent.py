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

rnaseq_agent = Agent(
    'openai:gpt-4o-mini',
    deps_type=RNAseqData,
    system_prompt="""
    You are an expert RNAseq data analyst. Your task is to analyze RNAseq data using a series of bioinformatics tools.

    Follow these general principles:
    1. Work systematically through the RNA-seq analysis workflow
    2. Validate inputs at each step
    3. Provide clear explanations of what's happening
    4. Handle errors gracefully
    5. Generate appropriate visualizations when needed
    6. Be comprehensive, both in the analysis steps but also more routine steps. For example, if you ncannot find a file, ensure you check other common file extensions.
    """
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
       2. Checks for the existence of a normalized counts CSV file (generated from DESeq2 analysis).
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
         console.log(f"[bold blue]Tool Called:[/] run_gsea_analysis with contrast_name: {contrast_name}")
         if hasattr(ctx, "message_history"):
             console.log(f"[bold magenta]Message History:[/] {ctx.message_history}")
         console.log(f"[bold blue]Context.deps details:[/]\n{vars(ctx.deps)}")
         console.log(f"[bold cyan]Fastq Directory:[/] {ctx.deps.fastq_dir}")
         if hasattr(ctx, "message_history"):
             console.log(f"[bold magenta]Message History:[/] {ctx.message_history}")

         # Check if we have normalized expression data
         norm_counts_file = os.path.join(ctx.deps.output_dir, "DESeq2_normalized_counts.csv")
         if not os.path.exists(norm_counts_file):
             msg = "Error: Normalized counts file not found. Please run DESeq2 analysis first."
             console.log(f"[bold red]Tool Error:[/] {msg}")
             return msg

         console.log(f"[bold yellow]Progress:[/] Loaded expression data with shape: {expr_data.shape}")
         # Load normalized counts
         expr_data = pd.read_csv(norm_counts_file, index_col=0)

         # Get contrast details
         contrast = ctx.deps.contrast_groups[contrast_name]
         group_a = contrast['numerator']
         group_b = contrast['denominator']

         # Create class vector
         class_vector = []
         for sample in expr_data.columns:
             if sample in ctx.deps.metadata_df[ctx.deps.metadata_df[ctx.deps.merged_column] == group_a].index:
                 class_vector.append(group_a)
             else:
                 class_vector.append(group_b)

         # Run GSEA
         gs_res = gp.gsea(
             data=expr_data,
             gene_sets='MSigDB_Hallmark_2020',
             cls=class_vector,
             permutation_type='phenotype',
             permutation_num=1000,
             outdir=os.path.join(ctx.deps.output_dir, f"gsea_{contrast_name}"),
             method='signal_to_noise',
             threads=4
         )

         # Generate plots
         terms = gs_res.res2d.Term
         gs_res.plot(
             terms=terms[:5],
             show_ranking=True,
             ofname=os.path.join(ctx.deps.output_dir, f"gsea_{contrast_name}_top5.png")
         )

         console.log(f"[bold yellow]Progress:[/] Generated GSEA plots in directory: {os.path.join(ctx.deps.output_dir, f'gsea_{contrast_name}')}")
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
         console.log(f"[bold green]Tool Completed:[/] run_gsea_analysis for contrast: {contrast_name}")
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
async def load_metadata(ctx: RunContext[RNAseqData]) -> str:
    """
    Load, validate, and store metadata from the file specified in ctx.deps.metadata_path into the RNAseqData context for later use in the analysis.

    This tool automatically uses the metadata_path from the dependency context - you don't need to provide a path parameter.

    Inputs:
      - ctx: RunContext[RNAseqData]
          Contains the RNAseqData dependency where the loaded metadata DataFrame will be stored (in metadata_df).

    Process:
      1. Determines file format by examining the file extension (.csv, .tsv, .txt) and reads the file accordingly.
      2. Loads the metadata into a pandas DataFrame.
      3. Stores the DataFrame in ctx.deps.metadata_df for downstream analyses.
      4. Computes a list of "useful columns" (columns with more than one unique value) to help filter out trivial data.
      5. Returns a detailed summary string with the number of samples and columns, and lists the useful columns and a sample
         of the dataset.

    Output:
      A string message detailing:
         • Successful metadata loading.
         • DataFrame dimensions (number of samples and columns).
         • A list of columns with variability.
         • A preview (first few rows) of the metadata.

    Purpose in pipeline:
      Loading metadata is critical for linking sample IDs to experimental groups and merging with quantification results for
      downstream differential expression and pathway analyses.
    """
    try:
        log_tool_header("load_metadata", {"metadata_path": ctx.deps.metadata_path})

        # Load metadata based on file extension
        if ctx.deps.metadata_path.endswith('.csv'):
            metadata_df = pd.read_csv(ctx.deps.metadata_path)
        elif ctx.deps.metadata_path.endswith('.tsv') or ctx.deps.metadata_path.endswith('.txt'):
            metadata_df = pd.read_csv(ctx.deps.metadata_path, sep='\t')
        else:
            metadata_df = pd.read_csv(ctx.deps.metadata_path, sep=None, engine='python')

        # Store metadata and compute useful columns
        ctx.deps.metadata_df = metadata_df
        useful_cols = metadata_df.loc[:, metadata_df.nunique() > 1].columns.tolist()

        result = f"Metadata loaded: {metadata_df.shape[0]} samples, {metadata_df.shape[1]} columns. Useful columns: {', '.join(useful_cols)}."

        # Add a preview in verbose mode
        if CURRENT_LOG_LEVEL >= LogLevel.VERBOSE:
            result += f"\n\nPreview of metadata:\n{metadata_df.head().to_string()}"

        log_tool_result(result)
        return result
    except Exception as e:
        error_msg = f"Error loading metadata: {str(e)}"
        log(error_msg, style="bold red")
        log_tool_result(error_msg)
        return error_msg

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
    s = str(s).strip()  # Convert to string and remove leading/trailing whitespaces
    s = unidecode(s)  # Normalize special characters to ASCII
    s = s.replace(" ", "_")  # Replace spaces with underscores
    s = re.sub(r'[^\w]', '', s)  # Remove non-word characters
    return s

# ----------------------------
# Metadata Analysis Tools
# ----------------------------
@rnaseq_agent.tool
async def identify_analysis_columns(ctx: RunContext[RNAseqData]) -> str:
    """
    Analyze the metadata DataFrame to determine which columns should be used for grouping in the differential expression analysis.

    Inputs:
      - ctx: RunContext[RNAseqData]
          Must contain a loaded metadata DataFrame (ctx.deps.metadata_df) which provides experiment details.

    Process:
      1. Checks that metadata is loaded; if not, issues an error.
      2. Identifies variable columns (i.e. columns with more than one unique value).
      3. Evaluates each variable column to classify it as a potential biological factor or a technical factor
         using predefined keyword lists.
      4. Determines whether multiple biological factors are present, suggesting a merge of columns, or if one clear factor exists.
      5. Sets the ctx.deps.merged_column variable based on either merging the columns or selecting the most appropriate single column.
      6. Returns a detailed report that includes detected biological and technical factors, the variables considered,
         the recommendation, and whether merging is needed.

    Output:
      A detailed multiline string that outlines:
         • The list of variable columns.
         • The identified biological and technical factors.
         • The recommendation on merging and the final column decision.
         • An explanation of the grouping strategy.

    Purpose in pipeline:
      Correctly identifying and consolidating experimental group labels from metadata is essential before designing contrasts
      for downstream differential expression analysis.
    """
    try:
        console.log(f"[bold blue]Context.deps:[/] {ctx.deps}")
        if hasattr(ctx, "message_history"):
            console.log(f"[bold magenta]Message History:[/] {ctx.message_history}")
        if ctx.deps.metadata_df is None:
            return "Error: Metadata not loaded. Please run load_metadata first."

        # Get columns with variability
        metadata_df = ctx.deps.metadata_df
        variable_cols = metadata_df.loc[:, metadata_df.nunique() > 1].columns.tolist()

        # Analyze potential biological factors
        biological_factors = []
        technical_factors = []

        # Common keywords for biological factors
        bio_keywords = ['treatment', 'condition', 'genotype', 'disease', 'cell', 'tissue',
                         'time', 'dose', 'age', 'sex', 'gender', 'strain', 'group']

        # Common keywords for technical factors
        tech_keywords = ['batch', 'run', 'lane', 'library', 'seq', 'date', 'id', 'rep']

        for col in variable_cols:
            col_lower = col.lower()
            is_bio = any(keyword in col_lower for keyword in bio_keywords)
            is_tech = any(keyword in col_lower for keyword in tech_keywords)

            if is_bio:
                biological_factors.append(col)
            elif is_tech:
                technical_factors.append(col)

        # Determine if columns should be merged
        merge_needed = len(biological_factors) > 1
        if merge_needed:
            recommendation = f"Multiple biological factors detected: {', '.join(biological_factors)}. Consider merging these columns for analysis."
            cols_to_merge = biological_factors
        else:
            if len(biological_factors) == 1:
                recommendation = f"One clear biological factor detected: {biological_factors[0]}. This can be used directly."
                cols_to_merge = [biological_factors[0]]
            else:
                # If no obvious biological factors, suggest the columns with the fewest unique values
                n_unique = {col: metadata_df[col].nunique() for col in variable_cols}
                sorted_cols = sorted(n_unique.items(), key=lambda x: x[1])
                cols_to_merge = [sorted_cols[0][0]]
                if len(sorted_cols) > 1 and sorted_cols[1][1] < 10:  # Only suggest merging if second column has few unique values
                    cols_to_merge.append(sorted_cols[1][0])
                    merge_needed = True
                    recommendation = f"No clear biological factors detected. Suggesting to merge columns with fewest unique values: {', '.join(cols_to_merge)}."
                else:
                    recommendation = f"No clear biological factors detected. Suggesting to use the column with fewest unique values: {cols_to_merge[0]}."

        # Store the recommendation in the context
        if merge_needed:
            ctx.deps.merged_column = "merged_analysis_group"
        else:
            ctx.deps.merged_column = cols_to_merge[0]

        return f"""
Analysis of metadata columns:
- Biological factors: {', '.join(biological_factors) if biological_factors else 'None clearly identified'}
- Technical factors: {', '.join(technical_factors) if technical_factors else 'None clearly identified'}
- Variable columns: {', '.join(variable_cols)}

Recommendation: {recommendation}
Columns to use: {', '.join(cols_to_merge)}
Merge needed: {merge_needed}
        """
    except Exception as e:
        return f"Error analyzing metadata columns: {str(e)}"

@rnaseq_agent.tool
async def merge_metadata_columns(ctx: RunContext[RNAseqData], columns: List[str], new_column_name: str = "merged_analysis_group") -> str:
    """
    Merge one or more metadata columns into a single new column for downstream analysis.

    Inputs:
      - ctx: RunContext[RNAseqData]
          Provides access to the metadata DataFrame (ctx.deps.metadata_df) that will be modified.
      - columns (List[str]):
          A list of metadata column names to be merged. The function checks that each specified column exists.
      - new_column_name (str, default "merged_analysis_group"):
          The name for the new merged column that will be created in the metadata DataFrame.

    Process:
      1. Validates that metadata has been loaded and that specified columns exist.
      2. If only one column is provided, simply renames and cleans it; if multiple columns are provided, concatenates their
         cleaned (using clean_string) values with underscores as separators.
      3. Updates the metadata DataFrame and sets the ctx.deps.merged_column attribute.
      4. Computes a list of unique values in the newly created column.
      5. Returns a detailed message describing:
             • What columns were merged (or renamed) and to what new column.
             • The number of unique groups found.
             • A summary of the counts per group.

    Output:
      A multiline string message detailing the results of the column merge, including the unique group values and their counts.

    Purpose in pipeline:
      Merging columns is often necessary when multiple metadata fields contribute to the definition of experimental groups;
      this unified grouping is then used for designing contrasts for differential expression analysis.
    """
    try:
        if ctx.deps.metadata_df is None:
            return "Error: Metadata not loaded. Please run load_metadata first."

        console.log(f"[bold blue]Context.deps:[/] {ctx.deps}")
        if hasattr(ctx, "message_history"):
            console.log(f"[bold magenta]Message History:[/] {ctx.message_history}")
        metadata_df = ctx.deps.metadata_df

        # Check that all columns exist
        missing_cols = [col for col in columns if col not in metadata_df.columns]
        if missing_cols:
            return f"Error: Columns not found in metadata: {', '.join(missing_cols)}"

        if len(columns) == 1:
            # Just rename the column if only one provided
            metadata_df[new_column_name] = metadata_df[columns[0]].apply(lambda x: clean_string(ctx, x))
            result_message = f"Renamed and cleaned column {columns[0]} to {new_column_name}."
        else:
            # Create the merged column by concatenating values with underscores
            metadata_df[new_column_name] = metadata_df[columns].apply(
                lambda row: '_'.join([clean_string(ctx, val) for val in row.values.astype(str)]),
                axis=1
            )
            result_message = f"Merged columns {', '.join(columns)} into new column {new_column_name}."

        # Update the metadata
        ctx.deps.metadata_df = metadata_df
        ctx.deps.merged_column = new_column_name

        # Get unique values in the merged column
        unique_values = metadata_df[new_column_name].unique().tolist()

        return f"""
{result_message}

The merged column '{new_column_name}' contains {len(unique_values)} unique values:
{', '.join(unique_values)}

Sample counts per group:
{metadata_df[new_column_name].value_counts().to_string()}
        """
    except Exception as e:
        return f"Error merging metadata columns: {str(e)}"

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
      differential expression analysis (e.g., by DESeq2).
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
        console.log(f"[bold blue]Finding Kallisto index for organism:[/] {ctx.deps.organism}")
        organism = ctx.deps.organism.lower()
        index_dir = ctx.deps.kallisto_index_dir

        # Look for index files in the specified directory
        index_files = await find_files(ctx, ctx.deps.kallisto_index_dir, '.idx')

        if not index_files:
            return f"Error: No Kallisto index files found in {index_dir}"

        # Try to find an index matching the organism
        matching_indices = [idx for idx in index_files if organism in os.path.basename(os.path.dirname(idx)).lower()]

        if matching_indices:
            index_path = matching_indices[0]
            # Also find the transcript-to-gene mapping file if available
            tx2gene_files = await find_files(ctx, os.path.dirname(index_path), '.txt')
            if tx2gene_files:
                t2g_files = [f for f in tx2gene_files if any(x in os.path.basename(f).lower() for x in ['t2g', 'tx2gene'])]
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
                "--bootstrap-samples=10",  # Number of bootstrap samples
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
      6. Saves this mapping DataFrame to a CSV file (named "deseq2_analysis_samples.csv") in the output directory.

    Output:
      Returns a detailed multiline string that includes:
         • The number of samples successfully mapped.
         • A preview of the sample-to-metadata mapping.
         • A summary of the group counts based on the merged analysis column.

    Purpose in pipeline:
      This tool bridges the quantification and differential expression steps by ensuring that each
      sample's expression data is accurately linked to its experimental metadata, a prerequisite for DESeq2 analysis.
      It assumes Kallisto quantification has already been completed in a previous step.
    """
    try:
        log_tool_header("prepare_deseq2_analysis")

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
        sample_names = [os.path.basename(os.path.dirname(f)) for f in ctx.deps.abundance_files]
        log(f"Extracted {len(sample_names)} sample names from abundance files", level=LogLevel.VERBOSE)

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
                matches = metadata_df[metadata_df[col].astype(str) == sample_name].index.tolist()
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
                    matches = metadata_df[metadata_df[col].astype(str).str.contains(sample_name, case=False, na=False)].index.tolist()
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

        # Create a DataFrame for DESeq2 analysis
        analysis_df = pd.DataFrame(index=list(matched_samples.keys()))

        # Add the file paths
        analysis_df['abundance_file'] = [matched_samples[s]['abundance_file'] for s in analysis_df.index]

        # Add the metadata
        for col in metadata_df.columns:
            analysis_df[col] = [metadata_df.loc[matched_samples[s]['metadata_row'], col] for s in analysis_df.index]

        # Save the analysis dataframe for later use
        analysis_df_path = os.path.join(ctx.deps.output_dir, "edger_analysis_samples.csv")
        analysis_df.to_csv(analysis_df_path)
        log(f"Saved sample mapping to {analysis_df_path}", level=LogLevel.NORMAL)
        # Optionally, store in a registry for later reference:
        ctx.deps.file_registry = getattr(ctx.deps, 'file_registry', {})
        ctx.deps.file_registry['sample_mapping'] = analysis_df_path
        ctx.deps.sample_mapping = analysis_df
        log(f"Saved sample mapping to {analysis_df_path}", level=LogLevel.NORMAL)
        # Store the sample mapping DataFrame in the runtime data
        ctx.deps.sample_mapping = analysis_df

        result = f"""
Successfully prepared data for DESeq2 analysis with {len(analysis_df)} samples.

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
async def run_edger_analysis(
    ctx: RunContext[RNAseqData],
    sample_mapping_file: str,
    contrast_names: Optional[List[str]] = None,
) -> str:
    """
    Run edgeR differential expression analysis for one or more contrasts using dynamic R code.

    Parameters:
      - sample_mapping_file (str):
          Path for the sample mapping CSV, typically "<output_dir>/edger_analysis_samples.csv".
      - contrast_names (Optional[List[str]]):
          A list of contrast names to be analyzed; if None, all defined contrasts are used.

    This tool creates a base R script that builds an edgeR DGEList from tximport data,
    applies dynamic filtering via filterByExpr (using the merged group column), and estimates dispersions.
    It then writes a separate R script for each contrast (using the provided numerator and denominator).
    """
    try:
        log_tool_header("run_deseq2_analysis", {"contrast_names": contrast_names, "sample_mapping_file": sample_mapping_file})

        # Validate that contrast groups exist
        if not ctx.deps.contrast_groups:
            error_msg = "Error: No contrast groups defined. Please run design_contrasts first."
            log_tool_result(error_msg)
            return error_msg

        # If no specific contrasts provided, use all defined contrasts
        if contrast_names is None:
            contrast_names = list(ctx.deps.contrast_groups.keys())
            log(f"No specific contrasts provided. Analyzing all {len(contrast_names)} contrasts: {contrast_names}", level=LogLevel.NORMAL)

        # Check that the sample mapping file exists
        if not os.path.exists(sample_mapping_file):
            error_msg = "Error: Sample mapping file not found. Please run prepare_deseq2_analysis first."
            log_tool_result(error_msg)
            return error_msg

        # ----------------------------
        # Create the base R script for edgeR analysis
        # ----------------------------
        base_r_script_path = os.path.join(ctx.deps.output_dir, "run_edger_base.R")
        with open(base_r_script_path, "w") as f:
            f.write(f'''
# Load required libraries
suppressMessages(library(tximport))
suppressMessages(library(edgeR))
suppressMessages(library(ggplot2))
suppressMessages(library(pheatmap))

# Set working directory explicitly
setwd("{os.path.abspath(ctx.deps.output_dir)}")
cat("DEBUG: Working directory set to:", getwd(), "\\n")

# Diagnostic: List contents of working directory
cat("DEBUG: Files in working directory:", paste(list.files(), collapse=", "), "\\n")

# Check that sample mapping file exists (using its basename)
sample_file <- "{os.path.basename(sample_mapping_file)}"
if (file.exists(sample_file)) {{
  cat("SUCCESS: Sample mapping file found at:", sample_file, "\\n")
  sample_info <- read.csv(sample_file, row.names=1)
}} else {{
  cat("ERROR: Sample mapping file not found at:", sample_file, "\\n")
  cat("Directory contents:", paste(list.files(), collapse=", "), "\\n")
  stop(paste("File not found:", sample_file))
}}

# Define the grouping column
group_col <- "{ctx.deps.merged_column}"

# Get abundance files from sample_info
files <- sample_info$abundance_file
names(files) <- rownames(sample_info)

# Import Kallisto data using tximport, with optional tx2gene mapping if available
{ "tx2gene <- read.csv(\"" + ctx.deps.tx2gene_path + "\", header=FALSE, sep=\"\\t\"); colnames(tx2gene) <- c(\"TXNAME\", \"GENEID\");" if ctx.deps.tx2gene_path else "# No tx2gene file specified" }
txi <- tximport(files, type="kallisto", { "tx2gene=tx2gene" if ctx.deps.tx2gene_path else "txOut=TRUE" })

# Create an edgeR DGEList using the imported counts
dge <- DGEList(counts=txi$counts)

# Assign group information from sample_info
dge$samples$group <- sample_info[[group_col]]

# Apply filtering using filterByExpr
keep <- filterByExpr(dge, group=dge$samples$group)
dge <- dge[keep, , keep.lib.sizes=FALSE]

# Calculate normalization factors
dge <- calcNormFactors(dge)

# Create design matrix and estimate dispersions
design <- model.matrix(~0 + group, data=dge$samples)
colnames(design) <- levels(as.factor(dge$samples$group))
dge <- estimateDisp(dge, design)

# Save the edgeR object for downstream analysis
saveRDS(dge, file="{os.path.join(ctx.deps.output_dir, 'dge.rds')}")
''')
        os.chmod(base_r_script_path, 0o755)
        log(f"Executing base R script for edgeR analysis: {base_r_script_path}", level=LogLevel.NORMAL)
        process = subprocess.run(['Rscript', base_r_script_path], capture_output=True, text=True)
        if process.returncode != 0:
            error_msg = f"Error running edgeR base analysis:\nSTDOUT:\n{process.stdout}\nSTDERR:\n{process.stderr}"
            log_tool_result(error_msg)
            return error_msg
        else:
            log_tool_result("Base edgeR analysis completed successfully.\nOutput:\n" + process.stdout)

        # ----------------------------
        # Process each contrast using dynamically generated R scripts
        # ----------------------------
        all_results = []
        for contrast_name in contrast_names:
            contrast = ctx.deps.contrast_groups[contrast_name]
            numerator = contrast['numerator']
            denominator = contrast['denominator']
            contrast_dir = os.path.join(ctx.deps.output_dir, f"edger_{contrast_name}")
            os.makedirs(contrast_dir, exist_ok=True)
            results_file = os.path.join(contrast_dir, f"{contrast_name}_results.csv")
            ma_plot_file = os.path.join(contrast_dir, f"{contrast_name}_MAplot.png")
            summary_stats_file = os.path.join(contrast_dir, "summary_stats.txt")
            contrast_r_script_path = os.path.join(ctx.deps.output_dir, f"run_edger_{contrast_name}.R")
            with open(contrast_r_script_path, "w") as f:
                f.write(f'''
suppressMessages(library(edgeR))
suppressMessages(library(ggplot2))

# Load the saved edgeR object
dge <- readRDS("{os.path.join(ctx.deps.output_dir, 'dge.rds')}")

# Re-create the design matrix
design <- model.matrix(~0 + group, data=dge$samples)
colnames(design) <- levels(as.factor(dge$samples$group))

# Fit the model using quasi-likelihood methods
fit <- glmQLFit(dge, design)

# Define the contrast: numerator vs denominator
contrast_vector <- makeContrasts(contrasts = "{numerator} - {denominator}", levels=design)

# Perform the quasi-likelihood F-test
qlf <- glmQLFTest(fit, contrast=contrast_vector)

# Get full results table
res <- topTags(qlf, n=Inf)$table

# Write the results to CSV
write.csv(as.data.frame(res), file="{results_file}")

# Generate an MA plot (using plotMD, which is analogous in edgeR)
png(filename="{ma_plot_file}")
plotMD(qlf, main="MA Plot: {contrast_name}")
dev.off()

# Save summary statistics
sink("{summary_stats_file}")
cat("Contrast: {contrast_name}\\n")
cat("Total genes tested: ", nrow(res), "\\n")
cat("Significant genes (FDR < 0.05): ", sum(res$FDR < 0.05, na.rm=TRUE), "\\n")
cat("Up-regulated genes (logFC > 0): ", sum(res$FDR < 0.05 & res$logFC > 0, na.rm=TRUE), "\\n")
cat("Down-regulated genes (logFC < 0): ", sum(res$FDR < 0.05 & res$logFC < 0, na.rm=TRUE), "\\n")
sink()

quit(save="no", status=0)
''')
            os.chmod(contrast_r_script_path, 0o755)
            log(f"Executing R script for edgeR analysis of {contrast_name}: {contrast_r_script_path}", level=LogLevel.NORMAL)
            process = subprocess.run(['Rscript', contrast_r_script_path], capture_output=True, text=True)
            if process.returncode != 0:
                error_msg = f"Error running edgeR analysis for {contrast_name}:\n{process.stderr}"
                log_tool_result(error_msg)
                all_results.append(error_msg)
            else:
                if os.path.exists(summary_stats_file):
                    with open(summary_stats_file, 'r') as f:
                        summary_stats = f.read()
                else:
                    summary_stats = "No summary statistics available."
                result = f"""
----- edgeR analysis for contrast: {contrast_name} ({numerator} vs {denominator}) -----

Results file: {results_file}
MA plot: {ma_plot_file}

{summary_stats}
"""
                all_results.append(result)

        combined_results = "\n\n" + "="*80 + "\n\n".join(all_results) + "\n\n" + "="*80
        overall_summary = f"""
edgeR analysis completed for {len(contrast_names)} contrasts:
{', '.join(contrast_names)}

Base edgeR object saved as: {os.path.join(ctx.deps.output_dir, 'dge.rds')}
"""
        final_result = overall_summary + combined_results
        log_tool_result(final_result)
        return final_result

    except Exception as e:
        error_msg = f"Error running edgeR analysis: {str(e)}"
        log(error_msg, style="bold red")
        log_tool_result(error_msg)
        return error_msg
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
    Please analyze the RNA-seq data by calling the following tools in sequence:

    1. First, use the print_dependency_paths tool to see all the paths available in the context.

    2. Next, use the find_files tool to locate FASTQ files. DO NOT pass "ctx.deps.fastq_dir" as a string - instead, use the actual directory value . The suffix parameter should be 'fastq.gz'.

    3. Use the find_kallisto_index tool to locate the Kallisto index for human.

    4. Use the run_kallisto_quantification tool to perform quantification on the identified FASTQ files.

    5. Use the load_metadata tool to load the metadata file.

    6. Use the identify_analysis_columns tool to analyze the metadata columns.

    7. If needed, use the merge_metadata_columns tool with the columns identified in the previous step.

    8. Use the design_contrasts tool to set up contrasts for differential expression analysis.

    9. Use the prepare_edgeR_analysis tool to prepare for edgeR analysis.

    10. Use the run_edger_analysis tool to run edgeR analysis for the contrasts.

    After each step, provide a brief summary of what was done and what will be done next. Make sure to call each tool explicitly with the correct parameter values, not variable references.
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
@rnaseq_agent.tool
async def check_file_existence(ctx: RunContext[RNAseqData], filepath: str) -> str:
    """
    Check if a file exists at the given filepath. If not, report the contents of the directory.
    
    Inputs:
      - ctx: RunContext[RNAseqData] (provides context and directory info)
      - filepath (str): Absolute file path to check.
    
    Output:
      A string message indicating success or an error with diagnostic directory contents.
    """
    log_tool_header("check_file_existence", {"filepath": filepath})
    if os.path.exists(filepath):
        result = f"SUCCESS: File exists at {filepath}"
    else:
        dir_path = os.path.dirname(filepath)
        if os.path.exists(dir_path):
            contents = os.listdir(dir_path)
            result = f"ERROR: File {os.path.basename(filepath)} not found in {dir_path}.\nDirectory contents: {contents}"
        else:
            result = f"ERROR: Directory {dir_path} does not exist."
    log_tool_result(result)
    return result

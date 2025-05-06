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
from shared import AnalysisContext, AnalysisContext
from shared.workflow_logging import log_tool
import gseapy
from unidecode import unidecode
from openai import OpenAI
import matplotlib.pyplot as plt
import nest_asyncio
from concurrent.futures import ThreadPoolExecutor, as_completed

@dataclass
class MetadataContext:
    metadata_path: str
    metadata_df: Optional[pd.DataFrame] = None
    merged_column: Optional[str] = None
    unique_groups: Optional[List[str]] = None
    contrast_details: Optional[Dict[str, Any]] = None


#############################
# SECTION: Metadata Parsing Agent tools
#############################
# Try reading your system prompt, otherwise use the fallback prompt
system_prompt_path = "./main_workflow/prompts/analysis.txt"
try:
    with open(system_prompt_path, 'r') as f:
        system_prompt = f.read()

    # Print the first and last few lines for verification
    lines = system_prompt.split('\n')
    print(f"\n--- METADATA AGENT SYSTEM PROMPT FILE LOADED ---")
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
    - Removes columns with identical values
    - Removes columns with all unique values
    - Cleans column names and values for consistency

    Returns information about the processed data.
    """
    try:
        log_tool_header("process_metadata", {
                        "metadata_path": ctx.deps.metadata_path})

        # Load metadata based on file extension
        if ctx.deps.metadata_path.endswith('.csv'):
            df = pd.read_csv(ctx.deps.metadata_path)
        elif ctx.deps.metadata_path.endswith(('.tsv', '.txt')):
            df = pd.read_csv(ctx.deps.metadata_path, sep='\t')
        else:
            df = pd.read_csv(ctx.deps.metadata_path, sep=None, engine='python')

        # Print message of metadata dimensions
        log(f"Loaded metadata with {df.shape[0]} rows and {df.shape[1]} columns")
        # Remove Run and Experiment columns - these are added in a previous metadata processing step, so these are removed to more accurately emulate the original metadata
        # df = df.loc[:, ~df.columns.str.contains('Run|Experiment', case=False)]
        # Remove duplicate rows
        df = df.drop_duplicates()
        # Remove columns where all values are identical
        df = df.loc[:, df.nunique() > 1]
        # Remove columns where all values are different (i.e. all values are unique)
        # df = df.loc[:, df.nunique() < df.shape[0]]

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

        summary = f"""Metadata processed: {df.shape[0]} rows and {df.shape[1]} columns.
Columns: {', '.join(df.columns)}
"""
        log_tool_result(summary)

        # Return information about the processed data
        return {
            "message": summary,
            "columns": list(df.columns),
            "column_stats": column_stats,
            "shape": df.shape
        }
    except Exception as e:
        error_msg = f"Error processing metadata: {str(e)}"
        log(error_msg, style="bold red")
        log_tool_result(error_msg)
        return {"message": error_msg, "error": True}

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
        log_tool_header("process_metadata_with_agent")
        log(f"Preparing to process metadata at: {ctx.deps.metadata_path}",
            level=LogLevel.NORMAL)

        # Create a AnalysisContext instance specifically for the metadata agent
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
        log("Running metadata agent...", level=LogLevel.NORMAL)
        metadata_result = await metadata_agent.run(
            metadata_prompt,
            deps=metadata_deps,
            output_type=Contrasts
        )

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

            # Save contrasts to a CSV file for later use
            os.makedirs(ctx.deps.output_dir, exist_ok=True)
            contrast_path = os.path.join(ctx.deps.output_dir, "contrasts.csv")
            pd.DataFrame(contrast_data).to_csv(contrast_path, index=False)
            ctx.deps.contrast_path = contrast_path
            log(f"Saved contrasts to {contrast_path}", level=LogLevel.NORMAL)
        else:
            log("No contrasts were generated by the metadata agent", level=LogLevel.NORMAL, style="bold yellow")

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

        log_tool_result(summary)
        return summary

    except Exception as e:
        error_msg = f"Error processing metadata with agent: {str(e)}"
        log(error_msg, style="bold red")
        log_tool_result(error_msg)
        return error_msg

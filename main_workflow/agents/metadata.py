from __future__ import annotations
import logging, os, re, subprocess, json, glob, argparse, asyncio, pathlib, datetime
from typing import List, Optional, Dict, Any, Union, Tuple, Literal

import pandas as pd
import numpy as np
from dataclasses import dataclass
from dotenv import load_dotenv
from pydantic import BaseModel, Field, ConfigDict
from pydantic_ai import Agent, RunContext
from shared import AnalysisContext
from shared.workflow_logging import log_tool
from unidecode import unidecode
from openai import OpenAI
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)

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

load_dotenv()

metadata_agent = Agent(
    'openai:o4-mini',         # Use a powerful model
    deps_type=MetadataContext,
    system_prompt=system_prompt
)

# ----------------------------
# Utility function: Clean a string
# ----------------------------


@metadata_agent.tool
@log_tool
def clean_string(ctx: RunContext[MetadataContext], s: str) -> str:
    """
    Normalize and clean an input string by removing non-ASCII characters, extra whitespace, and unwanted symbols.
    """
    logger.info("🧹 Cleaning string: %s", s[:20] + "..." if len(str(s)) > 20 else s)
    if pd.isna(s):
        logger.debug("🧹 Input is NaN, returning 'NA'")
        return "NA"
    s = str(s).strip()
    s = unidecode(s)
    s = s.replace(" ", "_")
    s = re.sub(r'[^\w]', '', s)
    logger.debug("🧹 Cleaned result: %s", s)
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
        logger.info("📋 Processing metadata from %s", ctx.deps.metadata_path)

        # Load metadata based on file extension
        if ctx.deps.metadata_path.endswith('.csv'):
            df = pd.read_csv(ctx.deps.metadata_path)
        elif ctx.deps.metadata_path.endswith(('.tsv', '.txt')):
            df = pd.read_csv(ctx.deps.metadata_path, sep='\t')
        else:
            df = pd.read_csv(ctx.deps.metadata_path, sep=None, engine='python')

        # Print message of metadata dimensions
        logger.info("📊 Loaded metadata with %d rows and %d columns", df.shape[0], df.shape[1])

        # Remove duplicate rows
        initial_rows = df.shape[0]
        df = df.drop_duplicates()
        if df.shape[0] < initial_rows:
            logger.info("🧹 Removed %d duplicate rows", initial_rows - df.shape[0])

        # Remove columns where all values are identical
        initial_cols = df.shape[1]
        df = df.loc[:, df.nunique() > 1]
        if df.shape[1] < initial_cols:
            logger.info("🧹 Removed %d columns with identical values", initial_cols - df.shape[1])

        # Clean column names
        new_columns = {col: clean_string(ctx, col) for col in df.columns}
        df.rename(columns=new_columns, inplace=True)
        logger.info("✅ Cleaned column names")

        # Clean each cell value
        logger.info("🧹 Cleaning cell values...")
        for col in df.columns:
            df[col] = df[col].apply(
                lambda x: clean_string(ctx, x) if pd.notna(x) else x)

        # Store cleaned metadata in context
        ctx.deps.metadata_df = df
        logger.info("💾 Stored cleaned metadata in context")

        # Count unique values for each column
        column_stats = {}
        for col in df.columns:
            unique_vals = df[col].unique()
            column_stats[col] = {
                "unique_count": len(unique_vals),
                "values": list(unique_vals) if len(unique_vals) < 20 else list(unique_vals[:20])
            }

        summary = f"Metadata processed: {df.shape[0]} rows and {df.shape[1]} columns.\nColumns: {', '.join(df.columns)}"
        logger.info("✅ %s", summary)

        # Return information about the processed data
        return {
            "message": summary,
            "columns": list(df.columns),
            "column_stats": column_stats,
            "shape": df.shape
        }
    except Exception as e:
        error_msg = f"Error processing metadata: {str(e)}"
        logger.error("❌ %s", error_msg, exc_info=True)
        return {"message": error_msg, "error": True}

# ----------------------------
# Tool 2: Merge Analysis Columns
# ----------------------------


@metadata_agent.tool
@log_tool
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
        logger.info("🔀 Merging analysis columns from input: %s", columns_input)

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
                logger.info("📋 Parsed JSON input to: %s", candidate_cols)
            except json.JSONDecodeError:
                # Not JSON, treat as a single column name
                candidate_cols = [columns_input]
                logger.info("📋 Using input as single column name: %s", candidate_cols)

        # Case 2: List input
        elif isinstance(columns_input, list):
            candidate_cols = columns_input
            logger.info("📋 Using list input directly: %s", candidate_cols)

        # Case 3: Dict input with "columns" key
        elif isinstance(columns_input, dict) and "columns" in columns_input:
            candidate_cols = columns_input["columns"]
            logger.info("📋 Extracted columns from dict input: %s", candidate_cols)

        # Validate we have something to work with
        if not candidate_cols:
            msg = "No columns provided for merging."
            logger.warning("⚠️ %s", msg)
            return {"message": msg, "success": False, "merged_column": None}

        df = ctx.deps.metadata_df.copy()

        # Filter to ensure columns actually exist in the dataframe
        valid_cols = [col for col in candidate_cols if col in df.columns]
        if not valid_cols:
            msg = f"None of the specified columns {candidate_cols} exist in the metadata."
            logger.warning("⚠️ %s", msg)
            return {"message": msg, "success": False, "merged_column": None}

        if len(valid_cols) == 1:
            # Single column; use it directly
            ctx.deps.merged_column = valid_cols[0]
            msg = f"Single column '{valid_cols[0]}' selected as the analysis column."
            logger.info("✅ %s", msg)
        else:
            # Multiple columns; merge them
            merged_col = "merged_analysis_group"
            df[merged_col] = df[valid_cols].astype(str).apply(
                lambda row: "_".join(row.values), axis=1)
            ctx.deps.metadata_df = df  # update the DataFrame
            ctx.deps.merged_column = merged_col
            msg = f"Multiple columns {', '.join(valid_cols)} merged into column '{merged_col}'."
            logger.info("✅ %s", msg)

            # Show a preview of the merged values
            unique_merged = df[merged_col].unique().tolist()
            preview = unique_merged[:5] if len(unique_merged) > 5 else unique_merged
            logger.info("📊 Merged values preview: %s", preview)
            msg += f"\nMerged values (preview): {preview}"

        return {
            "message": msg,
            "success": True,
            "merged_column": ctx.deps.merged_column,
            "input_columns": valid_cols
        }
    except Exception as e:
        error_msg = f"Error in merge_analysis_columns: {str(e)}"
        logger.error("❌ %s", error_msg, exc_info=True)
        return {"message": error_msg, "success": False, "merged_column": None}



# ----------------------------
# Tool 3: Extract Unique Values
# ----------------------------


@metadata_agent.tool
@log_tool
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
        logger.info("🔍 Extracting unique values from analysis column")

        if ctx.deps.metadata_df is None:
            error_msg = "Error: Metadata has not been processed."
            logger.error("❌ %s", error_msg)
            return {"success": False, "message": error_msg, "unique_values": []}

        if not ctx.deps.merged_column:
            error_msg = "Error: Analysis column not defined. Please run merge_analysis_columns first."
            logger.error("❌ %s", error_msg)
            return {"success": False, "message": error_msg, "unique_values": []}

        df = ctx.deps.metadata_df
        analysis_col = ctx.deps.merged_column
        logger.info("📊 Using analysis column: %s", analysis_col)

        # Extract unique values from the analysis column
        unique_values = sorted(df[analysis_col].dropna().unique().tolist())
        ctx.deps.unique_groups = unique_values
        logger.info("✅ Extracted %d unique values: %s", len(unique_values), unique_values)

        # Return only the unique values with basic metadata
        return {
            "success": True,
            "column": analysis_col,
            "unique_values": unique_values,
            "count": len(unique_values)
        }

    except Exception as e:
        error_msg = f"Error extracting unique values: {str(e)}"
        logger.error("❌ %s", error_msg, exc_info=True)
        return {"success": False, "message": error_msg, "unique_values": []}

@log_tool
async def run_agent_async(prompt: str, deps: MetadataContext, usage=None):
    logger.info("📋 Metadata agent invoked – prompt: %s", prompt)
    result = await metadata_agent.run(prompt, deps=deps, usage=usage, output_type=output_type)

    # Log the agent's output
    logger.info("📄 Metadata agent output: %s", result.output)

    # Log any contrasts generated (if available)
    if hasattr(result, 'output') and hasattr(result.output, 'contrasts'):
        logger.info("📊 Metadata agent generated %d contrasts", len(result.output.contrasts))

    # If you want to log usage statistics
    if hasattr(result, 'usage') and result.usage:
        try:
            usage_stats = result.usage()
            logger.info("📊 Metadata agent usage: %s", usage_stats)
        except Exception as e:
            logger.debug("Could not get usage stats: %s", e)

    return result

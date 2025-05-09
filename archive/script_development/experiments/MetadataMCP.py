import argparse
import os
import re
import json
import pandas as pd
from unidecode import unidecode
from dataclasses import dataclass
from typing import List, Dict, Any

from pydantic_ai import Agent, RunContext
from rich.console import Console

console = Console()

@dataclass
class MetadataContext:
    """
    Holds runtime state for metadata processing.
    After process_metadata runs, metadata_df is a DataFrame.
    After extract_unique_values, unique_groups is populated.
    """
    metadata_path: str
    metadata_df: pd.DataFrame = None
    merged_column: str = None
    unique_groups: List[str] = None

# Initialize the agent
metadata_agent = Agent(
    'openai:gpt-4o',
    deps_type=MetadataContext,
    system_prompt="""
You are an expert RNA‐seq metadata processor. Use the provided tools
to clean strings, load/clean metadata, and extract unique group values.
""".strip()
)

@metadata_agent.tool
def clean_string(ctx: RunContext[MetadataContext], s: Any) -> str:
    """
    Normalize and clean an input string by:
      - Converting non-ASCII to ASCII
      - Stripping whitespace
      - Replacing spaces with underscores
      - Removing non-word characters
    """
    if pd.isna(s):
        return "NA"
    txt = str(s).strip()
    txt = unidecode(txt)
    txt = re.sub(r'\s+', '_', txt)
    txt = re.sub(r'[^\w]', '', txt)
    return txt

@metadata_agent.tool
def process_metadata(ctx: RunContext[MetadataContext]) -> dict:
    """
    Load and clean metadata:
      1. Read CSV/TSV (auto-detect by extension)
      2. Drop columns with no variation
      3. Clean column names & cell values via clean_string
      4. Store DataFrame in ctx.deps.metadata_df
    Returns a dict with summary, column list, per-column stats, and shape.
    """
    path = ctx.deps.metadata_path
    if not os.path.exists(path):
        return {"message": f"Error: file not found at {path}", "columns": [], "column_stats": {}, "shape": []}
    # load
    if path.lower().endswith(('.tsv', '.txt')):
        df = pd.read_csv(path, sep='\t')
    else:
        df = pd.read_csv(path)
    # drop constant columns
    df = df.loc[:, df.nunique() > 1]
    # clean column names
    rename_map = {col: clean_string(ctx, col) for col in df.columns}
    df.rename(columns=rename_map, inplace=True)
    # clean cell values
    for col in df.columns:
        df[col] = df[col].apply(lambda x: clean_string(ctx, x) if pd.notna(x) else x)
    # store
    ctx.deps.metadata_df = df
    # build stats
    column_stats: Dict[str, Dict[str, Any]] = {}
    for col in df.columns:
        vals = df[col].dropna().unique().tolist()
        column_stats[col] = {
            "unique_count": len(vals),
            "values": vals if len(vals) < 20 else vals[:20]
        }
    return {
        "message": f"Metadata processed: {df.shape[0]} rows, {df.shape[1]} columns",
        "columns": list(df.columns),
        "column_stats": column_stats,
        "shape": list(df.shape)
    }

@metadata_agent.tool
def extract_unique_values(ctx: RunContext[MetadataContext]) -> dict:
    """
    Extract unique values from ctx.deps.merged_column of the cleaned metadata.
    Populates ctx.deps.unique_groups.
    """
    df = ctx.deps.metadata_df
    col = ctx.deps.merged_column
    if df is None:
        return {"success": False, "message": "Error: metadata not yet processed", "unique_values": [], "count": 0}
    if not col or col not in df.columns:
        return {"success": False, "message": f"Error: merged_column '{col}' invalid", "unique_values": [], "count": 0}
    vals = sorted(df[col].dropna().unique().tolist())
    ctx.deps.unique_groups = vals
    return {
        "success": True,
        "message": f"Extracted {len(vals)} unique values from '{col}'",
        "unique_values": vals,
        "count": len(vals)
    }

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="MCP Test: Metadata Agent")
    parser.add_argument("--metadata_path", required=True,
                        help="Path to metadata file (CSV/TSV)")
    parser.add_argument("--merged_column", required=False,
                        help="Name of column to extract unique values from")
    args = parser.parse_args()

    # Create deps and run
    deps = MetadataContext(metadata_path=args.metadata_path)
    # 1. Clean & load metadata
    meta_out = metadata_agent.run_sync("process_metadata", deps=deps)
    console.print(meta_out)
    # 2. Optionally extract unique values
    if args.merged_column:
        deps.merged_column = args.merged_column
        uniq_out = metadata_agent.run_sync("extract_unique_values", deps=deps)
        console.print(uniq_out)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MCP Data Extraction Server for RNA-seq Analysis

This server provides low-level data extraction tools for accessing and querying
RNA-seq differential expression results. It focuses on raw data extraction from
CSV files and basic filtering operations.

Tools provided:
- get_deg_counts: Count differentially expressed genes for a contrast
- get_gene_info: Get expression data for specific genes
- find_common_degs: Find genes differentially expressed across multiple contrasts
- get_contrast_metadata: Get metadata about contrasts
- list_available_contrasts: List all available contrasts in the dataset
"""

import os
import sys
import json
import pandas as pd
import numpy as np
from typing import List, Dict, Optional, Tuple, Any
from pathlib import Path
import logging
from collections import defaultdict
import argparse

from mcp.server.fastmcp import FastMCP

# Initialize MCP server
mcp = FastMCP("data_extractor")
log = lambda m: print(f"[DATA_EXTRACTOR] {m}", file=sys.stderr, flush=True)

# Set up logging
logging.basicConfig(
    level=logging.WARNING,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    stream=sys.stderr
)
logger = logging.getLogger("mcp-data-extractor")

# Add a decorator to log tool execution times
def log_tool_execution(func):
    async def wrapper(*args, **kwargs):
        import time
        logger.info(f"Starting tool: {func.__name__}")
        start_time = time.time()
        try:
            result = await func(*args, **kwargs)
            logger.info(f"Tool {func.__name__} completed in {time.time() - start_time:.2f}s")
            return result
        except Exception as e:
            logger.error(f"Tool {func.__name__} failed after {time.time() - start_time:.2f}s: {e}")
            raise
    return wrapper

# Global variables
RESULTS_DIR = os.environ.get("RESULTS_DIR", "")
_cache = {}  # Simple cache for repeated queries
_contrast_cache = {}  # Cache for contrast metadata

def _get_deg_file_path(analysis_id: str, contrast_id: str) -> Optional[str]:
    """Find the DEG file path for a given analysis and contrast."""
    # Try multiple possible paths
    possible_paths = [
        os.path.join(RESULTS_DIR, analysis_id, "RNAseqAnalysis", contrast_id, "DEG.csv"),
        os.path.join(RESULTS_DIR, "RNAseqAnalysis", contrast_id, "DEG.csv"),
    ]
    
    for path in possible_paths:
        if os.path.exists(path):
            return path
    
    # If not found by exact match, search recursively
    for root, dirs, files in os.walk(RESULTS_DIR):
        if "DEG.csv" in files and contrast_id in root:
            return os.path.join(root, "DEG.csv")
    
    return None

def _load_deg_data(analysis_id: str, contrast_id: str) -> Optional[pd.DataFrame]:
    """Load DEG data with caching."""
    cache_key = f"deg_{analysis_id}_{contrast_id}"
    if cache_key in _cache:
        return _cache[cache_key]
    
    file_path = _get_deg_file_path(analysis_id, contrast_id)
    if not file_path:
        logger.warning(f"No DEG file found for {analysis_id}/{contrast_id}")
        return None
    
    try:
        df = pd.read_csv(file_path)
        # Standardize gene column name
        gene_col = None
        for col in ['Gene', 'gene', 'SYMBOL', 'symbol', 'gene_name', 'gene_id']:
            if col in df.columns:
                df = df.rename(columns={col: 'Gene'})
                break
        
        _cache[cache_key] = df
        return df
    except Exception as e:
        logger.error(f"Error loading DEG file {file_path}: {e}")
        return None

def _load_contrast_metadata() -> Dict[str, Dict[str, Any]]:
    """Load contrast metadata from various possible locations."""
    if _contrast_cache:
        return _contrast_cache
    
    contrast_info = {}
    
    # Look for contrasts.csv files
    possible_locations = [
        os.path.join(RESULTS_DIR, "metadata", "contrasts.csv"),
        os.path.join(RESULTS_DIR, "contrasts.csv"),
    ]
    
    # Also search in analysis subdirectories
    for root, dirs, files in os.walk(RESULTS_DIR):
        if "contrasts.csv" in files:
            possible_locations.append(os.path.join(root, "contrasts.csv"))
    
    for contrasts_file in possible_locations:
        if os.path.exists(contrasts_file):
            try:
                contrasts_df = pd.read_csv(contrasts_file)
                for _, row in contrasts_df.iterrows():
                    if 'name' in row:
                        contrast_key = row['name']
                        contrast_info[contrast_key] = {
                            'name': row['name'],
                            'description': row.get('description', ''),
                            'expression': row.get('expression', ''),
                            'formula': row.get('formula', row.get('expression', ''))
                        }
                logger.info(f"Loaded {len(contrasts_df)} contrasts from {contrasts_file}")
            except Exception as e:
                logger.warning(f"Error loading contrasts from {contrasts_file}: {e}")
    
    _contrast_cache.update(contrast_info)
    return _contrast_cache

@mcp.tool()
@log_tool_execution
async def list_available_contrasts() -> str:
    """
    List all available contrasts in the dataset.
    
    Returns:
        JSON string with list of available contrasts and their metadata.
    """
    try:
        contrasts = []
        
        # Find all DEG files
        for root, dirs, files in os.walk(RESULTS_DIR):
            if "DEG.csv" in files and "RNAseqAnalysis" in root:
                # Extract analysis and contrast IDs from path
                path_parts = root.split(os.sep)
                
                # Find the analysis ID (parent of RNAseqAnalysis)
                analysis_id = None
                contrast_id = None
                
                for i, part in enumerate(path_parts):
                    if part == "RNAseqAnalysis" and i > 0:
                        analysis_id = path_parts[i-1]
                        if i + 1 < len(path_parts):
                            contrast_id = path_parts[i+1]
                        break
                
                if not analysis_id:
                    # Try to extract from the path differently
                    if "RNAseqAnalysis" in path_parts:
                        idx = path_parts.index("RNAseqAnalysis")
                        if idx > 0:
                            analysis_id = path_parts[idx-1]
                        contrast_id = os.path.basename(root)
                
                if analysis_id and contrast_id:
                    # Count DEGs
                    deg_count = 0
                    try:
                        df = pd.read_csv(os.path.join(root, "DEG.csv"))
                        if 'adj.P.Val' in df.columns and 'logFC' in df.columns:
                            deg_count = ((df['adj.P.Val'] < 0.05) & (abs(df['logFC']) > 1.0)).sum()
                    except:
                        pass
                    
                    contrasts.append({
                        'analysis_id': analysis_id,
                        'contrast_id': contrast_id,
                        'deg_count': int(deg_count),
                        'path': root
                    })
        
        # Load contrast metadata
        contrast_metadata = _load_contrast_metadata()
        
        # Enhance contrast info with metadata
        for contrast in contrasts:
            contrast_id = contrast['contrast_id']
            if contrast_id in contrast_metadata:
                contrast['description'] = contrast_metadata[contrast_id].get('description', '')
                contrast['formula'] = contrast_metadata[contrast_id].get('formula', '')
        
        return json.dumps({
            'contrasts': contrasts,
            'total_count': len(contrasts)
        })
    
    except Exception as e:
        logger.error(f"Error in list_available_contrasts: {e}")
        return json.dumps({"error": str(e)})

@mcp.tool()
@log_tool_execution
async def get_deg_counts(analysis_id: str,
                        contrast_id: str,
                        p_value_threshold: float = 0.05,
                        lfc_threshold: float = 1.0) -> str:
    """
    Get the number of differentially expressed genes for a specific contrast.
    
    Args:
        analysis_id: Identifier for the analysis/dataset.
        contrast_id: Identifier for the contrast.
        p_value_threshold: Adjusted p-value threshold for significance.
        lfc_threshold: Log fold change threshold for significance.
    
    Returns:
        JSON string with count statistics.
    """
    try:
        df = _load_deg_data(analysis_id, contrast_id)
        if df is None:
            return json.dumps({"error": f"No DEG data found for {analysis_id}/{contrast_id}"})
        
        # Identify p-value and logFC columns
        p_value_col = None
        if 'adj.P.Val' in df.columns:
            p_value_col = 'adj.P.Val'
        elif 'padj' in df.columns:
            p_value_col = 'padj'
        elif 'FDR' in df.columns:
            p_value_col = 'FDR'
        elif 'P.Value' in df.columns:
            p_value_col = 'P.Value'
        
        lfc_col = None
        if 'logFC' in df.columns:
            lfc_col = 'logFC'
        elif 'log2FoldChange' in df.columns:
            lfc_col = 'log2FoldChange'
        
        if not p_value_col or not lfc_col:
            return json.dumps({
                "error": "Missing required columns",
                "available_columns": df.columns.tolist()
            })
        
        # Count DEGs
        sig_mask = (df[p_value_col] < p_value_threshold) & (abs(df[lfc_col]) > lfc_threshold)
        up_mask = sig_mask & (df[lfc_col] > 0)
        down_mask = sig_mask & (df[lfc_col] < 0)
        
        total_count = sig_mask.sum()
        up_count = up_mask.sum()
        down_count = down_mask.sum()
        
        # Get contrast description if available
        contrast_metadata = _load_contrast_metadata()
        description = contrast_metadata.get(contrast_id, {}).get('description', '')
        
        return json.dumps({
            "analysis_id": analysis_id,
            "contrast_id": contrast_id,
            "description": description,
            "total_count": int(total_count),
            "up_regulated": int(up_count),
            "down_regulated": int(down_count),
            "total_genes_tested": len(df),
            "thresholds": {
                "p_value": p_value_threshold,
                "lfc": lfc_threshold
            }
        })
    
    except Exception as e:
        logger.error(f"Error in get_deg_counts: {e}")
        return json.dumps({"error": str(e)})

@mcp.tool()
@log_tool_execution
async def get_gene_info(gene_id: str,
                       analysis_ids: Optional[List[str]] = None,
                       contrast_ids: Optional[List[str]] = None) -> str:
    """
    Get expression information for a specific gene across selected analyses and contrasts.
    
    Args:
        gene_id: Gene identifier/symbol.
        analysis_ids: List of analysis IDs to include (None = all).
        contrast_ids: List of contrast IDs to include (None = all).
    
    Returns:
        JSON string with gene expression data across contrasts.
    """
    try:
        gene_data = []
        
        # Get all available contrasts if not specified
        if not analysis_ids or not contrast_ids:
            contrasts_response = await list_available_contrasts()
            all_contrasts = json.loads(contrasts_response).get('contrasts', [])
            
            if not analysis_ids:
                analysis_ids = list(set(c['analysis_id'] for c in all_contrasts))
            if not contrast_ids:
                contrast_ids = list(set(c['contrast_id'] for c in all_contrasts))
        
        # Collect gene data across contrasts
        for analysis_id in analysis_ids:
            for contrast_id in contrast_ids:
                df = _load_deg_data(analysis_id, contrast_id)
                if df is None or 'Gene' not in df.columns:
                    continue
                
                # Find the gene
                gene_row = df[df['Gene'] == gene_id]
                if gene_row.empty:
                    continue
                
                row = gene_row.iloc[0]
                
                # Extract relevant columns
                data_point = {
                    'analysis_id': analysis_id,
                    'contrast_id': contrast_id,
                    'gene': gene_id
                }
                
                # Add expression values
                if 'logFC' in df.columns:
                    data_point['logFC'] = float(row['logFC'])
                if 'adj.P.Val' in df.columns:
                    data_point['adj_p_value'] = float(row['adj.P.Val'])
                elif 'P.Value' in df.columns:
                    data_point['p_value'] = float(row['P.Value'])
                
                # Add other statistics if available
                if 'AveExpr' in df.columns:
                    data_point['average_expression'] = float(row['AveExpr'])
                if 't' in df.columns:
                    data_point['t_statistic'] = float(row['t'])
                
                gene_data.append(data_point)
        
        if not gene_data:
            return json.dumps({
                "gene": gene_id,
                "message": "Gene not found in any of the specified contrasts",
                "contrasts_checked": len(analysis_ids) * len(contrast_ids)
            })
        
        return json.dumps({
            "gene": gene_id,
            "data": gene_data,
            "contrasts_found": len(gene_data)
        })
    
    except Exception as e:
        logger.error(f"Error in get_gene_info: {e}")
        return json.dumps({"error": str(e)})

@mcp.tool()
@log_tool_execution
async def find_common_degs(analysis_ids: List[str],
                          contrast_ids: List[str],
                          min_frequency: int = 2,
                          p_value_threshold: float = 0.05,
                          lfc_threshold: float = 1.0,
                          max_genes: int = 50) -> str:
    """
    Find genes that are differentially expressed across multiple contrasts.
    
    Args:
        analysis_ids: List of analysis identifiers.
        contrast_ids: List of contrast identifiers.
        min_frequency: Minimum number of contrasts a gene must appear in.
        p_value_threshold: Adjusted p-value threshold for significance.
        lfc_threshold: Log fold change threshold for significance.
        max_genes: Maximum number of genes to return.
    
    Returns:
        JSON string with common DEGs and their properties.
    """
    try:
        # Track genes across contrasts
        gene_occurrences = defaultdict(int)
        gene_data = defaultdict(lambda: {
            "contrasts": {},
            "max_abs_lfc": 0,
            "consistent_direction": True,
            "directions": []
        })
        
        contrasts_analyzed = 0
        
        # Analyze each contrast
        for analysis_id in analysis_ids:
            for contrast_id in contrast_ids:
                df = _load_deg_data(analysis_id, contrast_id)
                if df is None:
                    continue
                
                contrasts_analyzed += 1
                
                # Identify columns
                p_value_col = None
                if 'adj.P.Val' in df.columns:
                    p_value_col = 'adj.P.Val'
                elif 'P.Value' in df.columns:
                    p_value_col = 'P.Value'
                
                lfc_col = 'logFC' if 'logFC' in df.columns else None
                
                if not all([p_value_col, lfc_col, 'Gene' in df.columns]):
                    continue
                
                # Filter for significant genes
                sig_df = df[(df[p_value_col] < p_value_threshold) & 
                           (abs(df[lfc_col]) > lfc_threshold)]
                
                # Count occurrences
                for _, row in sig_df.iterrows():
                    gene = row['Gene']
                    contrast_key = f"{analysis_id}:{contrast_id}"
                    
                    gene_occurrences[gene] += 1
                    gene_data[gene]["contrasts"][contrast_key] = {
                        "lfc": float(row[lfc_col]),
                        "p_value": float(row[p_value_col])
                    }
                    gene_data[gene]["max_abs_lfc"] = max(
                        gene_data[gene]["max_abs_lfc"], 
                        abs(row[lfc_col])
                    )
                    gene_data[gene]["directions"].append(np.sign(row[lfc_col]))
        
        # Check direction consistency
        for gene, data in gene_data.items():
            directions = data["directions"]
            if len(set(directions)) > 1:
                data["consistent_direction"] = False
        
        # Filter by frequency
        common_genes = []
        for gene, frequency in gene_occurrences.items():
            if frequency >= min_frequency:
                common_genes.append({
                    "gene": gene,
                    "frequency": frequency,
                    "percentage": round(frequency / contrasts_analyzed * 100, 1),
                    "contrasts": gene_data[gene]["contrasts"],
                    "max_abs_lfc": round(gene_data[gene]["max_abs_lfc"], 3),
                    "consistent_direction": gene_data[gene]["consistent_direction"],
                    "mean_lfc": round(
                        np.mean([c["lfc"] for c in gene_data[gene]["contrasts"].values()]), 
                        3
                    )
                })
        
        # Sort by frequency and max LFC
        common_genes.sort(key=lambda x: (x["frequency"], x["max_abs_lfc"]), reverse=True)
        
        # Limit number of returned genes
        common_genes = common_genes[:max_genes]
        
        return json.dumps({
            "genes": common_genes,
            "total_found": len(common_genes),
            "contrasts_analyzed": contrasts_analyzed,
            "thresholds": {
                "min_frequency": min_frequency,
                "p_value": p_value_threshold,
                "lfc": lfc_threshold
            }
        })
    
    except Exception as e:
        logger.error(f"Error in find_common_degs: {e}")
        return json.dumps({"error": str(e)})

@mcp.tool()
@log_tool_execution
async def get_contrast_metadata(analysis_id: str, contrast_id: str) -> str:
    """
    Get metadata information about a specific contrast.
    
    Args:
        analysis_id: Analysis identifier.
        contrast_id: Contrast identifier.
    
    Returns:
        JSON string with contrast metadata.
    """
    try:
        # Load contrast metadata
        contrast_metadata = _load_contrast_metadata()
        
        # Get basic info
        metadata = {
            "analysis_id": analysis_id,
            "contrast_id": contrast_id
        }
        
        # Add contrast description if available
        if contrast_id in contrast_metadata:
            metadata.update(contrast_metadata[contrast_id])
        
        # Try to get DEG counts
        deg_response = await get_deg_counts(analysis_id, contrast_id)
        deg_data = json.loads(deg_response)
        if "error" not in deg_data:
            metadata["deg_summary"] = {
                "total_degs": deg_data["total_count"],
                "up_regulated": deg_data["up_regulated"],
                "down_regulated": deg_data["down_regulated"],
                "total_genes": deg_data["total_genes_tested"]
            }
        
        # Check if files exist
        deg_path = _get_deg_file_path(analysis_id, contrast_id)
        metadata["files"] = {
            "deg_file": deg_path is not None,
            "deg_path": deg_path if deg_path else "Not found"
        }
        
        # Look for additional files
        if deg_path:
            contrast_dir = os.path.dirname(deg_path)
            additional_files = []
            for file in ["volcano_plot.png", "ma_plot.png", "heatmap_top50.png"]:
                if os.path.exists(os.path.join(contrast_dir, file)):
                    additional_files.append(file)
            metadata["files"]["plots"] = additional_files
        
        return json.dumps(metadata)
    
    except Exception as e:
        logger.error(f"Error in get_contrast_metadata: {e}")
        return json.dumps({"error": str(e)})

def main():
    """Main entry point for the MCP server."""
    parser = argparse.ArgumentParser(description="MCP Data Extractor Server")
    parser.add_argument("command", choices=["server"], help="Run as MCP server")
    parser.add_argument("--results-dir", help="Path to results directory")
    args = parser.parse_args()

    global RESULTS_DIR
    if args.results_dir:
        RESULTS_DIR = args.results_dir
        log(f"Using results directory from command line: {RESULTS_DIR}")

    logger.info(f"Starting data extractor MCP server with results_dir: {RESULTS_DIR}")
    if not RESULTS_DIR:
        logger.warning("RESULTS_DIR environment variable not set")
    elif not os.path.exists(RESULTS_DIR):
        logger.warning(f"RESULTS_DIR does not exist: {RESULTS_DIR}")
    
    # Run the MCP server
    mcp.run(transport="stdio")

if __name__ == "__main__":
    main()
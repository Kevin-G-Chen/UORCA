#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
UORCA Data MCP Server
--------------------
A Model Context Protocol server that provides access to UORCA RNA-seq analysis results.
Exposes tools for querying contrasts, genes, and fold changes from analysis data.
"""

import os
import sys
import json
import signal
import logging
import argparse
from pathlib import Path
from typing import Any, Dict, List, Optional, Set

import pandas as pd
from mcp.server.fastmcp import FastMCP

# Add the reporting directory to the path so we can import ResultsIntegrator
sys.path.insert(0, str(Path(__file__).parent.parent))
from ResultsIntegration import ResultsIntegrator

# -----------------------------------------------------------------------------
# MCP server setup & logging
# -----------------------------------------------------------------------------
mcp = FastMCP("uorca_data")

# Setup logging to stderr
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s – %(name)s – %(levelname)s – %(message)s",
    stream=sys.stderr
)
logger = logging.getLogger("uorca-data-server")

def log(message: str):
    """Simple logging helper"""
    print(f"UORCA MCP: {message}", file=sys.stderr, flush=True)

# -----------------------------------------------------------------------------
# Data Access Layer
# -----------------------------------------------------------------------------
class UORCADataAccess:
    """
    Provides access to UORCA analysis results data.
    """

    def __init__(self, results_dir: Optional[str] = None):
        """Initialize with results directory."""
        self.results_dir = results_dir or os.environ.get("UORCA_RESULTS_DIR", "/UORCA_results")
        self.integrator = None
        self._initialize_integrator()

    def _initialize_integrator(self):
        """Initialize the ResultsIntegrator."""
        try:
            if not os.path.exists(self.results_dir):
                raise FileNotFoundError(f"Results directory not found: {self.results_dir}")

            self.integrator = ResultsIntegrator(self.results_dir)
            self.integrator.load_analysis_info()
            self.integrator.load_data()
            log(f"Successfully initialized UORCA data access for: {self.results_dir}")

        except Exception as e:
            log(f"Failed to initialize UORCA data access: {e}")
            self.integrator = None

    def get_contrasts(self) -> List[str]:
        """Get list of all available contrasts."""
        if not self.integrator:
            return []

        contrasts = []
        for analysis_id, analysis_contrasts in self.integrator.deg_data.items():
            for contrast_id in analysis_contrasts.keys():
                # Create a full contrast identifier
                full_contrast = f"{analysis_id}:{contrast_id}"
                contrasts.append(full_contrast)

        return sorted(contrasts)

    def get_genes(self, limit: int = 1000) -> List[str]:
        """Get list of genes across all analyses (limited for performance)."""
        if not self.integrator:
            return []

        genes = set()
        count = 0

        # Collect genes from DEG data
        for analysis_id, analysis_contrasts in self.integrator.deg_data.items():
            for contrast_id, deg_df in analysis_contrasts.items():
                if 'Gene' in deg_df.columns:
                    for gene in deg_df['Gene'].dropna().unique():
                        genes.add(str(gene))
                        count += 1
                        if count >= limit:
                            break
                if count >= limit:
                    break
            if count >= limit:
                break

        return sorted(list(genes))

    def get_log_fold_change(self, contrast: str, gene: str) -> Optional[float]:
        """Get log fold change for a specific gene in a specific contrast."""
        if not self.integrator:
            return None

        try:
            # Parse the contrast identifier
            if ":" not in contrast:
                # If no analysis ID provided, search all analyses
                for analysis_id, analysis_contrasts in self.integrator.deg_data.items():
                    if contrast in analysis_contrasts:
                        deg_df = analysis_contrasts[contrast]
                        break
                else:
                    return None
            else:
                analysis_id, contrast_id = contrast.split(":", 1)
                if analysis_id not in self.integrator.deg_data:
                    return None
                if contrast_id not in self.integrator.deg_data[analysis_id]:
                    return None
                deg_df = self.integrator.deg_data[analysis_id][contrast_id]

            # Find the gene in the dataframe
            if 'Gene' not in deg_df.columns:
                return None

            gene_row = deg_df[deg_df['Gene'] == gene]
            if gene_row.empty:
                return None

            # Look for log fold change column (common names)
            lfc_columns = ['logFC', 'log2FoldChange', 'log2FC', 'LFC', 'LogFC', 'logfoldchange']
            lfc_col = None

            for col in lfc_columns:
                if col in deg_df.columns:
                    lfc_col = col
                    break

            if lfc_col is None:
                return None

            lfc_value = gene_row[lfc_col].iloc[0]
            return float(lfc_value) if not pd.isna(lfc_value) else None

        except Exception as e:
            log(f"Error getting LFC for {gene} in {contrast}: {e}")
            return None

    def get_analysis_info(self) -> Dict[str, Any]:
        """Get information about available analyses."""
        if not self.integrator:
            return {}

        info = {
            "analyses": list(self.integrator.analysis_info.keys()),
            "total_contrasts": len(self.get_contrasts()),
            "results_dir": self.results_dir
        }

        # Add analysis details
        info["analysis_details"] = {}
        for analysis_id, analysis_data in self.integrator.analysis_info.items():
            info["analysis_details"][analysis_id] = {
                "contrasts": list(self.integrator.deg_data.get(analysis_id, {}).keys()),
                "metadata": analysis_data
            }

        return info

# Global data access instance
_data_access: Optional[UORCADataAccess] = None

def get_data_access() -> UORCADataAccess:
    """Get or create the data access instance."""
    global _data_access
    if _data_access is None:
        results_dir = os.environ.get("UORCA_RESULTS_DIR")
        _data_access = UORCADataAccess(results_dir)
    return _data_access

# -----------------------------------------------------------------------------
# MCP Tools
# -----------------------------------------------------------------------------

@mcp.tool()
async def list_contrasts() -> str:
    """
    List all available contrasts across all UORCA analyses.

    Returns:
        JSON string containing list of contrast identifiers in format "analysis_id:contrast_id"
    """
    try:
        data_access = get_data_access()
        contrasts = data_access.get_contrasts()
        return json.dumps({"contrasts": contrasts}, ensure_ascii=False)
    except Exception as e:
        log(f"ERROR list_contrasts: {e}")
        return json.dumps({"contrasts": [], "error": str(e)}, ensure_ascii=False)

@mcp.tool()
async def list_genes(limit: int = 1000) -> str:
    """
    List genes available in the UORCA analyses.

    Args:
        limit: Maximum number of genes to return (default: 1000)

    Returns:
        JSON string containing list of gene symbols
    """
    try:
        data_access = get_data_access()
        genes = data_access.get_genes(limit=limit)
        return json.dumps({"genes": genes[:limit]}, ensure_ascii=False)
    except Exception as e:
        log(f"ERROR list_genes: {e}")
        return json.dumps({"genes": [], "error": str(e)}, ensure_ascii=False)

@mcp.tool()
async def get_lfc(contrast: str, gene: str) -> str:
    """
    Get log fold change for a specific gene in a specific contrast.

    Args:
        contrast: Contrast identifier (format: "analysis_id:contrast_id")
        gene: Gene symbol

    Returns:
        JSON string with log fold change value or null if not found
    """
    try:
        data_access = get_data_access()
        lfc = data_access.get_log_fold_change(contrast, gene)

        return json.dumps({
            "gene": gene,
            "contrast": contrast,
            "log_fold_change": lfc
        }, ensure_ascii=False)
    except Exception as e:
        log(f"ERROR get_lfc: {e}")
        return json.dumps({
            "gene": gene,
            "contrast": contrast,
            "log_fold_change": None,
            "error": str(e)
        }, ensure_ascii=False)

@mcp.tool()
async def get_analysis_info() -> str:
    """
    Get information about available UORCA analyses.

    Returns:
        JSON string with analysis metadata and structure
    """
    try:
        data_access = get_data_access()
        info = data_access.get_analysis_info()
        return json.dumps(info, ensure_ascii=False)
    except Exception as e:
        log(f"ERROR get_analysis_info: {e}")
        return json.dumps({"error": str(e)}, ensure_ascii=False)

# -----------------------------------------------------------------------------
# Signal handling and CLI
# -----------------------------------------------------------------------------

def _signal_handler(sig, frame):
    log(f"Received signal {sig} - exiting")
    sys.exit(0)

signal.signal(signal.SIGINT, _signal_handler)
signal.signal(signal.SIGTERM, _signal_handler)

def main():
    parser = argparse.ArgumentParser(description="UORCA Data MCP Server")
    parser.add_argument(
        "--transport", default="stdio", choices=["stdio", "sse"],
        help="Transport protocol (stdio or server-sent events)"
    )
    parser.add_argument(
        "--debug", action="store_true", help="Enable debug logging"
    )
    parser.add_argument(
        "--results-dir", help="Path to UORCA results directory"
    )
    args = parser.parse_args()

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)

    # Set results directory if provided
    if args.results_dir:
        os.environ["UORCA_RESULTS_DIR"] = args.results_dir
        log(f"Using results directory: {args.results_dir}")

    # Initialize data access to validate setup
    try:
        data_access = get_data_access()
        if data_access.integrator is None:
            log("WARNING: Failed to initialize data access - tools may not work properly")
        else:
            contrasts = data_access.get_contrasts()
            log(f"Successfully initialized with {len(contrasts)} contrasts")
    except Exception as e:
        log(f"ERROR during initialization: {e}")

    log("Starting UORCA Data MCP Server...")
    mcp.run(transport=args.transport)

if __name__ == "__main__":
    main()

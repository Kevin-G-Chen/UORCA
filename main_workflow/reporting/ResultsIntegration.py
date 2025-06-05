#!/usr/bin/env python3
"""
ResultsIntegration.py - Integrates results from multiple UORCA RNA-seq analyses

This script collates differential expression results and normalized expression data
across experiments and contrasts, producing interactive visualizations
that highlight important genes and patterns.

The script performs three main functions:
1. Identifying differentially expressed genes (DEGs) across contrasts
   - Finds genes differentially expressed in multiple contrasts
   - Identifies contrast-specific genes with large fold changes

2. Creating an interactive heatmap of log fold changes
   - Shows gene expression patterns across contrasts
   - Clusters genes and contrasts to reveal patterns

3. Generating expression plots for selected genes
   - Shows expression levels across samples grouped by condition
   - Uses sample information from edger_analysis_samples.csv

The main outputs are:
1. Interactive heatmap of important genes with LFC values (clustered on both axes)
2. Interactive boxplot/violin plots of gene expression across samples
3. HTML report integrating all visualizations

Usage:
  python ResultsIntegration.py --results_dir /path/to/results --output_dir /path/to/output
"""

import os
import glob
import argparse
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import pdist
from typing import List, Dict, Tuple, Optional, Union, Set
import logging
import re
from pathlib import Path
from datetime import datetime
import json
import math

# Set up logging
logging.basicConfig(
    level=logging.WARNING,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class ResultsIntegrator:
    """
    Integrates and visualizes results from multiple RNA-seq analyses.
    """

    def __init__(self,
                 results_dir: str,
                 output_dir: str = None,
                 pvalue_threshold: float = 0.05,
                 lfc_threshold: float = 1.0,
                 top_n_genes: int = 50):
        """
        Initialize the ResultsIntegrator.

        Parameters:
        -----------
        results_dir : str
            Directory containing analysis results
        output_dir : str, optional
            Directory where integrated results will be saved
        pvalue_threshold : float
            P-value threshold for considering genes as differentially expressed
        lfc_threshold : float
            Log fold change threshold for considering genes as differentially expressed
        top_n_genes : int
            Number of top genes to include in visualizations by default
        """
        self.results_dir = os.path.abspath(results_dir)
        if not os.path.isdir(self.results_dir):
            raise FileNotFoundError(f"Results directory not found: {self.results_dir}")

        # Store requested output directory but don't create it yet
        self._output_dir_requested = output_dir
        self.output_dir = None

        # Set threshold parameters
        self.pvalue_threshold = pvalue_threshold
        self.lfc_threshold = lfc_threshold
        self.top_n_genes = top_n_genes

        # Initialize containers for data
        self.analysis_info = {}   # Analysis metadata
        self.deg_data = {}        # DEG data for each contrast
        self.cpm_data = {}        # CPM data for each analysis
        self.contrast_info = {}   # Information about contrasts

        logger.info(f"Initialized ResultsIntegrator with results_dir: {self.results_dir}")

    def find_analysis_folders(self) -> List[str]:
        """Find all analysis folders in the results directory."""
        # Look for directories that contain RNAseqAnalysis folders
        analysis_folders = []

        # Check if the results directory itself contains RNAseqAnalysis
        rnaseq_dir = os.path.join(self.results_dir, "RNAseqAnalysis")
        if os.path.isdir(rnaseq_dir):
            analysis_folders.append(self.results_dir)
            logger.info(f"Found analysis folder: {self.results_dir}")

        # If no analysis folders found in the current directory, look for subdirectories
        if not analysis_folders:
            for item in os.listdir(self.results_dir):
                item_path = os.path.join(self.results_dir, item)

                # Check if it's a directory
                if os.path.isdir(item_path):
                    # Check if it contains an RNAseqAnalysis subdirectory
                    rnaseq_dir = os.path.join(item_path, "RNAseqAnalysis")
                    if os.path.isdir(rnaseq_dir):
                        analysis_folders.append(item_path)
                        logger.info(f"Found analysis folder: {item_path}")

        logger.info(f"Found {len(analysis_folders)} analysis folders")
        return analysis_folders

    def load_analysis_info(self) -> Dict[str, Dict]:
        """
        Load analysis information from JSON files in analysis folders.

        Returns:
        --------
        Dictionary mapping analysis IDs to their metadata
        """
        analysis_folders = self.find_analysis_folders()
        analysis_info = {}

        for folder in analysis_folders:
            analysis_id = os.path.basename(folder)
            # Check both old and new locations for analysis_info.json
            info_file_locations = [
                os.path.join(folder, "metadata", "analysis_info.json"),  # New location
                os.path.join(folder, "analysis_info.json")  # Old location for backward compatibility
            ]
            
            info_file = None
            for location in info_file_locations:
                if os.path.isfile(location):
                    info_file = location
                    break

            if info_file:
                try:
                    with open(info_file, 'r') as f:
                        info = json.load(f)

                    analysis_info[analysis_id] = info
                    logger.info(f"Loaded analysis info for {analysis_id}")
                except Exception as e:
                    logger.warning(f"Error loading analysis info from {info_file}: {str(e)}")

        self.analysis_info = analysis_info
        logger.info(f"Loaded analysis info for {len(analysis_info)} analyses")
        return analysis_info

    def find_data_files(self) -> Tuple[Dict, Dict]:
        """
        Find all DEG and CPM files across all analyses.

        Returns:
        --------
        Tuple of (deg_files, cpm_files) dictionaries
        """
        analysis_folders = self.find_analysis_folders()
        deg_files = {}
        cpm_files = {}

        for folder in analysis_folders:
            analysis_id = os.path.basename(folder)

            # Find DEG files (in contrast subdirectories)
            deg_files[analysis_id] = {}
            rnaseq_dir = os.path.join(folder, "RNAseqAnalysis")

            # Look for contrast subdirectories
            for item in os.listdir(rnaseq_dir):
                contrast_dir = os.path.join(rnaseq_dir, item)
                if os.path.isdir(contrast_dir) and item != "logs":
                    # Check for DEG.csv file
                    deg_file = os.path.join(contrast_dir, "DEG.csv")
                    if os.path.isfile(deg_file):
                        deg_files[analysis_id][item] = deg_file

            # Find CPM file
            cpm_file = os.path.join(rnaseq_dir, "CPM.csv")
            if os.path.isfile(cpm_file):
                cpm_files[analysis_id] = cpm_file

        logger.info(f"Found DEG files for {len(deg_files)} analyses")
        logger.info(f"Found CPM files for {len(cpm_files)} analyses")

        return deg_files, cpm_files

    def load_data(self):
        """
        Load all DEG and CPM files into memory.
        """
        self.dataset_info = {}  # Store dataset_info.txt content
        # Always load analysis info first
        self.load_analysis_info()

        # Find data files
        deg_files, cpm_files = self.find_data_files()

        # Load contrast information if available
        self.contrast_info = {}

        # Try to load contrasts from each analysis directory as well as global locations
        analysis_folders = self.find_analysis_folders()

        # First try to load from analysis-specific directories
        for analysis_folder in analysis_folders:
            analysis_id = os.path.basename(analysis_folder)
            for contrasts_file in [
                os.path.join(analysis_folder, "metadata", "contrasts.csv"),  # New primary location
                os.path.join(analysis_folder, "contrasts.csv")  # Old location for backward compatibility
            ]:
                if os.path.isfile(contrasts_file):
                    try:
                        contrasts_df = pd.read_csv(contrasts_file)
                        for _, row in contrasts_df.iterrows():
                            if 'name' in row and 'description' in row:
                                contrast_key = row['name']
                                # Store the extracted sections
                                self.contrast_info[contrast_key] = {
                                    'name': row['name'],  # Store the original name
                                    'description': row['description'],
                                    'expression': row.get('expression', '')
                                }
                        logger.info(f"Loaded contrasts from {contrasts_file}")
                    except Exception as e:
                        logger.warning(f"Error loading contrasts file {contrasts_file}: {str(e)}")

        # Then try global locations as fallback
        for contrasts_file in [
            os.path.join(self.results_dir, "metadata", "contrasts.csv"),
            os.path.join(self.results_dir, "contrasts.csv"),
            os.path.join(os.path.dirname(self.results_dir), "contrasts.csv")
        ]:
            if os.path.isfile(contrasts_file):
                try:
                    contrasts_df = pd.read_csv(contrasts_file)
                    for _, row in contrasts_df.iterrows():
                        if 'name' in row and 'description' in row:
                            contrast_key = row['name']
                            if contrast_key not in self.contrast_info:
                                self.contrast_info[contrast_key] = {
                                    'description': row['description'],
                                    'expression': row.get('expression', '')
                                }
                    logger.info(f"Loaded contrasts from global file {contrasts_file}")
                except Exception as e:
                    logger.warning(f"Error loading contrasts file {contrasts_file}: {str(e)}")

        logger.info(f"Loaded information for {len(self.contrast_info)} total contrasts")

        # Load DEG data
        for analysis_id, contrasts in deg_files.items():
            self.deg_data[analysis_id] = {}

            for contrast_id, deg_file in contrasts.items():
                try:
                    df = pd.read_csv(deg_file)

                    # Ensure gene column is correctly identified
                    gene_col = None
                    for col in ['Gene', 'gene', 'SYMBOL', 'symbol', 'gene_name', 'gene_id']:
                        if col in df.columns:
                            gene_col = col
                            break

                    if gene_col is None:
                        # If no specific gene column found, assume first column
                        gene_col = df.columns[0]
                        logger.warning(f"No standard gene column found in {deg_file}, using {gene_col}")

                    # Standardize column name to 'Gene'
                    df = df.rename(columns={gene_col: 'Gene'})

                    # Store in dictionary
                    self.deg_data[analysis_id][contrast_id] = df
                    logger.debug(f"Loaded DEG data for {analysis_id}/{contrast_id} with {len(df)} genes")

                except Exception as e:
                    logger.error(f"Error loading DEG file {deg_file}: {str(e)}")

        # Load CPM data
        for analysis_id, cpm_file in cpm_files.items():
            try:
                df = pd.read_csv(cpm_file)

                # Ensure gene column is correctly identified
                gene_col = None
                for col in ['Gene', 'gene', 'SYMBOL', 'symbol', 'gene_name', 'gene_id']:
                    if col in df.columns:
                        gene_col = col
                        break

                if gene_col is None:
                    # If no specific gene column found, assume first column
                    gene_col = df.columns[0]
                    logger.warning(f"No standard gene column found in {cpm_file}, using {gene_col}")

                # Standardize column name to 'Gene'
                df = df.rename(columns={gene_col: 'Gene'})

                # Store in dictionary
                self.cpm_data[analysis_id] = df
                logger.debug(f"Loaded CPM data for {analysis_id} with {len(df)} genes and {len(df.columns)-1} samples")

            except Exception as e:
                logger.error(f"Error loading CPM file {cpm_file}: {str(e)}")

        # Load dataset_info.txt files if available
        self._load_dataset_info()

        logger.info(f"Data loading complete: {len(self.deg_data)} DEG datasets, {len(self.cpm_data)} CPM datasets")

    def identify_important_genes(self,
                               top_frequent: int = 20,
                               top_unique: int = 10,
                               max_contrasts_for_unique: int = 2,
                               min_unique_per_contrast: int = 1,
                               p_value_threshold: float = None,
                               lfc_threshold: float = None) -> Set[str]:
        """
        Identify important genes across all analyses using two criteria:
        1. Top frequently occurring genes across all contrasts
        2. Genes that appear in few contrasts but have large fold changes

        Parameters:
        -----------
        top_frequent : int
            Number of top frequently occurring genes to select
        top_unique : int
            Number of top unique genes to select per contrast based on fold change
        max_contrasts_for_unique : int
            Maximum number of contrasts a gene can appear in to be considered "unique"

        Returns:
        --------
        Set of important gene identifiers
        """
        important_genes = set()

        # Dictionary to track gene occurrences across contrasts
        genes_by_occurrence = {}

        # Dictionary to track genes and their log fold changes in each contrast
        genes_by_contrast = {}

        # Use provided thresholds or fall back to instance defaults
        p_thresh = p_value_threshold if p_value_threshold is not None else self.pvalue_threshold
        lfc_thresh = lfc_threshold if lfc_threshold is not None else self.lfc_threshold

        # Collect significant genes and their properties
        logger.info(f"Identifying important genes across {len(self.deg_data)} analyses")
        for analysis_id, contrasts in self.deg_data.items():
            logger.info(f"Processing {len(contrasts)} contrasts in {analysis_id}")
            for contrast_id, df in contrasts.items():
                contrast_key = f"{analysis_id}_{contrast_id}"
                logger.info(f"Evaluating contrast: {contrast_key} with {len(df)} genes")

                # Skip if Gene column doesn't exist
                if 'Gene' not in df.columns:
                    logger.warning(f"No 'Gene' column in {contrast_key}")
                    continue

                # Use exact column names from the DEG.csv file
                # Find adjusted p-value column first, fall back to unadjusted if necessary
                p_value_col = None
                if 'adj.P.Val' in df.columns:
                    p_value_col = 'adj.P.Val'  # Preferred: use adjusted p-value for multiple testing correction
                    logger.info(f"Using adjusted p-value column 'adj.P.Val' for {contrast_key}")
                elif 'P.Value' in df.columns:
                    p_value_col = 'P.Value'  # Fall back to unadjusted p-value if adjusted not available
                    logger.info(f"Using unadjusted p-value column 'P.Value' for {contrast_key}")

                # Use exact column name from the DEG.csv file
                lfc_col = 'logFC' if 'logFC' in df.columns else None

                if p_value_col is None or lfc_col is None:
                    logger.warning(f"Missing required columns in {contrast_key}. Available columns: {df.columns.tolist()}")
                    continue

                logger.info(f"Using columns for {contrast_key}: p-value={p_value_col}, logFC={lfc_col}")

                # Filter significant genes
                sig_df = df[(df[p_value_col] < p_thresh) &
                           (abs(df[lfc_col]) > lfc_thresh)]

                # Store significant genes for this contrast with their fold changes
                genes_by_contrast[contrast_key] = {}
                for _, row in sig_df.iterrows():
                    gene = row['Gene']
                    # Store the gene with its fold change
                    genes_by_contrast[contrast_key][gene] = abs(row[lfc_col])

                # Count occurrences of each gene in significant results
                for gene in sig_df['Gene']:
                    if gene not in genes_by_occurrence:
                        genes_by_occurrence[gene] = 1
                    else:
                        genes_by_occurrence[gene] += 1

        # 1. Get the top frequently occurring genes
        sorted_by_occurrence = sorted(genes_by_occurrence.items(), key=lambda x: x[1], reverse=True)
        top_frequent_genes = [gene for gene, _ in sorted_by_occurrence[:top_frequent]]
        important_genes.update(top_frequent_genes)

        logger.info(f"Selected {len(top_frequent_genes)} top frequently occurring genes")

        # 2. Find unique/rare genes with large fold changes in each contrast
        for contrast_key, genes in genes_by_contrast.items():
            # Find genes that appear in few contrasts (unique or rare)
            unique_genes = {gene: lfc for gene, lfc in genes.items()
                           if genes_by_occurrence.get(gene, 0) <= max_contrasts_for_unique}

            # Sort by fold change and take top ones
            sorted_unique_genes = sorted(unique_genes.items(), key=lambda x: x[1], reverse=True)
            top_unique_genes = [gene for gene, _ in sorted_unique_genes[:top_unique]]

            if len(top_unique_genes) >= min_unique_per_contrast:
                important_genes.update(top_unique_genes)
                logger.info(f"Selected {len(top_unique_genes)} unique genes with large fold changes from {contrast_key}")
            else:
                logger.warning(f"Not enough unique genes found for {contrast_key}: {len(top_unique_genes)} < {min_unique_per_contrast}")

        logger.info(f"Identified {len(important_genes)} important genes in total")
        return important_genes

    def _ensure_output_dir(self):
        """Create output_dir only on first use."""
        if self.output_dir is None:  # not created yet
            if self._output_dir_requested:
                self.output_dir = os.path.abspath(self._output_dir_requested)
            else:
                timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                self.output_dir = os.path.join(self.results_dir, f"integrated_results_{timestamp}")

            os.makedirs(self.output_dir, exist_ok=True)
            logger.info(f"Output will be saved to: {self.output_dir}")

        return self.output_dir

    def create_lfc_heatmap(self,
                         genes: List[str] = None,
                         contrasts: List[Tuple[str, str]] = None,
                         output_file: str = None,
                         p_value_threshold: float = None,
                         lfc_threshold: float = None,
                         hide_empty_rows_cols: bool = False,
                         font_size: int = 12,
                         show_grid_lines: bool = True,
                         grid_opacity: float = 0.3) -> go.Figure:
        """
        Create an interactive heatmap of log fold changes for selected genes across contrasts.

        Parameters:
        -----------
        genes : List[str], optional
            List of genes to include in the heatmap. If None, important genes are identified.
        contrasts : List[Tuple[str, str]], optional
            List of (analysis_id, contrast_id) tuples to include. If None, all contrasts are used.
        output_file : str, optional
            Path to save the interactive HTML output. If None, no file is saved (useful for Streamlit).
        p_value_threshold : float, optional
            P-value threshold for significance filtering. If None, uses instance default.
        lfc_threshold : float, optional
            Log fold change threshold for significance filtering. If None, uses instance default.
        hide_empty_rows_cols : bool, optional
            If True, hide genes/contrasts where no significant values exist.
        font_size : int, optional
            Font size for labels and annotations
        show_grid_lines : bool, optional
            Whether to show grid lines
        grid_opacity : float, optional
            Opacity of grid lines

        Returns:
        --------
        plotly.graph_objects.Figure
        """
        # Use provided thresholds or fall back to instance defaults
        p_thresh = p_value_threshold if p_value_threshold is not None else self.pvalue_threshold
        lfc_thresh = lfc_threshold if lfc_threshold is not None else self.lfc_threshold

        # If genes not provided, identify important genes
        if genes is None:
            try:
                genes = list(self.identify_important_genes(
                    top_frequent=20,
                    top_unique=10,
                    max_contrasts_for_unique=2,
                    min_unique_per_contrast=1
                ))
                if not genes:
                    logger.warning("No important genes identified")
                    return None
                logger.info(f"Automatically identified {len(genes)} important genes")
            except Exception as e:
                logger.error(f"Error identifying important genes: {str(e)}")
                return None

        # Prepare data for heatmap
        heatmap_data = []
        logger.info(f"Creating LFC heatmap with {len(genes)} genes and {len(contrasts) if contrasts else 'all'} contrasts")

        # If contrasts not specified, use all
        if contrasts is None:
            contrasts = []
            for analysis_id, contrast_dict in self.deg_data.items():
                for contrast_id in contrast_dict.keys():
                    contrasts.append((analysis_id, contrast_id))
            logger.info(f"No contrasts specified, using all {len(contrasts)} contrasts")

        # Create simplified contrast labels for display
        def simplify_contrast_label(analysis_id: str, contrast_id: str) -> str:
            """Create simplified label for display"""
            # Try to get original name from contrast_info if available
            if hasattr(self, "contrast_info") and contrast_id in self.contrast_info:
                original_name = self.contrast_info[contrast_id].get('name', contrast_id)
            else:
                original_name = contrast_id

            # Remove dataset prefix and keep it short
            simplified = original_name.split(":", 1)[-1].split(" –")[0].split(" - ")[0][:25]
            return simplified

        # Create both full and simplified contrast labels
        contrast_labels = [f"{a_id}_{c_id}" for a_id, c_id in contrasts]
        simplified_labels = [simplify_contrast_label(a_id, c_id) for a_id, c_id in contrasts]
        logger.info(f"Using {len(contrast_labels)} contrast labels: {', '.join(simplified_labels[:5])}{' ...' if len(simplified_labels) > 5 else ''}")

        # Extract log fold change values for each gene across contrasts
        for gene in genes:
            row = {'Gene': gene}

            for (analysis_id, contrast_id), contrast_label in zip(contrasts, contrast_labels):
                if analysis_id in self.deg_data and contrast_id in self.deg_data[analysis_id]:
                    df = self.deg_data[analysis_id][contrast_id]

                    if 'Gene' not in df.columns:
                        logger.warning(f"No 'Gene' column found in {analysis_id}/{contrast_id}")
                        row[contrast_label] = 0  # Skip this contrast
                        continue

                    # Find log fold change column
                    lfc_col = None
                    for col in ['logFC', 'log2FoldChange', 'log2FC', 'LogFC']:
                        if col in df.columns:
                            lfc_col = col
                            break

                    if lfc_col is None:
                        logger.warning(f"No log fold change column found in {analysis_id}/{contrast_id}")
                        row[contrast_label] = 0  # Skip this contrast
                        continue

                    # Use exact column names from the DEG.csv file
                    # Find adjusted p-value column first, fall back to unadjusted if necessary
                    p_value_col = None
                    if 'adj.P.Val' in df.columns:
                        p_value_col = 'adj.P.Val'  # Preferred: use adjusted p-value for multiple testing correction
                        logger.debug(f"Using adjusted p-value column 'adj.P.Val' for {analysis_id}/{contrast_id}")
                    elif 'P.Value' in df.columns:
                        p_value_col = 'P.Value'  # Fall back to unadjusted p-value if adjusted not available
                        logger.debug(f"Using unadjusted p-value column 'P.Value' for {analysis_id}/{contrast_id}")

                    if p_value_col is None:
                        logger.warning(f"No p-value column found in {analysis_id}/{contrast_id}")
                        # Continue anyway, we'll use the log fold change without significance filtering

                    # Find the gene and get its logFC, but set to 0 if not significant
                    gene_row = df[df['Gene'] == gene]
                    if not gene_row.empty:
                        if p_value_col is not None:
                            p_value = gene_row.iloc[0][p_value_col]
                            abs_lfc = abs(gene_row.iloc[0][lfc_col])
                            if p_value >= p_thresh or abs_lfc <= lfc_thresh:
                                row[contrast_label] = 0  # Not significant, set to white
                            else:
                                row[contrast_label] = gene_row.iloc[0][lfc_col]
                        else:
                            row[contrast_label] = gene_row.iloc[0][lfc_col]
                    else:
                        # Gene not found in this contrast
                        row[contrast_label] = 0  # Use 0 for non-significant/missing
                else:
                    logger.debug(f"Contrast {analysis_id}/{contrast_id} not found in DEG data")
                    row[contrast_label] = 0  # Skip this contrast

            # Only add the row if there's at least one logFC value
            if len(row) > 1:  # More than just the Gene column
                heatmap_data.append(row)
            else:
                logger.warning(f"No logFC values found for gene {gene}, skipping")

        # Convert to DataFrame
        heatmap_df = pd.DataFrame(heatmap_data)

        # Log the shape of the data
        if len(heatmap_data) > 0:
            logger.info(f"Heatmap data shape: {len(heatmap_df)} rows (genes) × {len(heatmap_df.columns)} columns (including Gene column)")
            if len(heatmap_df.columns) > 1:
                logger.info(f"Heatmap columns: {', '.join(heatmap_df.columns[:10])}{' ...' if len(heatmap_df.columns) > 10 else ''}")

        if len(heatmap_df) == 0:
            logger.warning("No data for heatmap visualization - no matching genes found. Try selecting more genes or adjusting thresholds.")
            return None
        elif len(heatmap_df.columns) <= 1:
            logger.warning("Not enough data for heatmap visualization - no contrast data found")
            return None

        # Count non-zero values to ensure we have meaningful data
        if 'Gene' in heatmap_df.columns:
            heatmap_df_values = heatmap_df.drop(columns=['Gene'])
        else:
            heatmap_df_values = heatmap_df

        if heatmap_df_values.shape[1] == 0:
            logger.warning("No contrast columns found in heatmap data. Please select at least one valid contrast.")
            return None

        # Count non-zero entries
        non_zero_count = (heatmap_df_values != 0).sum().sum()
        total_cells = heatmap_df_values.size
        non_zero_percent = (non_zero_count / total_cells) * 100 if total_cells > 0 else 0

        logger.info(f"Heatmap contains {non_zero_count}/{total_cells} non-zero values ({non_zero_percent:.1f}%)")

        # Ensure we have at least some non-zero values
        if non_zero_count == 0:
            logger.warning("No non-zero values found in heatmap data - all genes may be non-significant")
            # Continue anyway, but log a warning

        # Prepare for heatmap - set Gene as index and fill NAs
        if 'Gene' in heatmap_df.columns:
            heatmap_df = heatmap_df.set_index('Gene')
        heatmap_df = heatmap_df.fillna(0)

        # Remove empty rows/columns if requested
        if hide_empty_rows_cols:
            # Remove genes (rows) where all values are 0
            non_zero_rows = (heatmap_df != 0).any(axis=1)
            heatmap_df = heatmap_df[non_zero_rows]

            # Remove contrasts (columns) where all values are 0
            non_zero_cols = (heatmap_df != 0).any(axis=0)
            heatmap_df = heatmap_df.loc[:, non_zero_cols]

            # Update contrasts list to match remaining columns
            if len(heatmap_df.columns) < len(contrasts):
                remaining_contrast_indices = [i for i, label in enumerate(contrast_labels) if label in heatmap_df.columns]
                contrasts = [contrasts[i] for i in remaining_contrast_indices]
                logger.info(f"After filtering empty rows/columns: {len(heatmap_df)} genes × {len(heatmap_df.columns)} contrasts")

        # Check if we still have data after filtering
        if len(heatmap_df) == 0:
            logger.warning("No genes remain after filtering empty rows")
            return None
        elif len(heatmap_df.columns) == 0:
            logger.warning("No contrasts remain after filtering empty columns")
            return None

        # Perform hierarchical clustering on both axes if there are enough data points
        # Perform hierarchical clustering on both axes
        # 1. Cluster genes (rows)
        if len(heatmap_df) > 1:  # Only cluster if we have multiple genes
            try:
                # Only cluster if we have non-zero values
                if non_zero_count > 0:
                    gene_linkage = sch.linkage(pdist(heatmap_df.values), method='complete')
                    gene_order = sch.leaves_list(gene_linkage)
                    clustered_genes = heatmap_df.index[gene_order].tolist()
                else:
                    # Skip clustering if all values are zero
                    logger.warning("Skipping gene clustering due to lack of non-zero values")
                    clustered_genes = heatmap_df.index.tolist()
            except Exception as e:
                logger.warning(f"Error clustering genes: {str(e)}. Using original order.")
                clustered_genes = heatmap_df.index.tolist()
        else:
            clustered_genes = heatmap_df.index.tolist()

        # 2. Cluster contrasts (columns)
        if len(heatmap_df.columns) > 1:  # Only cluster if we have multiple contrasts
            try:
                # Only cluster if we have non-zero values
                if non_zero_count > 0:
                    contrast_linkage = sch.linkage(pdist(heatmap_df.values.T), method='complete')
                    contrast_order = sch.leaves_list(contrast_linkage)
                    clustered_contrasts = [heatmap_df.columns[i] for i in contrast_order]
                else:
                    # Skip clustering if all values are zero
                    logger.warning("Skipping contrast clustering due to lack of non-zero values")
                    clustered_contrasts = heatmap_df.columns.tolist()
            except Exception as e:
                logger.warning(f"Error clustering contrasts: {str(e)}. Using original order.")
                clustered_contrasts = heatmap_df.columns.tolist()
        else:
            clustered_contrasts = heatmap_df.columns.tolist()

        # Reorder the DataFrame according to clustering
        heatmap_df = heatmap_df.loc[clustered_genes, clustered_contrasts]
        logger.info(f"Final heatmap dimensions: {len(clustered_genes)} genes × {len(clustered_contrasts)} contrasts")

        # Create heatmap using plotly
        max_abs_value = max(abs(heatmap_df.values.min()), abs(heatmap_df.values.max()))
        if max_abs_value == 0:
            max_abs_value = 1.0  # Avoid division by zero
            logger.warning("All values in heatmap are zero - using default color scale range of -1 to 1")

        # Create the heatmap
        fig = px.imshow(
            heatmap_df,
            color_continuous_scale='RdBu_r',  # Red-Blue color scale, reversed to match biology convention
            labels=dict(x="Contrast", y="Gene", color="Log2FC"),
            title="Differential Expression Heatmap (Log2 Fold Change)",
            aspect="auto",  # Maintain aspect ratio based on data
            height=max(600, min(2000, len(heatmap_df) * 22)),  # Taller height based on number of genes (capped at 2000px)
            width=max(600, min(1500, len(heatmap_df.columns) * 80)),  # Dynamic width based on number of contrasts
            zmin=-max_abs_value,  # Symmetrical color scale
            zmax=max_abs_value,
            color_continuous_midpoint=0  # Ensure white is at 0
        )

        # Improve layout
        fig.update_layout(
            margin=dict(l=250, r=20, t=60, b=120),
            coloraxis_colorbar=dict(title="Log2FC"),
            font_family="Inter, sans-serif",
            font_size=font_size
        )

        # Configure grid lines
        if show_grid_lines:
            fig.update_yaxes(showgrid=True, gridwidth=0.5, gridcolor=f"rgba(200,200,200,{grid_opacity})")
            fig.update_xaxes(showgrid=True, gridwidth=0.5, gridcolor=f"rgba(200,200,200,{grid_opacity})")
        else:
            fig.update_yaxes(showgrid=False)
            fig.update_xaxes(showgrid=False)

        # We'll create custom hover data with descriptions in the update_traces call below

        # Create a 3D array with gene logFC, contrast description, accession info, and original contrast name
        # First dimension: rows (genes)
        # Second dimension: columns (contrasts)
        # Third dimension: data types (0: contrast description, 1: accession, 2: original contrast name)
        hover_info = np.empty((*heatmap_df.shape, 3), dtype=object)

        # Keep track of original contrast names for both hover and axis labels
        original_contrast_names = []

        # Fill with contrast descriptions and accession info
        for i, gene in enumerate(heatmap_df.index):
            for j, col_idx in enumerate(heatmap_df.columns):
                # Get contrast information from the original contrasts list
                analysis_id, contrast_id = contrasts[j]

                # Get contrast description using the helper method
                description = self._get_contrast_description(analysis_id, contrast_id)

                # Get accession ID
                accession = self.analysis_info.get(analysis_id, {}).get('accession', analysis_id)

                # Get simplified name for hover consistency
                simplified_name = simplified_labels[j]
                if j == 0 and i == 0:  # Only collect once per contrast
                    original_contrast_names.append(simplified_name)

                # Store description, accession, and simplified contrast name
                hover_info[i, j, 0] = description
                hover_info[i, j, 1] = accession
                hover_info[i, j, 2] = simplified_name

        # Update hover template with conditional formatting, accession, and description
        hover_template = (
            "<b>Gene:</b> %{y}<br>" +
            "<b>Contrast:</b> %{customdata[2]}<br>" +  # Use original contrast name
            "<b>Log2FC:</b> %{z:.2f}<br>" +
            "<b>Accession:</b> %{customdata[1]}<br>" +
            "<b>Description:</b> %{customdata[0]}<br>" +
            "<extra></extra>"  # Hide secondary info
        )
        fig.update_traces(
            hovertemplate=hover_template,
            customdata=hover_info,
            # Use custom color for hover
            hoverlabel=dict(bgcolor="white", font_size=12, font_family="Arial")
        )

        # Update x-axis with simplified contrast names
        if simplified_labels:
            # Log simplified contrast names for debugging
            logger.debug(f"Original contrast IDs: {[cid for _, cid in contrasts]}")
            logger.debug(f"Using simplified contrast names: {simplified_labels}")

            # Ensure we have labels for all contrasts
            assert len(simplified_labels) == len(contrasts), "Mismatch between number of contrasts and simplified names"

            # Build new tick positions after clustering (need to map to clustered order)
            clustered_simplified_labels = [simplified_labels[contrast_labels.index(col)] for col in clustered_contrasts]

            # Build new tick positions after clustering
            tick_vals = list(range(len(clustered_simplified_labels)))
            fig.update_xaxes(
                tickmode="array",
                tickvals=tick_vals,
                ticktext=clustered_simplified_labels,
                tickangle=45,
                title="Biological contrast"
            )

        # Save to file if specified
        if output_file is not None:
            self._ensure_output_dir()
            fig.write_html(output_file)
            logger.info(f"Saved LFC heatmap to {output_file}")

        return fig

    def create_expression_plots(self,
                              genes: List[str],
                              plot_type: str = 'violin',
                              analyses: List[str] = None,
                              output_file: str = None,
                              hide_x_labels: bool = False,
                              page_number: int = 1,
                              facet_font_size: int = 10,
                              lock_y_axis: bool = False,
                              show_raw_points: bool = True,
                              legend_position: str = "Bottom",
                              show_grid_lines: bool = True,
                              grid_opacity: float = 0.3) -> go.Figure:
        """
        Create interactive expression plots for selected genes across samples.

        Parameters:
        -----------
        genes : List[str]
            List of genes to plot expression for
        plot_type : str
            Type of plot: 'violin', 'box', or 'both'
        analyses : List[str], optional
            List of analysis IDs to include. If None, all analyses are used.
        output_file : str, optional
            Path to save the interactive HTML output. If None, no file is saved (useful for Streamlit).
        facet_font_size : int, optional
            Font size for facet titles
        lock_y_axis : bool, optional
            Whether to use same y-axis range across all facets
        show_raw_points : bool, optional
            Whether to show individual data points
        legend_position : str, optional
            Legend position: "Bottom", "Right", or "Top"
        show_grid_lines : bool, optional
            Whether to show grid lines
        grid_opacity : float, optional
            Opacity of grid lines

        Returns:
        --------
        plotly.graph_objects.Figure
        """
        if not self.cpm_data:
            logger.error("No CPM data loaded. Cannot create expression plots.")
            return None

        # If analyses not specified, use all
        if analyses is None:
            analyses = list(self.cpm_data.keys())

        # Create a tidy (long) dataframe for all expression data
        all_long_data = []

        # Maximum number of genes to process in one page
        genes_per_page = 30
        # Save original gene list
        gene_list = genes.copy() if genes is not None else []

        # If no genes are provided, return None
        if not gene_list:
            logger.warning("No genes provided for expression plots")
            return None

        gene_pages = [gene_list[i:i + genes_per_page] for i in range(0, len(gene_list), genes_per_page)]

        if len(gene_pages) > 1:
            logger.info(f"Splitting {len(gene_list)} genes into {len(gene_pages)} pages for visualization")
            # Process first page only
            genes = gene_pages[0]

        # Process each analysis and build tidy data
        for analysis_id in analyses:
            if analysis_id not in self.cpm_data:
                logger.warning(f"No CPM data for analysis {analysis_id}")
                continue

            # Get CPM data for this analysis
            cpm_df = self.cpm_data[analysis_id].copy()

            # Skip if Gene column doesn't exist
            if 'Gene' not in cpm_df.columns:
                logger.warning(f"No 'Gene' column in CPM data for {analysis_id}")
                continue

            # Filter to only the genes we're working with
            cpm_filtered = cpm_df[cpm_df['Gene'].isin(genes)]

            if cpm_filtered.empty:
                logger.warning(f"No matching genes found in CPM data for {analysis_id}")
                continue

            # Get the sample groups from the edger_analysis_samples.csv file
            sample_to_group = {}

            # Helper function to find first existing file
            def first_existing(paths):
                for path in paths:
                    if os.path.isfile(path):
                        return path
                return None

            # Look inside the current analysis directory first
            analysis_dir = os.path.join(self.results_dir, analysis_id)
            sample_file_locations = [
                os.path.join(analysis_dir, "metadata", "edger_analysis_samples.csv"),  # New primary location
                os.path.join(analysis_dir, "edger_analysis_samples.csv"),  # Old location for backward compatibility
                os.path.join(self.results_dir, "edger_analysis_samples.csv"),
                os.path.join(os.path.dirname(self.results_dir), "edger_analysis_samples.csv")
            ]

            edger_samples_file = first_existing(sample_file_locations)
            if edger_samples_file:
                logger.info(f"Found sample mapping file at {edger_samples_file}")

            if edger_samples_file:
                try:
                    # First try loading the file to check its structure
                    sample_mapping_df = pd.read_csv(edger_samples_file, nrows=0)

                    # If the first column has no name or starts with 'Unnamed', it's probably an index
                    first_col = sample_mapping_df.columns[0]
                    has_index = first_col == '' or first_col.startswith('Unnamed')

                    # Now load the full file with appropriate index_col setting
                    sample_mapping_df = pd.read_csv(edger_samples_file, index_col=0 if has_index else None)
                    logger.info(f"Loaded sample mapping with columns: {', '.join(sample_mapping_df.columns)}")

                    # Find the group column
                    group_col = None
                    if self.analysis_info and analysis_id in self.analysis_info and 'analysis_column' in self.analysis_info[analysis_id]:
                        group_col = self.analysis_info[analysis_id]['analysis_column']
                        logger.info(f"Using analysis column '{group_col}' from analysis_info")
                    elif 'merged_analysis_group' in sample_mapping_df.columns:
                        group_col = 'merged_analysis_group'
                        logger.info(f"Using 'merged_analysis_group' column")
                    else:
                        for col in sample_mapping_df.columns:
                            if 'group' in col.lower() or 'condition' in col.lower():
                                group_col = col
                                logger.info(f"Using '{col}' column as group")
                                break

                    if group_col and group_col in sample_mapping_df.columns:
                        # Look for a 'Sample1', 'Sample2', etc. naming pattern in CPM columns
                        cpm_columns = self.cpm_data[analysis_id].columns.tolist()
                        sample_cols = [col for col in cpm_columns if col.startswith('Sample') and col != 'Gene']

                        if sample_cols:
                            # Create a mapping from sample name to group
                            for i in range(len(sample_mapping_df)):
                                if i < len(sample_cols):  # Ensure we don't go out of bounds
                                    sample_key = sample_cols[i]
                                    sample_to_group[sample_key] = sample_mapping_df.iloc[i][group_col]
                            logger.info(f"Created mapping for {len(sample_to_group)} samples using Sample# pattern")
                        else:
                            # Fall back to index-based mapping
                            for i in range(len(sample_mapping_df)):
                                # Use sample index (i) as the key (starting with Sample1)
                                sample_key = f"Sample{i+1}"
                                sample_to_group[sample_key] = sample_mapping_df.iloc[i][group_col]
                            logger.info(f"Created mapping for {len(sample_to_group)} samples using index")
                    else:
                        logger.warning(f"No group column found in sample mapping file")
                except Exception as e:
                    logger.warning(f"Error loading sample mapping file: {str(e)}")
            else:
                logger.warning("Could not find edger_analysis_samples.csv file. Using default sample names.")

            # Convert to long format
            melted_df = cpm_filtered.melt(
                id_vars='Gene',
                var_name='Sample',
                value_name='Expression'
            )

            # Add analysis ID
            melted_df['Analysis'] = analysis_id

            # Create a more descriptive sample label for hover info
            melted_df['SampleLabel'] = melted_df['Analysis'] + ":" + melted_df['Sample']

            # Map samples to groups if available
            if sample_to_group:
                melted_df['Group'] = melted_df['Sample'].map(sample_to_group)
                # Fill any missing mappings
                melted_df['Group'] = melted_df['Group'].fillna(melted_df['Sample'])
            else:
                # Use analysis ID as group if no mapping available
                melted_df['Group'] = analysis_id

            # Append to all data
            all_long_data.append(melted_df)

        # Combine all data
        if not all_long_data:
            logger.warning("No expression data found for specified genes and analyses")
            return None

        # Concatenate all data
        plot_df = pd.concat(all_long_data, ignore_index=True)

        if plot_df.empty:
            logger.warning("No expression data after filtering")
            return None

        # Number of genes actually found in data
        genes_found = plot_df['Gene'].nunique()

        # Determine facet column wrapping (4 columns max)
        facet_col_wrap = min(4, genes_found)

        # Calculate rows for height
        n_rows = math.ceil(genes_found / facet_col_wrap)

        # Apply log transformation for better visualization if needed
        if plot_df['Expression'].min() < 0 or (plot_df['Expression'].min() >= 0 and plot_df['Expression'].max() <= 30):
            # Data appears to be already log-transformed
            plot_df['LogExpression'] = plot_df['Expression']
            y_axis_title = "Log2(CPM)"
        else:
            # Apply log2 transformation (after adding 1 to avoid log(0))
            plot_df['LogExpression'] = np.log2(plot_df['Expression'] + 1)
            y_axis_title = "Log2(CPM+1)"

        # Use simple consistent spacing values that work well across different panel counts
        facet_row_spacing = 0.03
        facet_col_spacing = 0.03

        # Add sample label for hover info
        if 'SampleLabel' not in plot_df.columns:
            plot_df['SampleLabel'] = plot_df['Analysis'] + ":" + plot_df['Sample']

        # Create violin plot using Plotly Express with adjusted spacing
        fig = px.violin(
            plot_df,
            x='Group',
            y='LogExpression',
            color='Analysis',
            facet_col='Gene',
            facet_col_wrap=facet_col_wrap,
            hover_data=['SampleLabel', 'Group'],
            box=True,
            points="all",
            title=f"Gene Expression Across Samples (Page {page_number} of {len(gene_pages)})",
            labels={'LogExpression': y_axis_title, 'Group': 'Group', 'Analysis': 'Dataset'},
            facet_row_spacing=facet_row_spacing,
            facet_col_spacing=facet_col_spacing
        )

        # Layout adjustments optimized for readability
        fig.update_layout(
            height=350 * n_rows,  # 350px per row for more vertical space
            width=320 * facet_col_wrap,  # 320px per column
            legend_title_text="Dataset",
            margin=dict(l=40, r=20, t=80, b=40),
            font_family="Inter, sans-serif"
        )

        # Configure legend position
        if legend_position == "Bottom":
            fig.update_layout(legend=dict(orientation="h", y=-0.35, x=0.5, xanchor="center", yanchor="top"))
        elif legend_position == "Right":
            fig.update_layout(legend=dict(orientation="v", x=1.02, y=0.5, xanchor="left", yanchor="middle"))
        elif legend_position == "Top":
            fig.update_layout(legend=dict(orientation="h", y=1.15, x=0.5, xanchor="center", yanchor="bottom"))
        else:  # Default to bottom
            fig.update_layout(legend=dict(orientation="h", y=-0.35, x=0.5, xanchor="center", yanchor="top"))

        # Improve axis appearance
        fig.update_xaxes(categoryorder="category ascending")
        fig.update_layout(xaxis_title="")

        # Configure grid lines
        if show_grid_lines:
            fig.update_yaxes(showgrid=True, gridwidth=0.5, gridcolor=f"rgba(200,200,200,{grid_opacity})")
        else:
            fig.update_yaxes(showgrid=False)

        # Configure points visibility and transparency
        if show_raw_points:
            fig.update_traces(marker=dict(size=4, opacity=0.5), jitter=0.3)
        else:
            fig.update_traces(marker=dict(size=0), jitter=0.3)

        # Update facet formatting (remove "Gene=" prefix and adjust font size)
        fig.for_each_annotation(lambda a: a.update(
            text=a.text.split("=")[1],
            font=dict(size=facet_font_size, color="#333", family="Inter, sans-serif")
        ))

        # Lock y-axis if requested
        if lock_y_axis and not plot_df.empty:
            global_min = plot_df['LogExpression'].min()
            global_max = plot_df['LogExpression'].max()
            # Add some padding
            padding = (global_max - global_min) * 0.05
            fig.update_yaxes(range=[global_min - padding, global_max + padding])

        # Hide x-axis labels if requested
        if hide_x_labels:
            fig.update_xaxes(showticklabels=False)

        # Add note about number of genes
        if len(genes) < len(gene_list):
            fig.add_annotation(
                text=f"Showing {genes_found} genes (page {page_number} of {len(gene_pages)})",
                x=0.5, y=1.05,
                xref="paper", yref="paper",
                showarrow=False,
                font=dict(size=12),
                align="center"
            )
        elif genes_found < len(genes):
            fig.add_annotation(
                text=f"Only {genes_found} of {len(genes)} genes found in data",
                x=0.5, y=1.05,
                xref="paper", yref="paper",
                showarrow=False,
                font=dict(size=12, color="red"),
                align="center"
            )

        # Save to file if specified
        if output_file is not None:
            self._ensure_output_dir()
            fig.write_html(output_file)
            logger.info(f"Saved expression plots to {output_file}")

        return fig

    def _get_contrast_description(self, analysis_id: str, contrast_id: str) -> str:
        """Get the description for a contrast from the loaded contrast info."""
        # Try exact match first
        if hasattr(self, "contrast_info") and contrast_id in self.contrast_info:
            return self.contrast_info[contrast_id].get('description', f"Contrast: {contrast_id}")

        # Try simplified pattern matching (without dataset/accession prefix)
        if "_" in contrast_id:
            # Try parts after potential prefixes
            parts = contrast_id.split("_")
            for i in range(1, len(parts)):
                possible_id = "_".join(parts[i:])
                if hasattr(self, "contrast_info") and possible_id in self.contrast_info:
                    return self.contrast_info[possible_id].get('description', f"Contrast: {possible_id}")

        # Try to look in an analysis-specific contrasts file (check multiple locations)
        for contrasts_file in [
            os.path.join(self.results_dir, analysis_id, "metadata", "contrasts.csv"),  # New primary location
            os.path.join(self.results_dir, analysis_id, "contrasts.csv"),  # Old location for backward compatibility
            os.path.join(self.results_dir, "contrasts.csv")  # Also check global location
        ]:
            if os.path.isfile(contrasts_file):
                try:
                    contrasts_df = pd.read_csv(contrasts_file)
                    # Try exact match on 'name' field
                    match = contrasts_df[contrasts_df['name'] == contrast_id]
                    if not match.empty:
                        return match.iloc[0].get('description', f"Contrast: {contrast_id}")

                    # Try parts of the contrast ID for partial matches
                    if "_" in contrast_id:
                        parts = contrast_id.split("_")
                        for i in range(1, len(parts)):
                            possible_id = "_".join(parts[i:])
                            match = contrasts_df[contrasts_df['name'] == possible_id]
                            if not match.empty:
                                return match.iloc[0].get('description', f"Contrast: {possible_id}")
                except Exception as e:
                    logger.debug(f"Error loading contrasts from {contrasts_file}: {str(e)}")

        return f"Contrast: {contrast_id}"

    def _load_dataset_info(self):
        """
        Load dataset_info.txt files for each analysis and extract title, summary, and design sections.

        The expected format of dataset_info.txt is:
        - Title at the beginning
        - Summary in the middle section
        - Overall design at the end (usually starts with "Overall design:")

        The function attempts to parse these sections using both regex patterns and fallback line-by-line parsing.
        """
        import re  # Ensure re is imported for pattern matching
        analysis_folders = self.find_analysis_folders()

        for analysis_folder in analysis_folders:
            analysis_id = os.path.basename(analysis_folder)
            # Look for dataset_info.txt in various possible locations (though it's now stored in analysis_info.json)
            for info_path in [
                os.path.join(analysis_folder, "metadata", "dataset_info.txt"),
                os.path.join(analysis_folder, "dataset_info.txt"),
                os.path.join(self.results_dir, "metadata", "dataset_info.txt")
            ]:
                if os.path.isfile(info_path):
                    try:
                        with open(info_path, 'r') as f:
                            content = f.read()

                        # Initialize dataset info with empty values
                        self.dataset_info[analysis_id] = {
                            "title": "",
                            "summary": "",
                            "design": ""
                        }

                        # Try to extract sections using common patterns
                        import re

                        # Match common section patterns
                        title_match = re.search(r'^(.*?)(?=\n\s*\n|\n*Overall design:|\n*Summary:|\Z)', content, re.DOTALL | re.IGNORECASE)
                        design_match = re.search(r'(?:Overall design:)(.*?)(?=\n\s*\n|\Z)', content, re.DOTALL | re.IGNORECASE)

                        # If we can't find the design section explicitly, look for it at the end
                        if not design_match:
                            design_match = re.search(r'(?:Overall design\s*[-:])?(.*?)(?=\Z)', content, re.DOTALL | re.IGNORECASE)

                        # Extract the summary (everything between title and design)
                        if title_match:
                            title_end = title_match.end()
                            design_start = design_match.start() if design_match else len(content)
                            summary_text = content[title_end:design_start].strip()

                            # Clean up summary text - remove any "Summary:" prefix
                            summary_text = re.sub(r'^Summary:\s*', '', summary_text, flags=re.IGNORECASE)

                            # Store the extracted sections (remove "Title:" prefix if present)
                            title_text = title_match.group(1).strip()
                            if title_text.startswith("Title:"):
                                title_text = title_text[6:].strip()
                            self.dataset_info[analysis_id]["title"] = title_text
                            self.dataset_info[analysis_id]["summary"] = summary_text

                            if design_match:
                                design_text = design_match.group(1).strip()
                                self.dataset_info[analysis_id]["design"] = design_text
                        else:
                            # Fallback to simple parsing if regex approach fails
                            lines = content.split('\n')
                            title_lines = []
                            summary_lines = []
                            design_lines = []

                            section = "title"
                            for line in lines:
                                if section == "title":
                                    if line.strip() == "":
                                        section = "summary"
                                    else:
                                        title_lines.append(line)
                                elif section == "summary":
                                    if (line.lower().startswith("overall design:") or
                                        line.lower().startswith("overall design") or
                                        line.lower().strip() == "overall design"):
                                        section = "design"
                                        # Add the "Overall design:" part to design
                                        design_lines.append(line)
                                    else:
                                        summary_lines.append(line)
                                else:  # design section
                                    design_lines.append(line)

                            # Store the extracted sections
                            title_text = '\n'.join(title_lines).strip()
                            if title_text.startswith("Title:"):
                                title_text = title_text[6:].strip()
                            self.dataset_info[analysis_id]["title"] = title_text
                            self.dataset_info[analysis_id]["summary"] = '\n'.join(summary_lines).strip()
                            self.dataset_info[analysis_id]["design"] = '\n'.join(design_lines).strip()

                        logger.info(f"Loaded and parsed dataset info for {analysis_id}")
                        break
                    except Exception as e:
                        logger.warning(f"Error loading dataset info from {info_path}: {str(e)}")

        logger.info(f"Loaded dataset info for {len(self.dataset_info)} analyses")

        # Log contrast info for debugging
        if hasattr(self, "contrast_info") and self.contrast_info:
            logger.info(f"Available contrast IDs: {', '.join(self.contrast_info.keys())}")

    def create_integrated_report(self,
                               top_frequent: int = 20,
                               top_unique: int = 10,
                               plot_type: str = 'violin',
                               gene_list: List[str] = None,
                               max_genes: int = 100,
                               min_unique_per_contrast: int = 1,
                               p_value_threshold: float = None,
                               lfc_threshold: float = None,
                               max_contrasts_for_unique: int = 2,
                               hide_x_labels: bool = True,
                               output_dir: str = None):
        """
        Create a comprehensive integrated report with all visualizations.

        Parameters:
        -----------
        top_frequent : int
            Number of top frequently occurring genes to include
        top_unique : int
            Number of top unique genes to include per contrast
        plot_type : str
            Type of expression plot: 'violin', 'box', or 'both'
        gene_list : List[str], optional
            Specific list of genes to plot (overrides automatic selection)
        max_genes : int
            Maximum number of genes to include in visualizations
        min_unique_per_contrast : int
            Minimum number of unique genes to select from each contrast
        p_value_threshold : float, optional
            Override the default p-value threshold for gene selection
        lfc_threshold : float, optional
            Override the default log fold change threshold for gene selection
        max_contrasts_for_unique : int
            Maximum number of contrasts a gene can appear in to be considered unique
        hide_x_labels : bool
            Whether to hide x-axis labels in expression plots
        output_dir : str, optional
            Override the default output directory

        Returns:
        --------
        str
            Path to the output directory containing the report
        """
        # Override output directory if specified
        if output_dir is not None:
            self._output_dir_requested = output_dir
            self.output_dir = None

        # Now ensure the output directory exists
        self._ensure_output_dir()

        # Load data if not already done
        if not self.deg_data:
            self.load_data()

        # Identify important genes if not provided
        if gene_list is None:
            gene_list = list(self.identify_important_genes(
                top_frequent=top_frequent,
                top_unique=top_unique,
                max_contrasts_for_unique=max_contrasts_for_unique,
                min_unique_per_contrast=min_unique_per_contrast,
                # Use provided thresholds or default to 0.01 for automatic selection to highlight stronger signals
                p_value_threshold=p_value_threshold if p_value_threshold is not None else 0.01,
                lfc_threshold=lfc_threshold
            ))

            # Limit the number of genes to avoid performance issues
            if len(gene_list) > max_genes:
                logger.warning(f"Too many genes selected ({len(gene_list)}), limiting to {max_genes}")
                gene_list = gene_list[:max_genes]

        # Create LFC heatmap
        try:
            # Adjust max_genes for heatmap to allow for more genes
            max_heatmap_genes = min(100, len(gene_list))
            output_path = os.path.join(self.output_dir, "lfc_heatmap.html")
            lfc_heatmap = self.create_lfc_heatmap(
                genes=gene_list[:max_heatmap_genes],
                output_file=output_path
            )
        except Exception as e:
            logger.error(f"Error creating LFC heatmap: {str(e)}")
            lfc_heatmap = None

        # Create expression plots (first page)
        try:
            output_path = os.path.join(self.output_dir, f"expression_{plot_type}_plots.html")
            expr_plot = self.create_expression_plots(
                genes=gene_list[:30],  # Only show first page in main report
                plot_type=plot_type,
                output_file=output_path,
                hide_x_labels=hide_x_labels  # Use parameter value
            )
        except Exception as e:
            logger.error(f"Error creating expression plots: {str(e)}")
            expr_plot = None

        # If we have many genes, create additional pages (but limit to 2 additional pages to avoid long processing)
        if len(gene_list) > 30:
            # Maximum number of genes per page
            genes_per_page = 30
            num_pages = min(3, (len(gene_list) + genes_per_page - 1) // genes_per_page)

            # Create a directory for additional pages
            plots_dir = os.path.join(self.output_dir, "expression_plots")
            os.makedirs(plots_dir, exist_ok=True)

            # Create additional pages (limited to 2 more pages)
            for page in range(1, num_pages):
                start_idx = page * genes_per_page
                end_idx = min(start_idx + genes_per_page, len(gene_list))
                page_genes = gene_list[start_idx:end_idx]

                try:
                    output_path = os.path.join(plots_dir, f"expression_{plot_type}_page{page+1}.html")
                    self.create_expression_plots(
                        genes=page_genes,
                        plot_type=plot_type,
                        output_file=output_path,
                        hide_x_labels=hide_x_labels,  # Use parameter value
                        page_number=page+1  # Page number for this batch
                    )
                    logger.info(f"Created expression plot page {page+1}")
                except Exception as e:
                    logger.error(f"Error creating expression plots page {page+1}: {str(e)}")

        # Generate additional pages links if there are many genes
        additional_pages_html = ""
        if len(gene_list) > 30:
            page_links = []
            num_additional_pages = min(2, (len(gene_list) + 29) // 30 - 1)
            for i in range(num_additional_pages):
                page_num = i + 2
                page_links.append(f'<a href="expression_plots/expression_{plot_type}_page{page_num}.html" target="_blank" style="margin-right: 5px; font-size: 12px; padding: 5px 10px;">Page {page_num}</a>')
            additional_pages_html = f"""
                            <p><small>Additional pages:
                            {''.join(page_links)}
                            </small></p>"""

        # Generate contrast details HTML
        if not self.deg_data:
            contrast_details_html = '<p>No contrast information available.</p>'
        else:
            contrast_rows = []
            for analysis_id, contrasts in self.deg_data.items():
                for contrast_id in contrasts.keys():
                    contrast_rows.append(f'<tr><td>{analysis_id}</td><td>{contrast_id}</td><td>{self._get_contrast_description(analysis_id, contrast_id)}</td></tr>')
            contrast_details_html = f"""
                    <table>
                        <tr>
                            <th>Analysis</th>
                            <th>Contrast</th>
                            <th>Description</th>
                        </tr>
                        {''.join(contrast_rows)}
                    </table>"""

        # Generate HTML report
        html_content = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>UORCA Results Integration Report</title>
            <style>
                body {{ font-family: Arial, sans-serif; line-height: 1.6; margin: 0; padding: 20px; color: #333; }}
                h1, h2, h3 {{ color: #2c3e50; }}
                .container {{ max-width: 1200px; margin: 0 auto; }}
                .card {{ border: 1px solid #ddd; border-radius: 4px; padding: 20px; margin-bottom: 20px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
                .vis-link {{ display: inline-block; background: #3498db; color: white; padding: 10px 15px; text-decoration: none; border-radius: 4px; margin-top: 10px; }}
                .vis-link:hover {{ background: #2980b9; }}
                table {{ border-collapse: collapse; width: 100%; margin-top: 20px; }}
                th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
                th {{ background-color: #f2f2f2; }}
            </style>
        </head>
        <body>
            <div class="container">
                <h1>UORCA Results Integration Report</h1>
                <p>Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>

                <div class="card">
                    <h2>Visualizations</h2>
                    <p>The following interactive visualizations are available:</p>
                    <ul>
                        <li>
                            <strong>Log Fold Change Heatmap</strong>: Displays differential expression patterns across contrasts.
                            <br><a href="lfc_heatmap.html" class="vis-link" target="_blank">View Heatmap</a>
                        </li>
                        <li>
                            <strong>Gene Expression Plots</strong>: Shows expression levels (CPM) for selected genes across samples.
                            <br><a href="expression_{plot_type}_plots.html" class="vis-link" target="_blank">View Expression Plots</a>
                            {additional_pages_html}
                        </li>
                    </ul>
                </div>

                <div class="card">
                    <h2>Integration Summary</h2>
                    <p>This report integrates results from:</p>
                    <ul>
                        <li><strong>{len(self.deg_data)} datasets</strong> with a total of {sum(len(contrasts) for contrasts in self.deg_data.values())} contrasts</li>
                        <li><strong>{len(gene_list)} genes</strong> selected for visualization</li>
                    </ul>

                    <h3>Analyses Overview</h3>
                    <table>
                        <tr>
                            <th>Accession</th>
                            <th>Organism</th>
                            <th># Samples</th>
                            <th># Contrasts</th>
                        </tr>
                        {''.join(f'<tr><td>{info.get("accession", "Unknown")}</td><td>{info.get("organism", "Unknown")}</td><td>{info.get("number_of_samples", 0)}</td><td>{info.get("number_of_contrasts", 0)}</td></tr>' for analysis_id, info in self.analysis_info.items())}
                    </table>

                    <h3>Contrast Details</h3>
                    {contrast_details_html}
                </div>

                <div class="card">
                    <h2>Gene List</h2>
                    <table>
                        <tr>
                            <th>#</th>
                            <th>Gene</th>
                        </tr>
                        {''.join(f'<tr><td>{i+1}</td><td>{gene}</td></tr>' for i, gene in enumerate(gene_list))}
                    </table>
                </div>
            </div>
        </body>
        </html>
        """

        # Write HTML to file
        index_path = os.path.join(self.output_dir, "index.html")
        with open(index_path, 'w') as f:
            f.write(html_content)

        logger.info(f"Generated index HTML report at {index_path}")
        return self.output_dir


def main():
    """
    Main entry point for the script.

    Integrates results from multiple RNA-seq analyses and creates visualizations.
    """
    parser = argparse.ArgumentParser(
        description="Integrate and visualize results from multiple UORCA RNA-seq analyses"
    )
    parser.add_argument(
        "--results_dir",
        type=str,
        required=True,
        help="Directory containing analysis results"
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default=None,
        help="Directory where integrated results will be saved"
    )
    parser.add_argument(
        "--pvalue_threshold",
        type=float,
        default=0.05,
        help="P-value threshold for considering genes as differentially expressed (default: 0.05)"
    )
    parser.add_argument(
        "--lfc_threshold",
        type=float,
        default=1.0,
        help="Log fold change threshold for considering genes as differentially expressed (default: 1.0)"
    )
    parser.add_argument(
        "--top_frequent",
        type=int,
        default=20,
        help="Number of most frequently occurring genes to include (default: 20)"
    )
    parser.add_argument(
        "--top_unique",
        type=int,
        default=5,
        help="Number of unique genes with large fold changes to include per contrast (default: 10)"
    )
    parser.add_argument(
        "--max_genes",
        type=int,
        default=100,
        help="Maximum number of genes to include in visualizations (default: 100)"
    )
    parser.add_argument(
        "--plot_type",
        type=str,
        choices=["violin", "box", "both"],
        default="violin",
        help="Type of expression plot to generate (default: violin)"
    )
    parser.add_argument(
        "--gene_list",
        type=str,
        default=None,
        help="Path to text file with specific list of genes to plot (one gene per line)"
    )

    args = parser.parse_args()

    # Initialize the integrator
    integrator = ResultsIntegrator(
        results_dir=args.results_dir,
        output_dir=args.output_dir,
        pvalue_threshold=args.pvalue_threshold,
        lfc_threshold=args.lfc_threshold,
        top_n_genes=max(args.top_frequent, args.top_unique)
    )

    # Load gene list if provided
    gene_list = None
    if args.gene_list:
        try:
            with open(args.gene_list, 'r') as f:
                gene_list = [line.strip() for line in f if line.strip()]
            logger.info(f"Loaded {len(gene_list)} genes from {args.gene_list}")
        except Exception as e:
            logger.error(f"Error loading gene list from {args.gene_list}: {str(e)}")

    # Create the integrated report
    output_dir = integrator.create_integrated_report(
        top_frequent=args.top_frequent,
        top_unique=args.top_unique,
        plot_type=args.plot_type,
        gene_list=gene_list,
        max_genes=args.max_genes
    )

    print(f"\nResults integration complete!")
    print(f"Interactive report available at: {os.path.join(output_dir, 'index.html')}")


if __name__ == "__main__":
    main()

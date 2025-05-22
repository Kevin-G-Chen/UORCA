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
    level=logging.INFO,
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

        # Create output directory if not provided
        if output_dir is None:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            self.output_dir = os.path.join(self.results_dir, f"integrated_results_{timestamp}")
        else:
            self.output_dir = os.path.abspath(output_dir)

        # Create output directory
        os.makedirs(self.output_dir, exist_ok=True)

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
        logger.info(f"Output will be saved to: {self.output_dir}")

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
            info_file = os.path.join(folder, "analysis_info.json")

            if os.path.isfile(info_file):
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
        # Always load analysis info first
        self.load_analysis_info()

        # Find data files
        deg_files, cpm_files = self.find_data_files()

        # Load contrast information if available
        self.contrast_info = {}

        # Check multiple potential locations for the contrasts file
        contrast_file_locations = [
            os.path.join(self.results_dir, "contrasts.csv"),
            os.path.join(os.path.dirname(self.results_dir), "contrasts.csv"),
            os.path.join(self.results_dir, "metadata", "contrasts.csv")
        ]

        contrasts_file = None
        for file_path in contrast_file_locations:
            if os.path.isfile(file_path):
                contrasts_file = file_path
                logger.info(f"Found contrasts file at {file_path}")
                break

        if contrasts_file:
            try:
                contrasts_df = pd.read_csv(contrasts_file)
                for _, row in contrasts_df.iterrows():
                    if 'name' in row and 'description' in row:
                        self.contrast_info[row['name']] = {
                            'description': row['description'],
                            'expression': row.get('expression', '')
                        }
                logger.info(f"Loaded information for {len(self.contrast_info)} contrasts")
            except Exception as e:
                logger.warning(f"Error loading contrasts file: {str(e)}")

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

        logger.info(f"Data loading complete: {len(self.deg_data)} DEG datasets, {len(self.cpm_data)} CPM datasets")

    def identify_important_genes(self,
                              top_frequent: int = 20,
                              top_unique: int = 10,
                              max_contrasts_for_unique: int = 2) -> Set[str]:
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

                # Find p-value column
                p_value_col = None
                if 'adj.P.Val' in df.columns:
                    p_value_col = 'adj.P.Val'
                else:
                    for col in ['padj', 'FDR', 'q_value', 'PValue', 'pvalue', 'P.Value']:
                        if col in df.columns:
                            p_value_col = col
                            break

                # Find logFC column
                lfc_col = None
                if 'logFC' in df.columns:
                    lfc_col = 'logFC'
                else:
                    for col in ['log2FoldChange', 'log2FC', 'LogFC']:
                        if col in df.columns:
                            lfc_col = col
                            break

                if p_value_col is None or lfc_col is None:
                    logger.warning(f"Missing required columns in {contrast_key}. Available columns: {df.columns.tolist()}")
                    continue

                logger.info(f"Using columns for {contrast_key}: p-value={p_value_col}, logFC={lfc_col}")

                # Filter significant genes
                sig_df = df[(df[p_value_col] < self.pvalue_threshold) &
                           (abs(df[lfc_col]) > self.lfc_threshold)]

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

            important_genes.update(top_unique_genes)
            logger.info(f"Selected {len(top_unique_genes)} unique genes with large fold changes from {contrast_key}")

        logger.info(f"Identified {len(important_genes)} important genes in total")
        return important_genes

    def create_lfc_heatmap(self,
                         genes: List[str] = None,
                         contrasts: List[Tuple[str, str]] = None,
                         output_file: str = None) -> go.Figure:
        """
        Create an interactive heatmap of log fold changes for selected genes across contrasts.

        Parameters:
        -----------
        genes : List[str], optional
            List of genes to include in the heatmap. If None, important genes are identified.
        contrasts : List[Tuple[str, str]], optional
            List of (analysis_id, contrast_id) tuples to include. If None, all contrasts are used.
        output_file : str, optional
            Path to save the interactive HTML output. If None, a default path is used.

        Returns:
        --------
        plotly.graph_objects.Figure
        """
        # If genes not provided, identify important genes
        if genes is None:
            genes = list(self.identify_important_genes())

        # Prepare data for heatmap
        heatmap_data = []

        # If contrasts not specified, use all
        if contrasts is None:
            contrasts = []
            for analysis_id, contrast_dict in self.deg_data.items():
                for contrast_id in contrast_dict.keys():
                    contrasts.append((analysis_id, contrast_id))

        # Create contrast labels
        contrast_labels = [f"{a_id}_{c_id}" for a_id, c_id in contrasts]

        # Extract log fold change values for each gene across contrasts
        for gene in genes:
            row = {'Gene': gene}

            for (analysis_id, contrast_id), contrast_label in zip(contrasts, contrast_labels):
                if analysis_id in self.deg_data and contrast_id in self.deg_data[analysis_id]:
                    df = self.deg_data[analysis_id][contrast_id]

                    if 'Gene' not in df.columns:
                        continue

                    # Find log fold change column
                    lfc_col = None
                    for col in ['logFC', 'log2FoldChange', 'log2FC', 'LogFC']:
                        if col in df.columns:
                            lfc_col = col
                            break

                    if lfc_col is None:
                        continue

                    # Find p-value column
                    p_value_col = None
                    for col in ['adj.P.Val', 'padj', 'FDR', 'q_value', 'PValue', 'pvalue', 'P.Value']:
                        if col in df.columns:
                            p_value_col = col
                            break

                    # Find the gene and get its logFC, but set to 0 if not significant
                    gene_row = df[df['Gene'] == gene]
                    if not gene_row.empty:
                        if p_value_col is not None:
                            p_value = gene_row.iloc[0][p_value_col]
                            abs_lfc = abs(gene_row.iloc[0][lfc_col])
                            if p_value >= self.pvalue_threshold or abs_lfc <= self.lfc_threshold:
                                row[contrast_label] = 0  # Not significant, set to white
                            else:
                                row[contrast_label] = gene_row.iloc[0][lfc_col]
                        else:
                            row[contrast_label] = gene_row.iloc[0][lfc_col]
                    else:
                        # Gene not found in this contrast
                        row[contrast_label] = 0  # Use 0 for non-significant/missing

            # Only add the row if there's at least one logFC value
            if len(row) > 1:  # More than just the Gene column
                heatmap_data.append(row)

        # Convert to DataFrame
        heatmap_df = pd.DataFrame(heatmap_data)

        if len(heatmap_df) == 0 or len(heatmap_df.columns) <= 1:
            logger.warning("Not enough data for heatmap visualization")
            return None

        # Prepare for heatmap - set Gene as index and fill NAs
        heatmap_df = heatmap_df.set_index('Gene')
        heatmap_df = heatmap_df.fillna(0)

        # Perform hierarchical clustering on both axes if there are enough data points
        # 1. Cluster genes (rows)
        if len(heatmap_df) > 1:  # Only cluster if we have multiple genes
            try:
                gene_linkage = sch.linkage(pdist(heatmap_df.values), method='complete')
                gene_order = sch.leaves_list(gene_linkage)
                clustered_genes = heatmap_df.index[gene_order].tolist()
            except Exception as e:
                logger.warning(f"Error clustering genes: {str(e)}. Using original order.")
                clustered_genes = heatmap_df.index.tolist()
        else:
            clustered_genes = heatmap_df.index.tolist()

        # 2. Cluster contrasts (columns)
        if len(heatmap_df.columns) > 1:  # Only cluster if we have multiple contrasts
            try:
                contrast_linkage = sch.linkage(pdist(heatmap_df.values.T), method='complete')
                contrast_order = sch.leaves_list(contrast_linkage)
                clustered_contrasts = [heatmap_df.columns[i] for i in contrast_order]
            except Exception as e:
                logger.warning(f"Error clustering contrasts: {str(e)}. Using original order.")
                clustered_contrasts = heatmap_df.columns.tolist()
        else:
            clustered_contrasts = heatmap_df.columns.tolist()

        # Reorder the DataFrame according to clustering
        heatmap_df = heatmap_df.loc[clustered_genes, clustered_contrasts]

        # Create heatmap using plotly
        max_abs_value = max(abs(heatmap_df.values.min()), abs(heatmap_df.values.max()))
        if max_abs_value == 0:
            max_abs_value = 1.0  # Avoid division by zero

        # Create the heatmap
        fig = px.imshow(
            heatmap_df,
            color_continuous_scale='RdBu_r',  # Red-Blue color scale, reversed to match biology convention
            labels=dict(x="Contrast", y="Gene", color="Log2FC"),
            title="Differential Expression Heatmap (Log2 Fold Change)",
            aspect="auto",  # Maintain aspect ratio based on data
            height=max(400, len(heatmap_df) * 15),  # Dynamic height based on number of genes
            width=max(600, len(heatmap_df.columns) * 70),  # Dynamic width based on number of contrasts
            zmin=-max_abs_value,  # Symmetrical color scale
            zmax=max_abs_value,
            color_continuous_midpoint=0  # Ensure white is at 0
        )

        # Improve layout
        fig.update_layout(
            margin=dict(l=150, r=20, t=50, b=100),
            coloraxis_colorbar=dict(title="Log2FC"),
            xaxis=dict(tickangle=45)
        )

        # Update hover template
        hover_template = (
            "<b>Gene:</b> %{y}<br>" +
            "<b>Contrast:</b> %{x}<br>" +
            "<b>Log2FC:</b> %{z:.2f}<br>"
        )
        fig.update_traces(hovertemplate=hover_template)

        # Save to file if specified
        if output_file is None:
            output_file = os.path.join(self.output_dir, "lfc_heatmap.html")

        fig.write_html(output_file)
        logger.info(f"Saved LFC heatmap to {output_file}")

        return fig

    def create_expression_plots(self,
                              genes: List[str],
                              plot_type: str = 'violin',
                              analyses: List[str] = None,
                              output_file: str = None) -> go.Figure:
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
            Path to save the interactive HTML output. If None, a default path is used.

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
        gene_list = genes.copy()
        gene_pages = [genes[i:i + genes_per_page] for i in range(0, len(genes), genes_per_page)]

        if len(gene_pages) > 1:
            logger.info(f"Splitting {len(genes)} genes into {len(gene_pages)} pages for visualization")
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

            # Check multiple potential locations for the sample mapping file
            sample_file_locations = [
                os.path.join(self.results_dir, "edger_analysis_samples.csv"),
                os.path.join(os.path.dirname(self.results_dir), "edger_analysis_samples.csv")
            ]

            edger_samples_file = None
            for file_path in sample_file_locations:
                if os.path.isfile(file_path):
                    edger_samples_file = file_path
                    logger.info(f"Found sample mapping file at {file_path}")
                    break

            if edger_samples_file:
                try:
                    # Load the sample mapping file - try both with and without index
                    try:
                        sample_mapping_df = pd.read_csv(edger_samples_file)
                        # Check if first column should be index
                        if sample_mapping_df.columns[0] == '' or sample_mapping_df.columns[0].startswith('Unnamed'):
                            sample_mapping_df = pd.read_csv(edger_samples_file, index_col=0)
                    except Exception as e:
                        logger.warning(f"Error reading sample mapping file, trying alternative approach: {str(e)}")
                        sample_mapping_df = pd.read_csv(edger_samples_file, index_col=0)

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
                        # Create a mapping from sample index to group
                        for i, row in enumerate(sample_mapping_df.itertuples()):
                            # Use sample index (i) as the key (starting with Sample1)
                            sample_key = f"Sample{i+1}"
                            if hasattr(row, group_col):
                                sample_to_group[sample_key] = getattr(row, group_col)
                            elif group_col in sample_mapping_df.columns:
                                sample_to_group[sample_key] = sample_mapping_df.iloc[i][group_col]
                        logger.info(f"Created mapping for {len(sample_to_group)} samples")
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

        # Determine facet column wrapping (3 columns max)
        facet_col_wrap = min(3, genes_found)

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

        # Calculate appropriate facet_row_spacing based on number of rows
        # Must be less than 1/(rows-1) to avoid the error
        if n_rows <= 1:
            facet_row_spacing = 0.05  # Default for single row
        else:
            # Use 80% of the maximum allowed value to be safe
            max_spacing = 1 / (n_rows - 1)
            facet_row_spacing = max_spacing * 0.8
            logger.info(f"Using facet_row_spacing of {facet_row_spacing:.4f} for {n_rows} rows (max allowed: {max_spacing:.4f})")

        # Create the plot using Plotly Express with adjusted spacing
        if plot_type in ['violin', 'both']:
            fig = px.violin(
                plot_df,
                x='Group',
                y='LogExpression',
                color='Analysis',
                facet_col='Gene',
                facet_col_wrap=facet_col_wrap,
                hover_data=['Sample', 'Expression'],
                box=True,
                points="all",
                title=f"Gene Expression Across Samples (Page 1 of {len(gene_pages)})",
                labels={'LogExpression': y_axis_title, 'Group': 'Group', 'Analysis': 'Dataset'},
                facet_row_spacing=facet_row_spacing
            )
        else:  # box plot
            fig = px.box(
                plot_df,
                x='Group',
                y='LogExpression',
                color='Analysis',
                facet_col='Gene',
                facet_col_wrap=facet_col_wrap,
                hover_data=['Sample', 'Expression'],
                points="all",
                title=f"Gene Expression Across Samples (Page 1 of {len(gene_pages)})",
                labels={'LogExpression': y_axis_title, 'Group': 'Group', 'Analysis': 'Dataset'},
                facet_row_spacing=facet_row_spacing
            )

        # Layout adjustments
        fig.update_layout(
            height=min(250 * n_rows, 2000),  # Cap height (reduced from 300 to 250 per row)
            width=min(400 * facet_col_wrap, 1200),  # Cap width
            legend_title_text="Dataset",
            margin=dict(l=40, r=20, t=80, b=40)
        )

        # Update facet formatting (remove "Gene=" prefix)
        fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[1]))

        # Add note about number of genes
        if len(genes) < len(gene_list):
            fig.add_annotation(
                text=f"Showing {genes_found} genes (page 1 of {len(gene_pages)})",
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

        # Save to file
        if output_file is None:
            output_file = os.path.join(self.output_dir, f"expression_{plot_type}_plots.html")

        fig.write_html(output_file)
        logger.info(f"Saved expression plots to {output_file}")

        return fig

    def _get_contrast_description(self, analysis_id: str, contrast_id: str) -> str:
        """Get the description for a contrast from the loaded contrast info."""
        if contrast_id in self.contrast_info:
            return self.contrast_info[contrast_id]['description']

        return f"Contrast: {contrast_id}"

    def create_integrated_report(self,
                               top_frequent: int = 20,
                               top_unique: int = 10,
                               plot_type: str = 'violin',
                               gene_list: List[str] = None,
                               max_genes: int = 100):
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

        Returns:
        --------
        str
            Path to the output directory containing the report
        """
        # Load data if not already done
        if not self.deg_data:
            self.load_data()

        # Identify important genes if not provided
        if gene_list is None:
            gene_list = list(self.identify_important_genes(
                top_frequent=top_frequent,
                top_unique=top_unique
            ))

            # Limit the number of genes to avoid performance issues
            if len(gene_list) > max_genes:
                logger.warning(f"Too many genes selected ({len(gene_list)}), limiting to {max_genes}")
                gene_list = gene_list[:max_genes]

        # Create LFC heatmap
        try:
            lfc_heatmap = self.create_lfc_heatmap(
                genes=gene_list,
                output_file=os.path.join(self.output_dir, "lfc_heatmap.html")
            )
        except Exception as e:
            logger.error(f"Error creating LFC heatmap: {str(e)}")
            lfc_heatmap = None

        # Create expression plots (first page)
        try:
            expr_plot = self.create_expression_plots(
                genes=gene_list[:30],  # Only show first page in main report
                plot_type=plot_type,
                output_file=os.path.join(self.output_dir, f"expression_{plot_type}_plots.html")
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
                    self.create_expression_plots(
                        genes=page_genes,
                        plot_type=plot_type,
                        output_file=os.path.join(plots_dir, f"expression_{plot_type}_page{page+1}.html")
                    )
                    logger.info(f"Created expression plot page {page+1}")
                except Exception as e:
                    logger.error(f"Error creating expression plots page {page+1}: {str(e)}")

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
                            {f"""
                            <p><small>Additional pages:
                            {''.join([f'<a href="expression_plots/expression_{plot_type}_page{i+2}.html" target="_blank" style="margin-right: 5px; font-size: 12px; padding: 5px 10px;">Page {i+2}</a>' for i in range(min(2, (len(gene_list) + 29) // 30 - 1))])}
                            </small></p>
                            """ if len(gene_list) > 30 else ''}
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
                            <th>Analysis Column</th>
                            <th># Samples</th>
                            <th># Contrasts</th>
                        </tr>
                        {''.join(f'<tr><td>{info.get("accession", "Unknown")}</td><td>{info.get("organism", "Unknown")}</td><td>{info.get("analysis_column", "Unknown")}</td><td>{info.get("number_of_samples", 0)}</td><td>{info.get("number_of_contrasts", 0)}</td></tr>' for analysis_id, info in self.analysis_info.items())}
                    </table>

                    <h3>Contrast Details</h3>
                    {'<p>No contrast information available.</p>' if not self.deg_data else f"""
                    <table>
                        <tr>
                            <th>Analysis</th>
                            <th>Contrast</th>
                            <th>Description</th>
                        </tr>
                        {''.join([f'<tr><td>{analysis_id}</td><td>{contrast_id}</td><td>{self._get_contrast_description(analysis_id, contrast_id)}</td></tr>' for analysis_id, contrasts in self.deg_data.items() for contrast_id in contrasts.keys()])}
                    </table>
                    """}
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

#!/usr/bin/env python3
"""
ResultsIntegration.py - Integrates results from multiple UORCA RNA-seq analyses

This script collates differential expression results and normalized expression data
across multiple experiments and contrasts, producing interactive visualizations
that highlight important genes and patterns.

The main outputs are:
1. Interactive heatmap of important genes with LFC values (clustered on both axes)
2. Interactive boxplot/violin plots of gene expression across samples

These visualizations are combined into an HTML report that can be easily shared.
"""

import os
import glob
import argparse
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import pdist
from typing import List, Dict, Tuple, Optional, Union, Set
import logging
import re
from pathlib import Path
from datetime import datetime
import json

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
            Directory containing multiple analysis results (with subdirectories for each analysis)
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
        self.analysis_info = {}  # {analysis_id: dict with metadata}
        self.deg_files = {}      # {analysis_id: {contrast_id: path}}
        self.cpm_files = {}      # {analysis_id: path}
        self.metadata_files = {}  # {analysis_id: path}
        self.deg_data = {}       # {analysis_id: {contrast_id: DataFrame}}
        self.cpm_data = {}       # {analysis_id: DataFrame}
        self.metadata = {}       # {analysis_id: DataFrame}

        # Mapping to maintain consistent colors
        self.color_map = {}

        logger.info(f"Initialized ResultsIntegrator with results_dir: {self.results_dir}")
        logger.info(f"Output will be saved to: {self.output_dir}")

    def find_analysis_folders(self) -> List[str]:
        """Find all analysis folders in the results directory."""
        # Look for directories that contain RNAseqAnalysis folders
        analysis_folders = []
        for item in os.listdir(self.results_dir):
            item_path = os.path.join(self.results_dir, item)

            # Check if it's a directory
            if os.path.isdir(item_path):
                # Check if it contains an RNAseqAnalysis subdirectory
                rnaseq_dir = os.path.join(item_path, "RNAseqAnalysis")
                if os.path.isdir(rnaseq_dir):
                    analysis_folders.append(item_path)

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

    def find_data_files(self) -> Tuple[Dict, Dict, Dict]:
        """
        Find all DEG, CPM, and metadata files across all analyses.

        Returns:
        --------
        Tuple of (deg_files, cpm_files, metadata_files) dictionaries
        where keys are analysis IDs and values are paths or dicts of paths
        """
        analysis_folders = self.find_analysis_folders()

        deg_files = {}
        cpm_files = {}
        metadata_files = {}

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

            # Find metadata file
            metadata_file = os.path.join(folder, "metadata", "meta_long.csv")
            if os.path.isfile(metadata_file):
                metadata_files[analysis_id] = metadata_file

        logger.info(f"Found DEG files for {len(deg_files)} analyses")
        logger.info(f"Found CPM files for {len(cpm_files)} analyses")
        logger.info(f"Found metadata files for {len(metadata_files)} analyses")

        self.deg_files = deg_files
        self.cpm_files = cpm_files
        self.metadata_files = metadata_files

        return deg_files, cpm_files, metadata_files

    def load_data(self):
        """
        Load all DEG, CPM, and metadata files into memory.
        """
        # Load analysis info if not already done
        if not self.analysis_info:
            self.load_analysis_info()

        # Find data files if not already done
        if not self.deg_files:
            self.find_data_files()

        # Load DEG data
        for analysis_id, contrasts in self.deg_files.items():
            self.deg_data[analysis_id] = {}

            for contrast_id, deg_file in contrasts.items():
                # Load DEG file
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
        for analysis_id, cpm_file in self.cpm_files.items():
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

        # Load metadata
        for analysis_id, metadata_file in self.metadata_files.items():
            try:
                df = pd.read_csv(metadata_file)
                self.metadata[analysis_id] = df
                logger.debug(f"Loaded metadata for {analysis_id} with {len(df)} rows")
            except Exception as e:
                logger.error(f"Error loading metadata file {metadata_file}: {str(e)}")

        logger.info(f"Data loading complete: {len(self.deg_data)} DEG datasets, {len(self.cpm_data)} CPM datasets, {len(self.metadata)} metadata files")

    def identify_important_genes(self,
                                method: str = 'top_by_occurrence',
                                n_genes: int = None) -> Set[str]:
        """
        Identify important genes across all analyses based on specified method.

        Parameters:
        -----------
        method : str
            Method to identify important genes:
            - 'top_by_contrast': Top n genes by abs(logFC) from each contrast
            - 'top_by_occurrence': Top n genes that appear in the most contrasts
            - 'combined': Union of the above methods
        n_genes : int, optional
            Number of top genes to select (defaults to self.top_n_genes)

        Returns:
        --------
        Set of important gene identifiers
        """
        if n_genes is None:
            n_genes = self.top_n_genes

        important_genes = set()
        genes_by_occurrence = {}

        # Collect significant genes and their occurrence count
        for analysis_id, contrasts in self.deg_data.items():
            for contrast_id, df in contrasts.items():
                # Skip if Gene column doesn't exist
                if 'Gene' not in df.columns:
                    logger.warning(f"No 'Gene' column in {analysis_id}/{contrast_id}")
                    continue

                # Find differentially expressed genes
                # Check if we have p-value and logFC columns
                p_value_col = None
                for col in ['adj.P.Val', 'padj', 'FDR', 'q_value', 'PValue', 'pvalue', 'P.Value']:
                    if col in df.columns:
                        p_value_col = col
                        break

                lfc_col = None
                for col in ['logFC', 'log2FoldChange', 'log2FC', 'LogFC']:
                    if col in df.columns:
                        lfc_col = col
                        break

                if p_value_col is None or lfc_col is None:
                    logger.warning(f"Missing p-value or logFC column in {analysis_id}/{contrast_id}")
                    continue

                # Filter significant genes
                sig_genes = df[(df[p_value_col] < self.pvalue_threshold) &
                              (abs(df[lfc_col]) > self.lfc_threshold)]['Gene'].tolist()

                # Count occurrences of each gene
                for gene in sig_genes:
                    if gene not in genes_by_occurrence:
                        genes_by_occurrence[gene] = 0
                    genes_by_occurrence[gene] += 1

                # If using top_by_contrast method, add top genes from this contrast
                if method in ['top_by_contrast', 'combined']:
                    # Sort by absolute logFC and take top n
                    top_genes = df.sort_values(by=lfc_col, key=abs, ascending=False).head(n_genes)['Gene'].tolist()
                    important_genes.update(top_genes)

        # If using top_by_occurrence method, add genes that appear in the most contrasts
        if method in ['top_by_occurrence', 'combined']:
            # Sort genes by occurrence count (descending) and take top n
            sorted_genes = sorted(genes_by_occurrence.items(), key=lambda x: x[1], reverse=True)
            top_occurrence_genes = [gene for gene, count in sorted_genes[:n_genes]]
            important_genes.update(top_occurrence_genes)

        logger.info(f"Identified {len(important_genes)} important genes using method '{method}'")
        return important_genes

    def create_lfc_heatmap(self,
                          genes: List[str] = None,
                          contrasts: List[Tuple[str, str]] = None,
                          output_file: str = None) -> go.Figure:
        """
        Create an interactive heatmap of log fold changes for selected genes across contrasts.
        The heatmap is clustered on both axes using hierarchical clustering.

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
            genes = list(self.identify_important_genes(method='top_by_occurrence'))

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

                    # Find the gene and get its logFC
                    gene_row = df[df['Gene'] == gene]
                    if not gene_row.empty:
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

        # Perform hierarchical clustering on both axes
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
        fig = px.imshow(
            heatmap_df,
            color_continuous_scale='RdBu_r',  # Red-Blue color scale, reversed to match biology convention
            labels=dict(x="Contrast", y="Gene", color="Log2FC"),
            title="Differential Expression Heatmap (Log2 Fold Change)",
            aspect="auto",  # Maintain aspect ratio based on data
            height=max(400, len(heatmap_df) * 15),  # Dynamic height based on number of genes
            width=max(600, len(heatmap_df.columns) * 70)  # Dynamic width based on number of contrasts
        )

        # Improve layout
        fig.update_layout(
            margin=dict(l=150, r=20, t=50, b=100),
            coloraxis_colorbar=dict(title="Log2FC"),
            xaxis=dict(tickangle=45)
        )

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

        # Prepare data for plots
        plot_data = []

        for analysis_id in analyses:
            if analysis_id not in self.cpm_data:
                continue

            cpm_df = self.cpm_data[analysis_id]

            # Skip if Gene column doesn't exist
            if 'Gene' not in cpm_df.columns:
                logger.warning(f"No 'Gene' column in CPM data for {analysis_id}")
                continue

            # Use analysis_info to get grouping variable if available
            group_col = None
            if analysis_id in self.analysis_info and 'analysis_column' in self.analysis_info[analysis_id]:
                group_col = self.analysis_info[analysis_id]['analysis_column']
                logger.info(f"Using '{group_col}' as grouping variable for {analysis_id} from analysis_info")

                # Create sample to group mapping using analysis_info and metadata
                sample_to_group = {}

                # Check if we have metadata for this analysis
                if analysis_id in self.metadata:
                    metadata_df = self.metadata[analysis_id]

                    # Find columns that might contain sample IDs
                    sample_id_cols = []
                    for col in metadata_df.columns:
                        if col.lower() in ['gsm', 'sample', 'sample_id', 'geo_accession']:
                            sample_id_cols.append(col)

                    # Create mapping if we have a group column and sample ID columns
                    if sample_id_cols and group_col in metadata_df.columns:
                        for _, row in metadata_df.iterrows():
                            for col in sample_id_cols:
                                if pd.notna(row[col]):
                                    sample_id = str(row[col])
                                    group_val = row[group_col]
                                    sample_to_group[sample_id] = group_val

                    # Extract expression data for each gene
                    for gene in genes:
                        gene_rows = cpm_df[cpm_df['Gene'] == gene]

                        if gene_rows.empty:
                            logger.warning(f"Gene {gene} not found in CPM data for {analysis_id}")
                            continue

                        gene_row = gene_rows.iloc[0]

                        # Process each sample column
                        for col in cpm_df.columns:
                            if col != 'Gene':
                                # Column name is the sample ID
                                sample_id = col
                                expr_value = gene_row[col]

                                # Try to find the group for this sample
                                group_value = None

                                # Direct match
                                if sample_id in sample_to_group:
                                    group_value = sample_to_group[sample_id]
                                else:
                                    # Try partial matches
                                    for s_id, grp in sample_to_group.items():
                                        if s_id in sample_id or sample_id in s_id:
                                            group_value = grp
                                            break

                                # Use "Unknown" as fallback
                                if group_value is None:
                                    group_value = "Unknown"

                                plot_data.append({
                                    'Gene': gene,
                                    'Expression': expr_value,
                                    'Sample': sample_id,
                                    'Group': group_value,
                                    'Analysis': analysis_id
                                })
                else:
                    logger.warning(f"No metadata for {analysis_id}, cannot map samples to groups")
            else:
                logger.warning(f"No analysis_column in analysis_info for {analysis_id}, skipping")

        if not plot_data:
            logger.warning("No expression data found for specified genes and analyses")
            return None

        # Convert to DataFrame
        plot_df = pd.DataFrame(plot_data)

        # Create subplot grid based on number of genes
        n_genes = len(plot_df['Gene'].unique())
        n_cols = min(3, n_genes)
        n_rows = (n_genes + n_cols - 1) // n_cols  # Ceiling division

        fig = make_subplots(
            rows=n_rows,
            cols=n_cols,
            subplot_titles=[f"Gene: {gene}" for gene in plot_df['Gene'].unique()],
            vertical_spacing=0.1
        )

        # Add traces for each gene
        for i, gene in enumerate(plot_df['Gene'].unique()):
            gene_data = plot_df[plot_df['Gene'] == gene]

            row = i // n_cols + 1
            col = i % n_cols + 1

            # Add traces based on plot type
            if plot_type in ['violin', 'both']:
                for analysis in gene_data['Analysis'].unique():
                    analysis_data = gene_data[gene_data['Analysis'] == analysis]

                    fig.add_trace(
                        go.Violin(
                            x=analysis_data['Group'],
                            y=analysis_data['Expression'],
                            name=analysis,
                            box_visible=True,
                            meanline_visible=True,
                            points="all",
                            legendgroup=analysis,
                            showlegend=i==0,  # Only show in legend for first gene
                            line_color=self._get_color(analysis)
                        ),
                        row=row, col=col
                    )

            elif plot_type == 'box':
                for analysis in gene_data['Analysis'].unique():
                    analysis_data = gene_data[gene_data['Analysis'] == analysis]

                    fig.add_trace(
                        go.Box(
                            x=analysis_data['Group'],
                            y=analysis_data['Expression'],
                            name=analysis,
                            legendgroup=analysis,
                            showlegend=i==0,  # Only show in legend for first gene
                            marker_color=self._get_color(analysis)
                        ),
                        row=row, col=col
                    )

        # Update layout
        fig.update_layout(
            title="Gene Expression (CPM) Across Samples",
            height=300 * n_rows,
            width=400 * n_cols,
            boxmode='group',
            violinmode='group',
            margin=dict(l=50, r=50, t=50, b=50)
        )

        # Update y-axis title
        fig.update_yaxes(title_text="Expression (CPM)")

        # Save to file if specified
        if output_file is None:
            output_file = os.path.join(self.output_dir, f"expression_{plot_type}_plots.html")

        fig.write_html(output_file)
        logger.info(f"Saved expression plots to {output_file}")

        return fig

    def _get_color(self, key: str) -> str:
        """Get a consistent color for a given key."""
        if key not in self.color_map:
            # Default plotly colors
            colors = px.colors.qualitative.Plotly
            self.color_map[key] = colors[len(self.color_map) % len(colors)]
        return self.color_map[key]

    def create_integrated_report(self,
                               top_n_genes: int = None,
                               plot_type: str = 'violin',
                               contrasts_subset: List[Tuple[str, str]] = None,
                               analyses_subset: List[str] = None,
                               gene_list: List[str] = None):
        """
        Create a comprehensive integrated report with all visualizations.

        Parameters:
        -----------
        top_n_genes : int, optional
            Number of top genes to include (default from class initialization)
        plot_type : str
            Type of expression plot: 'violin', 'box', or 'both'
        contrasts_subset : List[Tuple[str, str]], optional
            List of (analysis_id, contrast_id) tuples to include in heatmap
        analyses_subset : List[str], optional
            List of analysis IDs to include in expression plots
        gene_list : List[str], optional
            Specific list of genes to plot (overrides automatic selection)
        """
        # Load data if not already done
        if not self.deg_data:
            self.load_data()

        if top_n_genes is None:
            top_n_genes = self.top_n_genes

        # Identify important genes if not provided
        if gene_list is None:
            gene_list = list(self.identify_important_genes(
                method='top_by_occurrence',
                n_genes=top_n_genes
            ))

        # Create LFC heatmap
        lfc_heatmap = self.create_lfc_heatmap(
            genes=gene_list,
            contrasts=contrasts_subset,
            output_file=os.path.join(self.output_dir, "lfc_heatmap.html")
        )

        # Create expression plots
        expr_plot = self.create_expression_plots(
            genes=gene_list,
            plot_type=plot_type,
            analyses=analyses_subset,
            output_file=os.path.join(self.output_dir, f"expression_{plot_type}_plots.html")
        )

        # Generate a summary HTML report
        self._generate_html_report(gene_list)

        logger.info(f"Created integrated report in {self.output_dir}")
        return self.output_dir

    def _generate_html_report(self, gene_list: List[str]):
        """Generate an HTML index page linking all visualizations."""
        # Calculate some summary statistics for the report
        total_contrasts = sum(len(contrasts) for contrasts in self.deg_data.values())

        # Create a dictionary to store genes and their occurrence counts
        gene_occurrence = {}
        for analysis_id, contrasts in self.deg_data.items():
            for contrast_id, df in contrasts.items():
                # Skip if Gene column doesn't exist
                if 'Gene' not in df.columns:
                    continue

                # Find p-value and logFC columns
                p_value_col = None
                for col in ['adj.P.Val', 'padj', 'FDR', 'q_value', 'PValue', 'pvalue', 'P.Value']:
                    if col in df.columns:
                        p_value_col = col
                        break

                lfc_col = None
                for col in ['logFC', 'log2FoldChange', 'log2FC', 'LogFC']:
                    if col in df.columns:
                        lfc_col = col
                        break

                if p_value_col is None or lfc_col is None:
                    continue

                # Filter significant genes
                sig_genes = df[(df[p_value_col] < self.pvalue_threshold) &
                              (abs(df[lfc_col]) > self.lfc_threshold)]['Gene'].tolist()

                # Count occurrences
                for gene in sig_genes:
                    if gene not in gene_occurrence:
                        gene_occurrence[gene] = 0
                    gene_occurrence[gene] += 1

        # Get top genes by occurrence for the report
        top_occurring_genes = []
        if gene_occurrence:
            sorted_genes = sorted(gene_occurrence.items(), key=lambda x: x[1], reverse=True)
            top_occurring_genes = [(gene, count) for gene, count in sorted_genes[:min(20, len(sorted_genes))]]

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
                tr:nth-child(even) {{ background-color: #f9f9f9; }}
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
                        </li>
                    </ul>
                </div>

                <div class="card">
                    <h2>Integration Summary</h2>
                    <p>This report integrates results from:</p>
                    <ul>
                        <li><strong>{len(self.deg_data)} datasets</strong> with a total of <strong>{total_contrasts} contrasts</strong></li>
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
                </div>

                <div class="card">
                    <h2>Top Differentially Expressed Genes</h2>
                    <p>Genes appearing as differentially expressed in the most contrasts:</p>
                    <table>
                        <tr>
                            <th>Rank</th>
                            <th>Gene</th>
                            <th># Contrasts</th>
                        </tr>
                        {''.join(f'<tr><td>{i+1}</td><td>{gene}</td><td>{count}</td></tr>' for i, (gene, count) in enumerate(top_occurring_genes))}
                    </table>
                </div>

                <div class="card">
                    <h2>Visualized Genes</h2>
                    <p>The following {len(gene_list)} genes were included in the visualizations:</p>
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


def main():
    """Main entry point for the script."""
    parser = argparse.ArgumentParser(
        description="Integrate and visualize results from multiple UORCA RNA-seq analyses"
    )
    parser.add_argument(
        "--results_dir",
        type=str,
        required=True,
        help="Directory containing multiple analysis results (with subdirectories for each analysis)"
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default=None,
        help="Directory where integrated results will be saved (default: results_dir/integrated_results_TIMESTAMP)"
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
        "--top_n_genes",
        type=int,
        default=50,
        help="Number of top genes to include in visualizations (default: 50)"
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
    parser.add_argument(
        "--contrasts_filter",
        type=str,
        default=None,
        help="Comma-separated list of analysis_id:contrast_id pairs to include (e.g., 'GSE123:A_vs_B,GSE456:C_vs_D')"
    )
    parser.add_argument(
        "--analysis_summary",
        action="store_true",
        help="Generate a separate summary of all analyses and their contrasts"
    )

    args = parser.parse_args()

    # Initialize and run the integrator
    integrator = ResultsIntegrator(
        results_dir=args.results_dir,
        output_dir=args.output_dir,
        pvalue_threshold=args.pvalue_threshold,
        lfc_threshold=args.lfc_threshold,
        top_n_genes=args.top_n_genes
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

    # Parse contrasts filter if provided
    contrasts_subset = None
    if args.contrasts_filter:
        try:
            contrasts_subset = []
            for pair in args.contrasts_filter.split(','):
                analysis_id, contrast_id = pair.strip().split(':')
                contrasts_subset.append((analysis_id, contrast_id))
            logger.info(f"Filtering to {len(contrasts_subset)} specified contrasts")
        except Exception as e:
            logger.error(f"Error parsing contrasts filter: {str(e)}")

    # Create the integrated report
    output_dir = integrator.create_integrated_report(
        top_n_genes=args.top_n_genes,
        plot_type=args.plot_type,
        contrasts_subset=contrasts_subset,
        gene_list=gene_list
    )

    # Generate analysis summary if requested
    if args.analysis_summary:
        # Create a summary of all analyses and their contrasts
        summary_data = []

        # Make sure data is loaded
        if not integrator.deg_data:
            integrator.load_data()

        # Collect summary information
        for analysis_id, contrasts in integrator.deg_data.items():
            analysis_info = integrator.analysis_info.get(analysis_id, {})

            for contrast_id in contrasts.keys():
                df = integrator.deg_data[analysis_id][contrast_id]

                # Count significant DEGs
                p_value_col = None
                for col in ['adj.P.Val', 'padj', 'FDR', 'q_value', 'PValue', 'pvalue', 'P.Value']:
                    if col in df.columns:
                        p_value_col = col
                        break

                lfc_col = None
                for col in ['logFC', 'log2FoldChange', 'log2FC', 'LogFC']:
                    if col in df.columns:
                        lfc_col = col
                        break

                if p_value_col is not None and lfc_col is not None:
                    sig_up = ((df[p_value_col] < args.pvalue_threshold) & (df[lfc_col] > args.lfc_threshold)).sum()
                    sig_down = ((df[p_value_col] < args.pvalue_threshold) & (df[lfc_col] < -args.lfc_threshold)).sum()
                else:
                    sig_up = sig_down = "N/A"

                summary_data.append({
                    "Analysis": analysis_id,
                    "Accession": analysis_info.get("accession", "Unknown"),
                    "Organism": analysis_info.get("organism", "Unknown"),
                    "Contrast": contrast_id,
                    "Total_DEGs": len(df),
                    "Up_DEGs": sig_up,
                    "Down_DEGs": sig_down
                })

        # Create a DataFrame and save to CSV
        if summary_data:
            summary_df = pd.DataFrame(summary_data)
            summary_file = os.path.join(integrator.output_dir, "analysis_summary.csv")
            summary_df.to_csv(summary_file, index=False)
            print(f"Analysis summary saved to: {summary_file}")

    print(f"\nResults integration complete!")
    print(f"Interactive report available at: {os.path.join(output_dir, 'index.html')}")


if __name__ == "__main__":
    main()

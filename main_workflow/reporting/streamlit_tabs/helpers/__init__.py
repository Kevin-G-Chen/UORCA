"""
Helper utilities for UORCA Explorer Streamlit tabs.

This module contains shared functions, classes, and utilities used across
multiple Streamlit tabs in the UORCA Explorer application.
"""

import os
import logging
import streamlit as st
import plotly.graph_objects as go
from pathlib import Path
from typing import Tuple, Optional, Dict, Any, Set, List
from datetime import datetime
import zipfile
import shutil
import io
import pandas as pd

# Import the main integrator
from ResultsIntegration import ResultsIntegrator

# Import streamlit logging utilities
from .streamlit_logging import (
    setup_streamlit_logging,
    log_streamlit_function,
    log_streamlit_agent,
    log_streamlit_tab,
    log_streamlit_event,
    log_streamlit_data_load,
    log_streamlit_user_action
)

# Import AI agent tool logging utilities
from .ai_agent_tool_logger import (
    start_ai_analysis_session,
    get_ai_tool_logs_for_display,
    clear_ai_tool_logs,
    get_ai_tool_logger,
    get_current_log_file,
    read_log_file_contents
)

logger = logging.getLogger(__name__)


class ModuleFilter(logging.Filter):
    """Filter for logging modules."""

    def __init__(self, names):
        super().__init__()
        self.names = names

    def filter(self, record):
        return any(record.name.startswith(n) for n in self.names)


def plotly_fig_to_pdf_bytes(fig: go.Figure, width: Optional[int] = None, height: Optional[int] = None, scale: float = 2.0) -> Optional[bytes]:
    """
    Convert a Plotly figure to PDF bytes for download.
    
    Args:
        fig: Plotly figure object
        width: Width of the PDF in pixels (defaults to figure width if None)
        height: Height of the PDF in pixels (defaults to figure height if None)
        scale: Scale factor for higher resolution (default 2.0 for 2x resolution)
    
    Returns:
        PDF bytes if successful, None if failed
    """
    try:
        import plotly.io as pio
        # Determine export dimensions from figure if not provided
        export_width = int(fig.layout.width) if width is None and fig.layout.width else width
        export_height = int(fig.layout.height) if height is None and fig.layout.height else height

        # Build args for to_image without passing None for width/height
        to_image_kwargs = dict(format='pdf', scale=scale)
        if export_width is not None:
            to_image_kwargs['width'] = export_width
        if export_height is not None:
            to_image_kwargs['height'] = export_height

        # Convert figure to PDF bytes using kaleido
        pdf_bytes = pio.to_image(fig, **to_image_kwargs)
        
        logger.debug(f"Successfully converted Plotly figure to PDF ({len(pdf_bytes)} bytes)")
        return pdf_bytes
        
    except ImportError as e:
        logger.error(f"PDF export requires kaleido package: {e}")
        st.error("PDF export requires the kaleido package. Please install it with: uv add kaleido")
        return None
    except Exception as e:
        logger.error(f"Failed to convert figure to PDF: {e}")
        st.error(f"Failed to generate PDF: {str(e)}")
        return None


def generate_plot_filename(plot_type: str, extension: str = "pdf") -> str:
    """
    Generate a timestamped filename for plot downloads.
    
    Args:
        plot_type: Type of plot (e.g., "heatmap", "expression_plots")
        extension: File extension (default "pdf")
    
    Returns:
        Formatted filename with timestamp
    """
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return f"uorca_{plot_type}_{timestamp}.{extension}"


def add_directory_to_zip(zip_file: zipfile.ZipFile, source_dir: str, archive_dir: str):
    """
    Recursively add a directory to a ZIP file.
    
    Args:
        zip_file: ZipFile object to add files to
        source_dir: Source directory path
        archive_dir: Directory name in the archive
    """
    if not os.path.exists(source_dir):
        logger.warning(f"Directory not found: {source_dir}")
        return
    
    for root, dirs, files in os.walk(source_dir):
        for file in files:
            file_path = os.path.join(root, file)
            arcname = os.path.join(archive_dir, os.path.relpath(file_path, source_dir))
            try:
                zip_file.write(file_path, arcname)
            except Exception as e:
                logger.warning(f"Could not add {file_path} to archive: {e}")


def create_dataset_download_package(
    ri: 'ResultsIntegrator',
    results_dir: str,
    dataset_id: str,
    dataset_accession: str,
    components: Optional[Dict[str, bool]] = None
) -> Optional[bytes]:
    """
    Create a comprehensive download package for a dataset analysis.
    
    Args:
        ri: ResultsIntegrator instance
        results_dir: Path to results directory
        dataset_id: Dataset ID (internal)
        dataset_accession: Dataset accession (e.g., GSE147834)
        components: Dictionary specifying which components to include:
                   - 'metadata': Include metadata files
                   - 'abundance': Include abundance data
                   - 'qc_plots': Include QC plots
                   - 'deg_analysis': Include DEG analysis data and plots
                   If None, all components are included (default)
    
    Returns:
        ZIP file bytes if successful, None if failed
    """
    # Default to including all components if not specified
    if components is None:
        components = {
            'metadata': True,
            'abundance': True,
            'qc_plots': True,
            'deg_analysis': True
        }
    
    try:
        # Create in-memory ZIP file
        zip_buffer = io.BytesIO()
        
        with zipfile.ZipFile(zip_buffer, mode='w', compression=zipfile.ZIP_DEFLATED) as zf:
            package_name = f"{dataset_accession}_analysis_package"
            dataset_path = os.path.join(results_dir, dataset_id)
            
            # 1. Add metadata files (if selected) - only GSE_metadata.csv and contrasts.csv
            if components.get('metadata', True):
                metadata_dir = os.path.join(dataset_path, "metadata")
                if os.path.exists(metadata_dir):
                    # Only include specific metadata files
                    metadata_files = [
                        f"{dataset_accession}_metadata.csv",
                        "contrasts.csv"
                    ]
                    for filename in metadata_files:
                        file_path = os.path.join(metadata_dir, filename)
                        if os.path.exists(file_path):
                            with open(file_path, 'rb') as f:
                                zf.writestr(f"{package_name}/metadata/{filename}", f.read())
            
            # 2. Add abundance directory (if selected and exists)
            if components.get('abundance', True):
                abundance_dir = os.path.join(dataset_path, "abundance")
                if os.path.exists(abundance_dir):
                    add_directory_to_zip(zf, abundance_dir, f"{package_name}/abundance")
            
            # 3. Add QC plots (if selected)
            if components.get('qc_plots', True):
                include_pdfs = components.get('include_pdfs', False)
                _add_qc_plots_to_zip(zf, ri, dataset_id, dataset_path, package_name, include_pdfs)
            
            # 4. Add DEG analysis data and plots (if selected)
            if components.get('deg_analysis', True):
                include_pdfs = components.get('include_pdfs', False)
                _add_deg_analysis_to_zip(zf, ri, dataset_id, dataset_path, package_name, include_pdfs)
                
                # Add CPM data as part of DEG analysis
                if dataset_id in ri.cpm_data:
                    cpm_df = ri.cpm_data[dataset_id]
                    cpm_csv = cpm_df.to_csv(index=False)
                    zf.writestr(f"{package_name}/DEG_analysis/CPM_data.csv", cpm_csv)
            
            # 6. Create and add README (always include, but customize content)
            readme_content = _generate_dataset_readme(ri, dataset_id, dataset_accession, components)
            zf.writestr(f"{package_name}/README.txt", readme_content)
        
        zip_buffer.seek(0)
        return zip_buffer.getvalue()
        
    except Exception as e:
        logger.error(f"Failed to create download package for {dataset_accession}: {e}")
        st.error(f"Failed to create download package: {str(e)}")
        return None


def _add_qc_plots_to_zip(
    zf: zipfile.ZipFile,
    ri: 'ResultsIntegrator',
    dataset_id: str,
    dataset_path: str,
    package_name: str,
    include_pdfs: bool = False
):
    """Add QC plots to the ZIP package."""
    rnaseq_dir = os.path.join(dataset_path, "RNAseqAnalysis")
    
    # Static QC plots (always include PNGs)
    static_plots = [
        ("MDS.png", "MDS.png"),
        ("filtering_density.png", "filtering_density.png"),
        ("normalization_boxplots.png", "normalization_boxplots.png"),
        ("voom_mean_variance.png", "voom_mean_variance.png"),
        ("sa_plot.png", "sa_plot.png")
    ]
    
    for source_file, dest_file in static_plots:
        source_path = os.path.join(rnaseq_dir, source_file)
        if os.path.exists(source_path):
            with open(source_path, 'rb') as f:
                zf.writestr(f"{package_name}/QC_plots/{dest_file}", f.read())
    
    # Only generate PCA plot PDF if requested and data available
    if include_pdfs and dataset_id in ri.cpm_data:
        try:
            # Import only if needed
            from single_analysis_plots import create_pca_plot, load_sample_groups
            
            cpm_df = ri.cpm_data[dataset_id]
            analysis_info = ri.analysis_info.get(dataset_id)
            groups = load_sample_groups(ri.results_dir, dataset_id, cpm_df, analysis_info)
            
            pca_fig = create_pca_plot(cpm_df, groups)
            if pca_fig:
                pdf_bytes = plotly_fig_to_pdf_bytes(pca_fig, width=1200, height=800, scale=2)
                if pdf_bytes:
                    zf.writestr(f"{package_name}/QC_plots/PCA_plot.pdf", pdf_bytes)
        except Exception as e:
            logger.debug(f"Could not generate PCA plot PDF: {e}")


def _add_deg_analysis_to_zip(
    zf: zipfile.ZipFile,
    ri: 'ResultsIntegrator',
    dataset_id: str,
    dataset_path: str,
    package_name: str,
    include_pdfs: bool = False
):
    """Add DEG analysis data and plots to the ZIP package."""
    rnaseq_dir = os.path.join(dataset_path, "RNAseqAnalysis")
    
    # Process each contrast
    if dataset_id in ri.deg_data:
        for contrast_id, deg_df in ri.deg_data[dataset_id].items():
            contrast_dir = f"{package_name}/DEG_analysis/contrast_{contrast_id}"
            
            # Always save DEG results as CSV
            deg_csv = deg_df.to_csv(index=False)
            zf.writestr(f"{contrast_dir}/DEG_results.csv", deg_csv)
            
            # Always add static PNG plots if they exist
            contrast_path = os.path.join(rnaseq_dir, contrast_id)
            static_plots = [
                ("volcano_plot.png", "volcano_plot.png"),
                ("ma_plot.png", "ma_plot.png"),
                ("heatmap_top50.png", "heatmap_top50.png")
            ]
            
            for source_file, dest_file in static_plots:
                source_path = os.path.join(contrast_path, source_file)
                if os.path.exists(source_path):
                    with open(source_path, 'rb') as f:
                        zf.writestr(f"{contrast_dir}/{dest_file}", f.read())
            
            # Only generate interactive plot PDFs if requested
            if include_pdfs:
                try:
                    from single_analysis_plots import create_volcano_plot, create_ma_plot, create_deg_heatmap, load_sample_groups
                    
                    # Volcano plot PDF
                    volcano_fig = create_volcano_plot(deg_df)
                    if volcano_fig:
                        pdf_bytes = plotly_fig_to_pdf_bytes(volcano_fig, width=1200, height=800, scale=2)
                        if pdf_bytes:
                            zf.writestr(f"{contrast_dir}/volcano_plot.pdf", pdf_bytes)
                    
                    # MA plot PDF
                    ma_fig = create_ma_plot(deg_df)
                    if ma_fig:
                        pdf_bytes = plotly_fig_to_pdf_bytes(ma_fig, width=1200, height=800, scale=2)
                        if pdf_bytes:
                            zf.writestr(f"{contrast_dir}/ma_plot.pdf", pdf_bytes)
                    
                    # Heatmap PDF (if CPM data available)
                    if dataset_id in ri.cpm_data:
                        cpm_df = ri.cpm_data[dataset_id]
                        analysis_info = ri.analysis_info.get(dataset_id)
                        groups = load_sample_groups(ri.results_dir, dataset_id, cpm_df, analysis_info)
                        
                        heatmap_fig = create_deg_heatmap(cpm_df, deg_df, groups)
                        if heatmap_fig:
                            pdf_bytes = plotly_fig_to_pdf_bytes(heatmap_fig, width=1200, height=1000, scale=2)
                            if pdf_bytes:
                                zf.writestr(f"{contrast_dir}/heatmap.pdf", pdf_bytes)
                                
                except Exception as e:
                    logger.debug(f"Could not generate interactive plot PDFs for {contrast_id}: {e}")


def _generate_dataset_readme(
    ri: 'ResultsIntegrator',
    dataset_id: str,
    dataset_accession: str,
    components: Optional[Dict[str, bool]] = None
) -> str:
    """Generate README content for the dataset package."""
    if components is None:
        components = {'metadata': True, 'abundance': True, 'qc_plots': True, 'deg_analysis': True}
    
    readme_lines = [
        f"UORCA Analysis Package for {dataset_accession}",
        "=" * 50,
        "",
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        "PACKAGE CONTENTS",
        "-" * 20,
        ""
    ]
    
    section_num = 1
    
    # Add sections based on what was included
    if components.get('metadata', False):
        readme_lines.extend([
            f"{section_num}. metadata/",
            f"   - {dataset_accession}_metadata.csv: Original dataset metadata with sample information",
            "   - contrasts.csv: Contrast definitions and descriptions",
            ""
        ])
        section_num += 1
    
    if components.get('abundance', False):
        readme_lines.extend([
            f"{section_num}. abundance/",
            "   - Sample-specific Kallisto quantification results",
            "   - Each sample directory contains:",
            "     * abundance.tsv: Transcript-level abundance estimates",
            "     * run_info.json: Kallisto run information",
            ""
        ])
        section_num += 1
    
    if components.get('qc_plots', False):
        readme_lines.extend([
            f"{section_num}. QC_plots/",
            "   - PCA_plot.pdf: Principal component analysis (if PDF option selected)",
            "   - MDS.png: Multi-dimensional scaling plot",
            "   - filtering_density.png: Gene expression filtering results",
            "   - normalization_boxplots.png: Before/after normalization",
            "   - voom_mean_variance.png: Mean-variance relationship",
            "   - sa_plot.png: Sigma vs average expression",
            ""
        ])
        section_num += 1
    
    if components.get('deg_analysis', False):
        readme_lines.extend([
            f"{section_num}. DEG_analysis/",
            "   - CPM_data.csv: Normalized expression values (counts per million)",
            "   - contrast_*/: Results for each contrast containing:",
            "     * DEG_results.csv: Differential expression statistics",
            "     * volcano_plot.pdf/png: Volcano plot visualization",
            "     * ma_plot.pdf/png: MA plot visualization",
            "     * heatmap.pdf/png: Top 50 DEGs heatmap",
            ""
        ])
    
    readme_lines.extend([
        "DATASET INFORMATION",
        "-" * 20,
        ""
    ])
    
    # Add dataset-specific information
    if dataset_id in ri.analysis_info:
        info = ri.analysis_info[dataset_id]
        readme_lines.extend([
            f"Organism: {info.get('organism', 'Unknown')}",
            f"Sample groups: {', '.join(info.get('unique_groups', []))}",
            f"Number of contrasts: {len(ri.deg_data.get(dataset_id, {}))}",
            ""
        ])
    
    # Add dataset title if available
    if hasattr(ri, 'dataset_info') and dataset_id in ri.dataset_info:
        title = ri.dataset_info[dataset_id].get('title', '')
        if title:
            if title.startswith('Title:'):
                title = title[6:].strip()
            readme_lines.extend([
                f"Title: {title}",
                ""
            ])
    
    readme_lines.extend([
        "NOTES",
        "-" * 20,
        "- All plots are static images (PNG or PDF format)",
        "- PDF files are high-resolution vector graphics suitable for publication",
        "- PNG files are raster images from the original analysis",
        "- DEG_results.csv contains full differential expression statistics",
        "- CPM_data.csv contains normalized expression values for all genes",
        "",
        "For questions or issues, please refer to the UORCA documentation."
    ])
    
    return "\n".join(readme_lines)


def setup_fragment_decorator():
    """Set up the fragment decorator."""
    from streamlit import fragment
    st.fragment = fragment


def _validate_results_dir(path: str) -> Tuple[bool, str]:
    """Simple validation for the user supplied results directory."""
    if not os.path.isdir(path):
        return False, "Directory does not exist"
    for root, dirs, files in os.walk(path):
        if "RNAseqAnalysis" in dirs or "DEG.csv" in files:
            return True, ""
    return False, "No RNAseqAnalysis data found"


@st.cache_resource
@log_streamlit_function
def get_integrator(path: str) -> Tuple[Optional[ResultsIntegrator], Optional[str]]:
    """Load and cache the ResultsIntegrator for the given path."""
    try:
        ri = ResultsIntegrator(results_dir=path)
        ri.load_data()
        return ri, None
    except Exception as e:
        return None, str(e)


@st.cache_data(show_spinner=False, ttl=3600)  # Cache for 1 hour
@log_streamlit_function
def cached_identify_important_genes(
    results_dir: str,
    selected_contrasts: List[Tuple[str, str]],
    top_frequent: int,
    p_value_threshold: float,
    lfc_threshold: float
) -> List[str]:
    """
    Cache the identify_important_genes computation for selected contrasts only.

    This function identifies genes that are commonly differentially expressed
    across the selected contrasts, not all contrasts in the dataset.

    Args:
        results_dir: Path to results directory
        selected_contrasts: List of (analysis_id, contrast_id) tuples to analyse
        top_frequent: Number of top frequently DE genes to return
        p_value_threshold: P-value threshold for significance
        lfc_threshold: Log fold change threshold for significance

    Returns:
        List of gene names sorted alphabetically
    """
    ri, error = get_integrator(results_dir)
    if error or not ri:
        return []

    try:
        # Sort contrasts to ensure deterministic caching regardless of order
        sorted_contrasts = sorted(selected_contrasts)

        # Count gene occurrences across selected contrasts only
        gene_counts = {}
        contrasts_processed = 0
        total_significant_genes = 0

        logger.info(f"Processing {len(sorted_contrasts)} selected contrasts for gene identification")
        logger.info(f"Using thresholds: P-value < {p_value_threshold}, |logFC| > {lfc_threshold}")

        for analysis_id, contrast_id in sorted_contrasts:
            if analysis_id in ri.deg_data and contrast_id in ri.deg_data[analysis_id]:
                deg_df = ri.deg_data[analysis_id][contrast_id]

                # Filter for significant genes
                if 'adj.P.Val' in deg_df.columns and 'logFC' in deg_df.columns and 'Gene' in deg_df.columns:
                    significant_genes = deg_df[
                        (deg_df['adj.P.Val'] < p_value_threshold) &
                        (abs(deg_df['logFC']) > lfc_threshold)
                    ]['Gene'].tolist()

                    contrasts_processed += 1
                    total_significant_genes += len(significant_genes)
                    logger.info(f"Contrast {analysis_id}:{contrast_id} - found {len(significant_genes)} significant genes")

                    # Count occurrences
                    for gene in significant_genes:
                        gene_counts[gene] = gene_counts.get(gene, 0) + 1
                else:
                    logger.warning(f"Missing required columns in {analysis_id}:{contrast_id}")
            else:
                logger.warning(f"Contrast {analysis_id}:{contrast_id} not found in DEG data")

        # Sort by frequency (descending) then by gene name (ascending)
        sorted_genes = sorted(gene_counts.items(), key=lambda x: (-x[1], x[0]))

        # Log summary statistics
        unique_genes_found = len(gene_counts)
        logger.info(f"Gene identification summary:")
        logger.info(f"  - Contrasts processed: {contrasts_processed}/{len(sorted_contrasts)}")
        logger.info(f"  - Total significant gene occurrences: {total_significant_genes}")
        logger.info(f"  - Unique genes found: {unique_genes_found}")
        logger.info(f"  - Requested top genes: {top_frequent}")

        # Return top N most frequent genes
        top_genes = [gene for gene, count in sorted_genes[:top_frequent]]
        actual_returned = len(top_genes)

        logger.info(f"  - Genes returned: {actual_returned}")
        if actual_returned > 0:
            # Show frequency distribution of top genes
            top_5_genes = sorted_genes[:min(5, len(sorted_genes))]
            logger.info(f"  - Top gene frequencies: {[(gene, count) for gene, count in top_5_genes]}")

        return top_genes

    except Exception as e:
        logger.error(f"Error identifying important genes: {e}")
        return []


# Add a cache for figure objects to improve performance
@st.cache_data(show_spinner=False, ttl=3600, hash_funcs={go.Figure: lambda _: None, ResultsIntegrator: lambda _: None})
@log_streamlit_function
def cached_figure_creation(
    func_name: str,
    results_dir: str,  # Use results_dir instead of RI instance for better caching
    *args,
    **kwargs) -> Optional[go.Figure]:
    """Cache figure objects to avoid recreating them. Cache expires after 1 hour."""
    ri, error = get_integrator(results_dir)
    if error or not ri:
        return None

    try:
        if func_name == "create_lfc_heatmap":
            try:
                return ri.create_lfc_heatmap(*args, **kwargs)
            except TypeError as e:
                # Backward compatibility: older RI may not accept 'cluster_genes'
                if "cluster_genes" in str(e):
                    safe_kwargs = dict(kwargs)
                    safe_kwargs.pop("cluster_genes", None)
                    return ri.create_lfc_heatmap(*args, **safe_kwargs)
                raise
        elif func_name == "create_expression_plots":
            # Map positional arguments to correct parameters for create_expression_plots
            if len(args) >= 15:
                genes, plot_type, analyses, output_file, hide_x_labels, page_number, facet_font_size, lock_y_axis, show_raw_points, legend_position, show_grid_lines, grid_opacity, selected_groups, plots_per_row, show_box = args[:15]
                return ri.create_expression_plots(
                    genes=genes,
                    plot_type=plot_type,
                    analyses=analyses,
                    output_file=output_file,
                    hide_x_labels=hide_x_labels,
                    page_number=page_number,
                    facet_font_size=facet_font_size,
                    lock_y_axis=lock_y_axis,
                    show_raw_points=show_raw_points,
                    legend_position=legend_position,
                    show_grid_lines=show_grid_lines,
                    grid_opacity=grid_opacity,
                    selected_groups=selected_groups,
                    plots_per_row=plots_per_row,
                    show_box=show_box
                )
            elif len(args) >= 14:
                genes, plot_type, analyses, output_file, hide_x_labels, page_number, facet_font_size, lock_y_axis, show_raw_points, legend_position, show_grid_lines, grid_opacity, selected_groups, plots_per_row = args[:14]
                return ri.create_expression_plots(
                    genes=genes,
                    plot_type=plot_type,
                    analyses=analyses,
                    output_file=output_file,
                    hide_x_labels=hide_x_labels,
                    page_number=page_number,
                    facet_font_size=facet_font_size,
                    lock_y_axis=lock_y_axis,
                    show_raw_points=show_raw_points,
                    legend_position=legend_position,
                    show_grid_lines=show_grid_lines,
                    grid_opacity=grid_opacity,
                    selected_groups=selected_groups,
                    plots_per_row=plots_per_row
                )
            elif len(args) >= 12:
                genes, plot_type, analyses, output_file, hide_x_labels, page_number, facet_font_size, lock_y_axis, show_raw_points, legend_position, show_grid_lines, grid_opacity = args[:12]
                return ri.create_expression_plots(
                    genes=genes,
                    plot_type=plot_type,
                    analyses=analyses,
                    output_file=output_file,
                    hide_x_labels=hide_x_labels,
                    page_number=page_number,
                    facet_font_size=facet_font_size,
                    lock_y_axis=lock_y_axis,
                    show_raw_points=show_raw_points,
                    legend_position=legend_position,
                    show_grid_lines=show_grid_lines,
                    grid_opacity=grid_opacity
                )
            else:
                return ri.create_expression_plots(*args, **kwargs)
    except Exception as e:
        logger.error(f"Error creating cached figure: {e}")
    return None


@st.cache_data(show_spinner=False, ttl=3600)
@log_streamlit_function
def cached_get_all_genes_from_integrator(results_dir: str) -> List[str]:
    """Extract all unique genes from all datasets in the integrator with caching."""
    ri, error = get_integrator(results_dir)
    if error or not ri:
        return []

    try:
        all_genes = set()
        for cpm_df in ri.cpm_data.values():
            if 'Gene' in cpm_df.columns:
                all_genes.update(cpm_df['Gene'].tolist())
        return sorted(all_genes)
    except Exception as e:
        logger.error(f"Error getting all genes: {e}")
        return []

@log_streamlit_function
def get_all_genes_from_integrator(ri: ResultsIntegrator) -> List[str]:
    """Extract all unique genes from all datasets in the integrator."""
    return cached_get_all_genes_from_integrator(ri.results_dir)





@log_streamlit_function
def calculate_pagination_info(gene_sel: List[str], genes_per_page: int = 30) -> Tuple[int, int, int, List[str]]:
    """
    Calculate pagination information for gene lists.

    Returns:
        total_pages, current_page, genes_per_page, current_genes
    """
    if not gene_sel:
        return 1, 1, genes_per_page, []

    total_pages = (len(gene_sel) + genes_per_page - 1) // genes_per_page
    current_page = st.session_state.get('page_num', 1)

    # Calculate gene slice for the current page
    start_idx = (current_page - 1) * genes_per_page
    end_idx = min(start_idx + genes_per_page, len(gene_sel))
    current_genes = gene_sel[start_idx:end_idx]

    return total_pages, current_page, genes_per_page, current_genes


def safe_rerun():
    """Call streamlit rerun."""
    st.rerun()


def check_ai_generating():
    """Check if AI is currently generating to skip expensive operations."""
    return st.session_state.get('ai_generating', False)


@log_streamlit_function
def add_custom_css():
    """Add custom CSS styles for the Streamlit app."""
    st.markdown("""
    <style>
      /* Fix for text wrapping in dataframes */
      .stDataFrame tbody tr td {
        white-space: normal !important;
        word-wrap: break-word !important;
        max-width: 300px;
      }
      .stDataFrame th {
        white-space: normal !important;
        word-wrap: break-word !important;
        max-width: 300px;
      }
      /* Ensure Description column has more space */
      .stDataFrame td:nth-child(3) {
        min-width: 250px;
        max-width: 500px;
      }
      /* Compact multiselect tags */
      span[data-baseweb="tag"] {
        font-size: 11px !important;
        padding: 0.25rem 0.5rem !important;
        height: 1.2rem !important;
      }
    </style>
    """, unsafe_allow_html=True)


@log_streamlit_function
def load_environment():
    """Load environment variables from .env file if available."""
    try:
        from dotenv import load_dotenv
        load_dotenv()
        # Also try loading from project root
        project_root = Path(__file__).parent.parent.parent.parent
        env_file = project_root / ".env"
        if env_file.exists():
            load_dotenv(env_file)
    except ImportError:
        # dotenv not available, continue without it
        pass


@log_streamlit_function
def get_organisms_from_datasets(ri: ResultsIntegrator, selected_datasets: List[str]) -> List[str]:
    """
    Get unique organisms from selected datasets.

    Args:
        ri: ResultsIntegrator instance
        selected_datasets: List of dataset IDs

    Returns:
        Sorted list of unique organism names
    """
    organisms = set()
    for dataset_id in selected_datasets:
        if dataset_id in ri.analysis_info:
            organism = ri.analysis_info[dataset_id].get('organism', 'Unknown')
            organisms.add(organism)
    return sorted(list(organisms))


@log_streamlit_function
def group_datasets_by_organism(ri: ResultsIntegrator, selected_datasets: List[str]) -> Dict[str, List[str]]:
    """
    Group datasets by their organism.

    Args:
        ri: ResultsIntegrator instance
        selected_datasets: List of dataset IDs

    Returns:
        Dictionary mapping organism names to lists of dataset IDs
    """
    organism_groups = {}
    for dataset_id in selected_datasets:
        if dataset_id in ri.analysis_info:
            organism = ri.analysis_info[dataset_id].get('organism', 'Unknown')
            if organism not in organism_groups:
                organism_groups[organism] = []
            organism_groups[organism].append(dataset_id)
    return organism_groups


@log_streamlit_function
def group_contrasts_by_organism(ri: ResultsIntegrator, selected_contrasts: List[Tuple[str, str]]) -> Dict[str, List[Tuple[str, str]]]:
    """
    Group contrasts by their dataset's organism.

    Args:
        ri: ResultsIntegrator instance
        selected_contrasts: List of (analysis_id, contrast_id) tuples

    Returns:
        Dictionary mapping organism names to lists of contrast tuples
    """
    organism_groups = {}
    for analysis_id, contrast_id in selected_contrasts:
        if analysis_id in ri.analysis_info:
            organism = ri.analysis_info[analysis_id].get('organism', 'Unknown')
            if organism not in organism_groups:
                organism_groups[organism] = []
            organism_groups[organism].append((analysis_id, contrast_id))
    return organism_groups


@log_streamlit_function
def filter_genes_by_organism(ri: ResultsIntegrator, genes: List[str], organism: str, selected_contrasts: List[Tuple[str, str]]) -> List[str]:
    """
    Filter genes to only include those found in datasets of a specific organism (for contrasts).

    Args:
        ri: ResultsIntegrator instance
        genes: List of gene names to filter
        organism: Target organism name
        selected_contrasts: List of (analysis_id, contrast_id) tuples for this organism

    Returns:
        List of genes that are present in the organism's datasets
    """
    organism_genes = set()

    # Collect all genes from datasets of this organism
    for analysis_id, contrast_id in selected_contrasts:
        if analysis_id in ri.analysis_info and ri.analysis_info[analysis_id].get('organism') == organism:
            # Check DEG data
            if analysis_id in ri.deg_data and contrast_id in ri.deg_data[analysis_id]:
                deg_df = ri.deg_data[analysis_id][contrast_id]
                if 'Gene' in deg_df.columns:
                    organism_genes.update(deg_df['Gene'].tolist())

            # Check CPM data
            if analysis_id in ri.cpm_data:
                cpm_df = ri.cpm_data[analysis_id]
                if 'Gene' in cpm_df.columns:
                    organism_genes.update(cpm_df['Gene'].tolist())

    # Return only genes that are in the input list AND found in this organism's data
    return [gene for gene in genes if gene in organism_genes]


@log_streamlit_function
def filter_genes_by_organism_datasets(ri: ResultsIntegrator, genes: List[str], organism: str, selected_datasets: List[str]) -> List[str]:
    """
    Filter genes to only include those found in datasets of a specific organism (for datasets).

    Args:
        ri: ResultsIntegrator instance
        genes: List of gene names to filter
        organism: Target organism name
        selected_datasets: List of dataset IDs for this organism

    Returns:
        List of genes that are present in the organism's datasets
    """
    organism_genes = set()

    # Collect all genes from datasets of this organism
    for analysis_id in selected_datasets:
        if analysis_id in ri.analysis_info and ri.analysis_info[analysis_id].get('organism') == organism:
            # Check CPM data (primary source for expression plots)
            if analysis_id in ri.cpm_data:
                cpm_df = ri.cpm_data[analysis_id]
                if 'Gene' in cpm_df.columns:
                    organism_genes.update(cpm_df['Gene'].tolist())

    # Return only genes that are in the input list AND found in this organism's data
    return [gene for gene in genes if gene in organism_genes]


@log_streamlit_function
def get_organism_display_name(organism: str) -> str:
    """
    Get a user-friendly display name for an organism.

    Args:
        organism: Scientific or common name of organism

    Returns:
        Cleaned display name for UI
    """
    if not organism or organism == 'Unknown':
        return 'Unknown Species'

    # Handle common organism name patterns
    organism = organism.strip()

    # Capitalize first letter of each word for better display
    if len(organism.split()) <= 2:  # Scientific names typically have 1-2 words
        return organism.title()
    else:
        return organism


def is_analysis_successful(ri, analysis_id: str) -> bool:
    """
    Check if an analysis was successful by checking analysis_info.json.

    Args:
        ri: ResultsIntegrator instance
        analysis_id: Analysis ID to check

    Returns:
        True if analysis was successful, False otherwise
    """
    if analysis_id not in ri.analysis_info:
        return False

    info = ri.analysis_info[analysis_id]
    return info.get("analysis_success", False)


def get_valid_contrasts_from_analysis_info(ri, analysis_id: str) -> List[Dict[str, str]]:
    """
    Get valid contrasts from analysis_info.json for a successful analysis.

    Args:
        ri: ResultsIntegrator instance
        analysis_id: Analysis ID to get contrasts for

    Returns:
        List of contrast dictionaries with name, description, etc.
    """
    if not is_analysis_successful(ri, analysis_id):
        return []

    info = ri.analysis_info[analysis_id]
    return info.get("contrasts", [])


def validate_contrast_has_deg_data(ri, analysis_id: str, contrast_name: str) -> bool:
    """
    Validate that a contrast has associated DEG data.

    Args:
        ri: ResultsIntegrator instance
        analysis_id: Analysis ID
        contrast_name: Contrast name to validate

    Returns:
        True if contrast has DEG data, False otherwise
    """
    return (analysis_id in ri.deg_data and
            contrast_name in ri.deg_data[analysis_id] and
            not ri.deg_data[analysis_id][contrast_name].empty)


def get_valid_contrasts_with_data(ri, analysis_ids: List[str] = None) -> List[Dict[str, Any]]:
    """
    Get all valid contrasts that have successful analysis and DEG data.

    This is the standardized way to get contrasts across all tabs.

    Args:
        ri: ResultsIntegrator instance
        analysis_ids: Optional list of analysis IDs to filter by. If None, use all.

    Returns:
        List of contrast information dictionaries with consistent names
    """
    from typing import List, Dict, Any

    if analysis_ids is None:
        analysis_ids = list(ri.analysis_info.keys())

    valid_contrasts = []

    for analysis_id in analysis_ids:
        # Step 1: Check if analysis was successful
        if not is_analysis_successful(ri, analysis_id):
            continue

        # Step 2: Get contrasts from analysis_info
        contrasts = get_valid_contrasts_from_analysis_info(ri, analysis_id)

        # Step 3: Validate each contrast has DEG data and use consistent names
        for contrast in contrasts:
            original_contrast_name = contrast["name"]
            if validate_contrast_has_deg_data(ri, analysis_id, original_contrast_name):
                info = ri.analysis_info[analysis_id]
                accession = info.get("accession", analysis_id)

                # Find the consistent contrast name from contrast_info
                consistent_name = original_contrast_name
                for contrast_key, contrast_data in ri.contrast_info.items():
                    if (contrast_data.get('original_name') == original_contrast_name and
                        contrast_data.get('analysis_id') == analysis_id):
                        consistent_name = contrast_data.get('name', original_contrast_name)
                        break

                valid_contrasts.append({
                    "analysis_id": analysis_id,
                    "accession": accession,
                    "contrast_name": consistent_name,
                    "original_contrast_name": original_contrast_name,
                    "description": contrast.get("description", ""),
                    "expression": contrast.get("expression", ""),
                    "justification": contrast.get("justification", "")
                })

    return valid_contrasts


# Export all functions for easy importing
__all__ = [
    # Main helper functions
    'get_integrator',
    'cached_identify_important_genes',
    'add_custom_css',
    'setup_fragment_decorator',
    'get_all_genes_from_integrator',
    'cached_get_all_genes_from_integrator',
    'calculate_pagination_info',
    'safe_rerun',
    'check_ai_generating',
    'cached_figure_creation',
    'load_environment',

    # Contrast validation functions
    'is_analysis_successful',
    'get_valid_contrasts_from_analysis_info',
    'validate_contrast_has_deg_data',
    'get_valid_contrasts_with_data',

    # Organism/species helper functions
    'get_organisms_from_datasets',
    'group_datasets_by_organism',
    'group_contrasts_by_organism',
    'filter_genes_by_organism',
    'filter_genes_by_organism_datasets',
    'get_organism_display_name',
    
    # PDF export helpers
    'plotly_fig_to_pdf_bytes',
    'generate_plot_filename',
    
    # Dataset download helpers
    'create_dataset_download_package',
    'add_directory_to_zip',

    # Streamlit logging functions
    'setup_streamlit_logging',
    'log_streamlit_function',
    'log_streamlit_agent',
    'log_streamlit_tab',
    'log_streamlit_event',
    'log_streamlit_data_load',
    'log_streamlit_user_action',

    # AI agent tool logging functions
    'start_ai_analysis_session',
    'get_ai_tool_logs_for_display',
    'clear_ai_tool_logs',
    'get_ai_tool_logger',
    'get_current_log_file',
    'read_log_file_contents',

    # Private utilities
    '_validate_results_dir',
    'ModuleFilter'
]

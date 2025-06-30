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


def short_label(full_label: str) -> str:
    """Create short labels for contrast multiselect display."""
    # keeps 'KO_vs_WT' out of 'GSE12345:KO_vs_WT – long sentence …'
    return full_label.split(":", 1)[-1].split(" –")[0].split(" - ")[0][:25]


def setup_fragment_decorator():
    """Set up the fragment decorator with fallbacks for older Streamlit versions."""
    # Check if fragment is available (Streamlit >=1.33.0)
    # If not, fallback to experimental_fragment
    try:
        from streamlit import fragment
        st.fragment = fragment
    except ImportError:
        try:
            from streamlit import experimental_fragment
            st.fragment = experimental_fragment
        except ImportError:
            # Fallback for very old Streamlit versions
            def fragment(func):
                """Fallback decorator when st.fragment is not available."""
                def wrapper(*args, **kwargs):
                    return func(*args, **kwargs)
                return wrapper
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
    top_frequent: int,
    top_unique: int,
    max_contrasts_for_unique: int,
    min_unique_per_contrast: int,
    p_value_threshold: float,
    lfc_threshold: float
) -> List[str]:
    """
    Cache the identify_important_genes computation to avoid recomputation on every rerun.

    This function is cached based on the input parameters, so it will only recompute
    when the parameters change. Cache expires after 1 hour to handle data updates.
    """
    ri, error = get_integrator(results_dir)
    if error or not ri:
        return []

    try:
        top_genes = ri.identify_important_genes(
            top_frequent=top_frequent,
            top_unique=top_unique,
            max_contrasts_for_unique=max_contrasts_for_unique,
            min_unique_per_contrast=min_unique_per_contrast,
            p_value_threshold=p_value_threshold,
            lfc_threshold=lfc_threshold
        )
        return sorted(top_genes)
    except Exception as e:
        logger.error(f"Error identifying important genes: {e}")
        return []


# Add a cache for figure objects to improve performance
@st.cache_data(show_spinner=False, ttl=1800, hash_funcs={go.Figure: lambda _: None, ResultsIntegrator: lambda _: None})
@log_streamlit_function
def cached_figure_creation(
    func_name: str,
    results_dir: str,  # Use results_dir instead of RI instance for better caching
    *args,
    **kwargs) -> Optional[go.Figure]:
    """Cache figure objects to avoid recreating them. Cache expires after 30 minutes."""
    ri, error = get_integrator(results_dir)
    if error or not ri:
        return None

    try:
        if func_name == "create_lfc_heatmap":
            return ri.create_lfc_heatmap(*args, **kwargs)
        elif func_name == "create_expression_plots":
            # Map positional arguments to correct parameters for create_expression_plots
            if len(args) >= 12:
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
def initialize_session_state(ri: ResultsIntegrator):
    """Initialize session state variables for dataset and contrast selections."""
    # Initialize session state for selections if not exists
    if 'selected_datasets' not in st.session_state:
        # Default to first 5 datasets
        all_dataset_ids = list(ri.cpm_data.keys())
        st.session_state['selected_datasets'] = set(all_dataset_ids[:5])

    if 'selected_contrasts' not in st.session_state:
        # Default to all contrasts for selected datasets
        selected_contrasts = set()
        for analysis_id in st.session_state['selected_datasets']:
            for c in ri.analysis_info.get(analysis_id, {}).get("contrasts", []):
                selected_contrasts.add((analysis_id, c["name"]))
        st.session_state['selected_contrasts'] = selected_contrasts

    # Initialize page number
    if 'page_num' not in st.session_state:
        st.session_state.page_num = 1


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
    """Safely call streamlit rerun with fallbacks for older versions."""
    try:
        st.rerun()
    except AttributeError:
        try:
            st.experimental_rerun()
        except AttributeError:
            # For very old versions, do nothing
            pass


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
        # Try loading from current directory first, then parent directories
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


# Export all functions for easy importing
__all__ = [
    # Main helper functions
    'get_integrator',
    'cached_identify_important_genes',
    'initialize_session_state',
    'add_custom_css',
    'setup_fragment_decorator',
    'get_all_genes_from_integrator',
    'cached_get_all_genes_from_integrator',
    'calculate_pagination_info',
    'safe_rerun',
    'check_ai_generating',
    'cached_figure_creation',
    'short_label',
    'load_environment',

    # Organism/species helper functions
    'get_organisms_from_datasets',
    'group_datasets_by_organism',
    'group_contrasts_by_organism',
    'filter_genes_by_organism',
    'filter_genes_by_organism_datasets',
    'get_organism_display_name',

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

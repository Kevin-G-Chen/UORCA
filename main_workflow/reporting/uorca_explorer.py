#!/usr/bin/env python3
"""
UORCA Explorer - Modular Streamlit App for Interactive Exploration of UORCA RNA-seq Results

This is the refactored modular version of the UORCA Explorer that separates functionality
into individual tab modules for better maintainability and development.

Usage:
  streamlit run uorca_explorer_modular.py

Or with custom port:
  streamlit run uorca_explorer_modular.py --server.port 8501
"""

import os
import sys
import logging
import streamlit as st
from pathlib import Path
from typing import Dict, Any, List, Tuple

# Add the current directory to the path for imports
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_dir)

# Import helper functions and setup
from streamlit_tabs.helpers import (
    _validate_results_dir,
    get_integrator,
    initialize_session_state,
    add_custom_css,
    load_environment,
    setup_fragment_decorator,
    setup_streamlit_logging,
    log_streamlit_function,
    log_streamlit_event,
    log_streamlit_data_load
)

# Import tab modules
from streamlit_tabs.sidebar_controls import render_sidebar_controls
from streamlit_tabs.data_selection_tab import render_data_selection_tab
from streamlit_tabs.heatmap_tab import render_heatmap_tab
from streamlit_tabs.expression_plots_tab import render_expression_plots_tab
from streamlit_tabs.analysis_plots_tab import render_analysis_plots_tab
from streamlit_tabs.datasets_info_tab import render_datasets_info_tab
from streamlit_tabs.contrasts_info_tab import render_contrasts_info_tab
from streamlit_tabs.ai_assistant_tab import render_ai_assistant_tab

# Import the main integrator
from ResultsIntegration import ResultsIntegrator

# Set up logging
logger = logging.getLogger(__name__)

# Load environment variables
load_environment()

# Set up fragment decorator
setup_fragment_decorator()

# Set page configuration
st.set_page_config(
    page_title="UORCA Explorer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Add custom CSS
add_custom_css()


@log_streamlit_function
def main():
    """Main application function."""
    # Initialize Streamlit-specific logging
    try:
        log_file = setup_streamlit_logging()
        log_streamlit_event(f"Streamlit logging initialized: {log_file}")
    except Exception as e:
        st.error(f"Failed to initialize logging: {e}")

    # Get the default results directory
    default_dir = os.getenv("UORCA_DEFAULT_RESULTS_DIR")
    if not default_dir:
        default_dir = os.path.join(os.path.dirname(os.path.dirname(script_dir)), "UORCA_results")
        if not os.path.exists(default_dir):
            default_dir = os.path.dirname(os.path.dirname(script_dir))

    # Determine initial results directory
    if os.path.exists('/UORCA_results'):
        initial_results_dir = '/UORCA_results'
    else:
        initial_results_dir = default_dir

    # Load and validate data
    ri, results_dir, validation_error = load_and_validate_data(initial_results_dir)

    if validation_error:
        render_error_state(validation_error)
        return

    if not ri or not ri.cpm_data:
        render_no_data_state()
        return

    # Store integrator and results directory in session state for AI assistant access
    st.session_state['results_integrator'] = ri
    st.session_state['results_dir'] = results_dir

    # Initialize session state
    initialize_session_state(ri)

    # Render sidebar controls and get parameters
    sidebar_params = render_sidebar_controls(ri, results_dir)

    # Render main interface
    render_main_interface(ri, results_dir, sidebar_params)


@log_streamlit_function
def load_and_validate_data(initial_results_dir: str) -> Tuple[ResultsIntegrator, str, str]:
    """
    Load and validate the results data.

    Returns:
        (ResultsIntegrator, results_dir, error_message)
    """
    # Get results directory from session state or use initial
    if 'results_dir' not in st.session_state:
        st.session_state.results_dir = initial_results_dir

    results_dir = st.session_state.results_dir

    # Validate directory
    valid_dir, validation_error = _validate_results_dir(results_dir)

    with st.sidebar.status("Loading data...", expanded=True) as status:
        if not valid_dir:
            log_streamlit_event(f"Directory validation failed: {validation_error}")
            status.update(label=f"Error: {validation_error}", state="error")
            return None, results_dir, validation_error

        log_streamlit_event(f"Loading data from: {results_dir}")
        ri, error = get_integrator(results_dir)

        if error:
            error_msg = f"Error loading data: {error}"
            log_streamlit_event(f"Data loading failed: {error}")
            status.update(label=error_msg, state="error")
            return None, results_dir, error_msg
        elif not ri or not ri.cpm_data:
            error_msg = "No data found. Please check the directory path."
            log_streamlit_event("No CPM data found in results")
            status.update(label=error_msg, state="error")
            return None, results_dir, error_msg
        else:
            log_streamlit_data_load("CPM datasets", len(ri.cpm_data))
            log_streamlit_data_load("DEG datasets", len(ri.deg_data))
            status.update(label=f"✅ Loaded {len(ri.cpm_data)} datasets", state="complete")

            # Check if this is the first time loading this directory
            if 'previous_results_dir' not in st.session_state or st.session_state.previous_results_dir != results_dir:
                st.session_state.previous_results_dir = results_dir
                log_streamlit_event(f"New results directory loaded: {results_dir}")

    return ri, results_dir, None


@log_streamlit_function
def render_main_interface(ri: ResultsIntegrator, results_dir: str, sidebar_params: Dict[str, Any]):
    """Render the main tabbed interface."""
    # Create main tabs
    tab_sel, tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs([
        "☑️ Select Data & Contrasts",
        "🌡️ Explore DEG Heatmap",
        "📈 Plot Gene Expression",
        "🧑‍🔬 Analyze Experiments",
        "📋 View Dataset Info",
        "🔍 View Contrast Info",
        "🤖 AI Assistant"
    ])

    # Get current selections from session state
    selected_contrasts = list(st.session_state.get('selected_contrasts', set()))
    selected_datasets = list(st.session_state.get('selected_datasets', set()))
    gene_sel = sidebar_params['gene_sel']

    log_streamlit_event(f"Rendering interface: {len(selected_datasets)} datasets, {len(selected_contrasts)} contrasts, {len(gene_sel)} genes")

    # Tab 1: Data Selection
    with tab_sel:
        render_data_selection_tab(
            ri=ri,
            pvalue_thresh=sidebar_params['pvalue_thresh'],
            lfc_thresh=sidebar_params['lfc_thresh']
        )

    # Tab 2: Heatmap
    with tab1:
        render_heatmap_tab(
            ri=ri,
            gene_sel=gene_sel,
            selected_contrasts=selected_contrasts,
            effective_pvalue_thresh=sidebar_params['effective_pvalue_thresh'],
            effective_lfc_thresh=sidebar_params['effective_lfc_thresh'],
            use_dynamic_filtering=sidebar_params['use_dynamic_filtering'],
            hide_empty_rows_cols=sidebar_params['hide_empty_rows_cols']
        )

    # Tab 3: Expression Plots
    with tab2:
        render_expression_plots_tab(
            ri=ri,
            gene_sel=gene_sel,
            selected_datasets=selected_datasets,
            hide_x_labels=sidebar_params['hide_x_labels']
        )

    # Tab 4: Analysis Plots
    with tab3:
        render_analysis_plots_tab(
            ri=ri,
            results_dir=results_dir
        )

    # Tab 5: Dataset Info
    with tab4:
        render_datasets_info_tab(ri=ri)

    # Tab 6: Contrast Info
    with tab5:
        render_contrasts_info_tab(
            ri=ri,
            pvalue_thresh=sidebar_params['pvalue_thresh'],
            lfc_thresh=sidebar_params['lfc_thresh']
        )

    # Tab 7: AI Assistant
    with tab6:
        render_ai_assistant_tab(
            ri=ri,
            results_dir=results_dir
        )


@log_streamlit_function
def render_error_state(error_message: str):
    """Render error state when data loading fails."""
    st.error(f"Failed to load data: {error_message}")
    render_help_info()


@log_streamlit_function
def render_no_data_state():
    """Render state when no data is found."""
    st.info("No data loaded. Please check your results directory path.")
    render_help_info()


@log_streamlit_function
def render_help_info():
    """Render help information when no data is loaded."""
    # Check if we're running in container mode
    container_mode = os.path.exists('/workspace') and os.path.exists('/UORCA_results')
    log_streamlit_event(f"Container mode detected: {container_mode}")

    if container_mode:
        st.info(
            """
            ## Welcome to UORCA Explorer (Modular)

            This app allows you to interactively explore RNA-seq results from UORCA analyses.

            ### Getting Started (Container Mode)
            1. The results directory should already be set to `/UORCA_results` in the sidebar
            2. If not, enter `/UORCA_results` as your results directory path
            3. The app will load your data and display interactive visualizations
            4. Use the sidebar controls to filter genes, datasets, and contrasts

            ### Container Path Mapping
            - **Inside container:** `/UORCA_results`
            - **Host system:** Your actual results directory path
            - Always use `/UORCA_results` when entering paths in this app

            ### Troubleshooting
            - Make sure the path contains valid UORCA analysis results
            - Each analysis should have a directory structure with:
              - RNAseqAnalysis/ directory with CPM.csv file
              - metadata/ directory with contrasts.csv and edger_analysis_samples.csv files

            ### New Modular Architecture
            This version of UORCA Explorer has been refactored into separate modules for:
            - Better maintainability and development
            - Improved caching of expensive operations (like gene identification)
            - Cleaner separation of concerns
            """
        )
    else:
        st.info(
            """
            ## Welcome to UORCA Explorer (Modular)

            This app allows you to interactively explore RNA-seq results from UORCA analyses.

            ### Getting Started
            1. Enter the path to your UORCA results directory in the sidebar
            2. The app will load your data and display interactive visualizations
            3. Use the sidebar controls to filter genes, datasets, and contrasts

            ### Troubleshooting
            - Make sure the path contains valid UORCA analysis results
            - Each analysis should have a directory structure with:
              - RNAseqAnalysis/ directory with CPM.csv file
              - metadata/ directory with contrasts.csv and edger_analysis_samples.csv files

            ### New Modular Architecture
            This version of UORCA Explorer has been refactored into separate modules for:
            - Better maintainability and development
            - Improved caching of expensive operations (like gene identification)
            - Cleaner separation of concerns

            ### Key Improvements
            - **Cached Gene Identification**: The expensive `identify_important_genes` operation is now cached and won't recompute on every app reload
            - **Modular Design**: Each tab is now in its own module for easier development
            - **Better Performance**: Fragment isolation prevents unnecessary recomputations
            - **Improved Structure**: Code is organized into logical modules with clear responsibilities
            """
        )


if __name__ == "__main__":
    main()

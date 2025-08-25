#!/usr/bin/env python3
"""
UORCA Explorer - Modular Streamlit App for Interactive Exploration of UORCA RNA-seq Results

This is the refactored modular version of the UORCA Explorer that separates functionality
into individual tab modules for better maintainability and development.

Usage:
  streamlit run uorca_explorer.py

Or with custom port:
  streamlit run uorca_explorer.py --server.port 8501
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
    add_custom_css,
    load_environment,
    setup_fragment_decorator,
    setup_streamlit_logging,
    log_streamlit_function,
    log_streamlit_event,
    log_streamlit_data_load,
    get_valid_contrasts_with_data
)

# Import tab modules
from streamlit_tabs.sidebar_controls import render_sidebar_controls

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
    layout="wide",
    initial_sidebar_state="expanded",
)

# Add custom CSS
add_custom_css()


@log_streamlit_function
def main():
    """Main application function."""
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

    # Render sidebar controls and get selected datasets
    sidebar_params = render_sidebar_controls(ri, results_dir)
    selected_datasets = sidebar_params.get('selected_datasets', [])

    # Render main interface
    render_main_interface(ri, results_dir, selected_datasets)


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

    if not valid_dir:
        log_streamlit_event(f"Directory validation failed: {validation_error}")
        st.sidebar.error(f"Error: {validation_error}")
        return None, results_dir, validation_error

    log_streamlit_event(f"Loading data from: {results_dir}")
    ri, error = get_integrator(results_dir)

    if error:
        error_msg = f"Error loading data: {error}"
        log_streamlit_event(f"Data loading failed: {error}")
        st.sidebar.error(error_msg)
        return None, results_dir, error_msg
    elif not ri or not ri.cpm_data:
        error_msg = "No data found. Please check the directory path."
        log_streamlit_event("No CPM data found in results")
        st.sidebar.error(error_msg)
        return None, results_dir, error_msg
    else:
        log_streamlit_data_load("CPM datasets", len(ri.cpm_data))

    # Use valid contrasts logic to get accurate counts, without logging full lists
    valid_contrasts = get_valid_contrasts_with_data(ri)
    n_datasets = len(ri.deg_data)
    n_valid_contrasts = len(valid_contrasts)

    # Simplified logging: show counts only
    log_streamlit_data_load(
        f"Valid DEG contrasts across {n_datasets} datasets",
        n_valid_contrasts
    )

    # Simple sidebar text instead of a status panel
    st.sidebar.caption(f"Loaded {len(ri.cpm_data)} datasets")

    # Check if this is the first time loading this directory
    if 'previous_results_dir' not in st.session_state or st.session_state.previous_results_dir != results_dir:
        st.session_state.previous_results_dir = results_dir
        log_streamlit_event(f"New results directory loaded: {results_dir}")

    return ri, results_dir, None


@log_streamlit_function
def render_main_interface(ri: ResultsIntegrator, results_dir: str, selected_datasets: List[str]):
    """Render the main tabbed interface."""

    # Create main tabs
    tab_ai, tab1, tab2, tab3, tab4, tab5 = st.tabs([
        "AI Assistant",
        "Explore DEG Heatmap",
        "Plot Gene Expression",
        "View Dataset Analyses",
        "View Dataset Info",
        "View Contrast Info"
    ])

    log_streamlit_event(f"Rendering interface with {len(selected_datasets)} selected datasets")

    # Tab 1: AI Assistant
    with tab_ai:
        render_ai_assistant_tab(
            ri=ri,
            results_dir=results_dir,
            selected_datasets=selected_datasets
        )

    # Tab 2: Heatmap (now handles its own contrast and gene selection)
    with tab1:
        render_heatmap_tab(
            ri=ri,
            selected_datasets=selected_datasets
        )

    # Tab 3: Expression Plots (now handles its own group and gene selection)
    with tab2:
        render_expression_plots_tab(
            ri=ri,
            selected_datasets=selected_datasets
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

    # Tab 6: Contrast Info (uses default thresholds since these are now tab-specific)
    with tab5:
        render_contrasts_info_tab(
            ri=ri,
            pvalue_thresh=0.05,  # Default threshold
            lfc_thresh=1.0       # Default threshold
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
            4. Use the sidebar to select datasets, then configure analysis parameters in each tab

            ### Container Path Mapping
            - **Inside container:** `/UORCA_results`
            - **Host system:** Your actual results directory path
            - Always use `/UORCA_results` when entering paths in this app

            ### New Modular Design
            **Dataset Selection:** Choose datasets in the sidebar
            **Heatmap Analysis:** Select contrasts and configure gene selection in the heatmap tab
            **Expression Analysis:** Choose sample groups and enter custom genes in the expression tab

            ### Troubleshooting
            - Make sure the path contains valid UORCA analysis results
            - Each analysis should have a directory structure with:
              - RNAseqAnalysis/ directory with CPM.csv file
              - metadata/ directory with contrasts.csv and edger_analysis_samples.csv files
            """
        )
    else:
        st.info(
            """
            ## Welcome to UORCA Explorer (Modular)

            This app allows you to interactively explore RNA-seq results from UORCA analyses.

            ### Getting Started
            1. Enter the path to your UORCA results directory in the sidebar
            2. Select datasets for analysis using the sidebar form
            3. Configure analysis parameters in each tab:
               - **Heatmap Tab:** Select contrasts and configure gene selection
               - **Expression Tab:** Choose sample groups and enter custom genes
            4. View interactive visualizations and results

            ### New Modular Architecture
            This version has been completely restructured for better usability:

            **Improved Workflow:**
            - **Step 1:** Dataset selection in sidebar (applies to all tabs)
            - **Step 2:** Tab-specific configuration (contrasts, genes, groups)
            - **Step 3:** Independent analysis in each tab

            **Key Improvements:**
            - **Independent Controls:** Each tab has its own parameter forms
            - **Flexible Gene Selection:** Heatmaps and expression plots use different gene selection methods
            - **Better Organization:** Clear separation between dataset selection and analysis configuration
            - **Enhanced Performance:** Cached operations and optimized rendering

            ### Troubleshooting
            - Make sure the path contains valid UORCA analysis results
            - Each analysis should have a directory structure with:
              - RNAseqAnalysis/ directory with CPM.csv file
              - metadata/ directory with contrasts.csv and edger_analysis_samples.csv files

            ### Tab-Specific Features
            - **Heatmap:** Contrast selection + gene selection (frequent DEGs or custom)
            - **Expression Plots:** Sample group selection + custom gene list only
            - **Analysis Plots:** Individual dataset QC and DE plots
            - **AI Assistant:** Intelligent analysis with automated contrast selection
            """
        )


if __name__ == "__main__":
    main()

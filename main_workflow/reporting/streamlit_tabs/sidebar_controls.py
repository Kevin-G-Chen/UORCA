"""
Sidebar Controls for UORCA Explorer - Dataset Selection Only.

This module handles dataset selection in the sidebar. Contrast and gene selection
have been moved to their respective tabs.
"""

import os
import logging
import streamlit as st
import pandas as pd
from typing import List, Dict, Any, Optional

from .helpers import (
    load_environment,
    log_streamlit_function,
    log_streamlit_event,
    log_streamlit_user_action
)
from ResultsIntegration import ResultsIntegrator

logger = logging.getLogger(__name__)

# Load environment variables
load_environment()


@log_streamlit_function
def render_sidebar_controls(ri: ResultsIntegrator, results_dir: str) -> Dict[str, Any]:
    """
    Render sidebar controls for dataset selection only.

    Args:
        ri: ResultsIntegrator instance
        results_dir: Path to results directory

    Returns:
        Dictionary containing selected datasets
    """
    st.sidebar.title("UORCA Explorer")

    # Help sections
    _render_tab_descriptions()
    _render_ai_assistant_help()
    _render_heatmap_help()
    _render_expression_plots_help()

    # Initialize return parameters with defaults
    params = {
        'selected_datasets': []
    }

    # Dataset Selection Form
    dataset_params = _render_dataset_selection_section(ri, results_dir)
    if dataset_params:
        params['selected_datasets'] = dataset_params.get('datasets', [])

    return params


@log_streamlit_function
def _render_dataset_selection_section(ri: ResultsIntegrator, results_dir: str) -> Optional[Dict[str, Any]]:
    """Render the dataset selection section."""
    with st.sidebar.expander("Dataset Selection", expanded=True):
        st.subheader("Select Datasets")
        st.markdown("Choose which datasets to include in your analysis.")
        st.markdown("Ensure you click **Apply Dataset Selection** if you make changes.")

        # Initialize session state for dataset selection
        if 'selected_datasets_from_sidebar' not in st.session_state:
            st.session_state['selected_datasets_from_sidebar'] = set()

        # Quick selection buttons (always available)
        col1, col2 = st.columns(2)
        with col1:
            if st.button("Select All", key="select_all_datasets_top"):
                all_datasets = set()
                dataset_data = _create_dataset_table_data(ri)
                if dataset_data:
                    all_datasets = set([d['Accession'] for d in dataset_data])
                st.session_state['selected_datasets_from_sidebar'] = all_datasets
                st.rerun()

        with col2:
            if st.button("Clear All", key="clear_all_datasets_top"):
                st.session_state['selected_datasets_from_sidebar'] = set()
                st.rerun()

        # Add search bar
        search_filter = st.text_input("Search datasets", "", key="sidebar_dataset_search")

        # Dataset Selection Form
        with st.form("dataset_selection_sidebar"):
            dataset_data = _create_dataset_table_data(ri)

            if dataset_data:
                df = pd.DataFrame(dataset_data)

                # Apply search filter
                if search_filter:
                    search_mask = df.apply(
                        lambda row: any(search_filter.lower() in str(val).lower() for val in row),
                        axis=1
                    )
                    df = df[search_mask]

                # Pre-select based on session state
                df['Select'] = df['Accession'].isin(st.session_state['selected_datasets_from_sidebar'])

                edited_df = st.data_editor(
                    df,
                    use_container_width=True,
                    hide_index=True,
                    column_config={
                        "Select": st.column_config.CheckboxColumn(
                            "",
                            help="Check to include this dataset",
                            default=False
                        ),
                        "Accession": st.column_config.TextColumn(
                            "Accession",
                            help="GEO accession number",
                            width=100
                        ),
                        "Organism": st.column_config.TextColumn(
                            "Organism",
                            help="Species/organism",
                            width=100
                        ),
                        "Samples": st.column_config.NumberColumn(
                            "Samples",
                            help="Number of samples",
                            width=85
                        ),
                        "Contrasts": st.column_config.NumberColumn(
                            "Contrasts",
                            help="Number of contrasts",
                            width=85
                        ),
                        "Title": st.column_config.TextColumn(
                            "Title",
                            help="Dataset title",
                            width=250
                        ),
                        "Description": st.column_config.TextColumn(
                            "Description",
                            help="Dataset description/summary",
                            width=300
                        )
                    },
                    key="dataset_selection_sidebar_table"
                )

                dataset_submitted = st.form_submit_button("Apply Dataset Selection", type="primary")

                if dataset_submitted:
                    # Get selected datasets
                    selected_datasets = set()
                    if not edited_df.empty:
                        selected_rows = edited_df[edited_df["Select"] == True]
                        selected_datasets = set(selected_rows["Accession"].tolist())

                    # Update session state
                    st.session_state['selected_datasets_from_sidebar'] = selected_datasets

                    log_streamlit_user_action(f"Dataset selection updated: {len(selected_datasets)} datasets selected")
                    st.rerun()
            else:
                st.info("No datasets available")
                st.form_submit_button("Apply Dataset Selection", disabled=True)

        # Show current selection summary (removed "No datasets selected" message)

        # Return results if we have selections
        if st.session_state['selected_datasets_from_sidebar']:
            return {
                'datasets': list(st.session_state['selected_datasets_from_sidebar'])
            }

    return None


@log_streamlit_function
def _create_dataset_table_data(ri: ResultsIntegrator) -> List[Dict[str, Any]]:
    """Create data for the dataset selection table for successful analyses only."""
    dataset_data = []

    # Only include datasets that have successful CPM data (indicator of successful analysis)
    for analysis_id in ri.cpm_data.keys():
        if analysis_id in ri.analysis_info:
            info = ri.analysis_info[analysis_id]

            # Get dataset title and description (same as dataset info tab)
            title = ""
            description = ""
            if hasattr(ri, "dataset_info") and analysis_id in getattr(ri, "dataset_info", {}):
                title = ri.dataset_info[analysis_id].get("title", "")
                if isinstance(title, str) and title.startswith("Title:"):
                    title = title[6:].strip()
                description = ri.dataset_info[analysis_id].get("summary", "")

            dataset_data.append({
                "Select": False,  # Default to unselected
                "Accession": info.get("accession", analysis_id),
                "Organism": info.get("organism", "Unknown"),
                "Samples": info.get("number_of_samples", 0),
                "Contrasts": info.get("number_of_contrasts", 0),
                "Title": title,
                "Description": description
            })

    return dataset_data


@log_streamlit_function
def _render_tab_descriptions():
    """Render expandable pane describing all tabs."""
    with st.sidebar.expander("Tab Guide", expanded=False):
        st.markdown("""
        **AI Assistant** - Use an AI agent to automatically analyse contrasts from your selected datasets. You will be able to see the steps the AI took to arrive at its findings.

        **Explore DEG Heatmap** - Create heatmaps showing log2 fold changes for selected genes across multiple contrasts.

        **Plot Gene Expression** - Generate boxplots displaying gene expression distributions across sample groups.

        **View Dataset Analyses** - Explore quality control and differential expression plots for individual datasets. View PCA plots, volcano plots, MA plots, and DEG heatmaps.

        **View Dataset Info** - Browse dataset information including study titles, summaries, experimental designs, and organism information.

        **View Contrast Info** - Browse all available contrasts across all datasets to get an overview of the comparisons available for analysis.
        """)


@log_streamlit_function
def _render_ai_assistant_help():
    """Render expandable pane describing AI assistant usage."""
    with st.sidebar.expander("Using AI Assistant", expanded=False):
        st.markdown("""
        1. Select datasets in the sidebar (e.g., GSE123456, GSE789012)
        2. Go to **AI Assistant** tab
        3. Enter your research question (e.g., "What contrasts are most relevant to T cell activation?")
        4. Click "**Run Complete AI Analysis**." Do note this will typically take a few minutes.
        5. You will be able to evaluate the agent's results and methodology

        **Note**: AI will only analyse contrasts from your selected datasets.
        """)


def _render_heatmap_help():
    """Render expandable pane describing heatmap creation."""
    with st.sidebar.expander("Creating Heatmaps", expanded=False):
        st.markdown("""
        1. Select datasets in the sidebar (e.g., GSE123456, GSE789012)
        2. Click "**Apply Dataset Selection**" - this will cause all contrasts associated with the datasets to appear in the **Explore DEG Heatmap** tab
        3. Go to **Explore DEG Heatmap** tab
        4. Choose contrasts from your selected datasets
        5. Configure gene selection (Frequent DEGs or Custom)
        6. Click "**Generate Heatmap Analysis**"
        """)


@log_streamlit_function
def _render_expression_plots_help():
    """Render expandable pane describing expression plot creation."""
    with st.sidebar.expander("Creating Expression Plots", expanded=False):
        st.markdown("""
        1. Select datasets in the sidebar (e.g., GSE123456, GSE789012)
        2. Click "**Apply Dataset Selection**" - this will cause all sample groups associated with the datasets to appear in the **Plot Gene Expression** tab
        3. Go to **Plot Gene Expression** tab
        4. Choose sample groups from your selected datasets
        5. Enter genes of interest (e.g., TP53, EGFR, MYC)
        6. Click "**Generate Expression Plots**"
        """)

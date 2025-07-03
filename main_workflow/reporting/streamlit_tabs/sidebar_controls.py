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
        st.markdown("**Select datasets for analysis**")
        st.info("Choose which datasets to include in your analysis. Contrast and gene selection are available in each tab.")

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

        # Dataset Selection Form
        with st.form("dataset_selection_sidebar"):
            dataset_data = _create_dataset_table_data(ri)

            if dataset_data:
                df = pd.DataFrame(dataset_data)
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

        # Show current selection summary
        if st.session_state['selected_datasets_from_sidebar']:
            st.success(f"✅ {len(st.session_state['selected_datasets_from_sidebar'])} datasets selected")
        else:
            st.info("No datasets selected")

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
            dataset_data.append({
                "Select": False,  # Default to unselected
                "Accession": info.get("accession", analysis_id),
                "Organism": info.get("organism", "Unknown"),
                "Samples": info.get("number_of_samples", 0),
                "Contrasts": info.get("number_of_contrasts", 0)
            })

    return dataset_data

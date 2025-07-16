"""
Contrasts Info Tab for UORCA Explorer.

This tab allows users to browse and filter contrast details.
"""

import os
import pandas as pd
import streamlit as st
from typing import List, Dict, Any, Set, Tuple

from .helpers import (
    check_ai_generating,
    setup_fragment_decorator,
    log_streamlit_tab,
    log_streamlit_function,
    log_streamlit_event,
    get_valid_contrasts_with_data
)
from ResultsIntegration import ResultsIntegrator

# Set up fragment decorator
setup_fragment_decorator()


@log_streamlit_tab("Contrasts Info")
def render_contrasts_info_tab(ri: ResultsIntegrator, pvalue_thresh: float, lfc_thresh: float):
    """
    Render the contrasts info tab.

    Args:
        ri: ResultsIntegrator instance
        pvalue_thresh: P-value threshold for DEG counting
        lfc_thresh: Log fold change threshold for DEG counting
    """
    st.header("View Contrast Info")
    st.markdown("Browse contrast details and filter by dataset or significance.")

    # Render the main contrasts interface
    _render_contrasts_interface(ri, pvalue_thresh, lfc_thresh)


@st.fragment
@log_streamlit_function
def _render_contrasts_interface(ri: ResultsIntegrator, pvalue_thresh: float, lfc_thresh: float):
    """Render the main contrasts interface using fragment isolation."""
    # Skip execution if AI is currently generating
    if check_ai_generating():
        return

    # Create contrast information DataFrame
    contrast_info = _create_contrast_info_dataframe(ri, pvalue_thresh, lfc_thresh)

    if not contrast_info:
        log_streamlit_event("No contrast information available")
        st.info("No contrast information available.")
        return

    df = pd.DataFrame(contrast_info)
    # Format DEG count column as integer
    if 'DEGs' in df.columns:
        df['DEGs'] = df['DEGs'].astype(int)

    # Add filtering options
    filtered_df = _render_filtering_controls(df)

    # Display the filtered contrast information
    if not filtered_df.empty:
        log_streamlit_event(f"Displaying {len(filtered_df)} contrasts")
        _render_contrast_table(filtered_df, pvalue_thresh, lfc_thresh)
        _render_selection_controls(filtered_df)
    else:
        log_streamlit_event("No contrasts match current filters")
        st.info("No contrasts match the current filters.")


@log_streamlit_function
def _create_contrast_info_dataframe(ri: ResultsIntegrator, pvalue_thresh: float, lfc_thresh: float) -> List[Dict[str, Any]]:
    """Create a list of contrast information dictionaries using standardized validation logic."""
    # Use the centralized validation function
    valid_contrasts = get_valid_contrasts_with_data(ri)

    contrast_info = []

    for contrast in valid_contrasts:
        analysis_id = contrast["analysis_id"]
        accession = contrast["accession"]
        contrast_name = contrast["contrast_name"]

        # Get dataset title for formatting
        dataset_title = ""
        if hasattr(ri, 'dataset_info') and analysis_id in ri.dataset_info:
            title = ri.dataset_info[analysis_id].get('title', '')
            if title and title.startswith('Title:'):
                title = title[6:].strip()
            dataset_title = title

        # Format dataset display name
        dataset_display = accession
        if dataset_title:
            dataset_display = f"{accession} - {dataset_title}"

        # Count DEGs for this contrast
        deg_count = 0
        deg_df = ri.deg_data[analysis_id][contrast_name]

        # Use exact column names from DEG.csv file - prefer adjusted p-value
        p_value_col = None
        if 'adj.P.Val' in deg_df.columns:
            p_value_col = 'adj.P.Val'  # Adjusted p-value (preferred)
        elif 'P.Value' in deg_df.columns:
            p_value_col = 'P.Value'  # Unadjusted p-value (fallback)

        lfc_col = 'logFC' if 'logFC' in deg_df.columns else None

        if p_value_col and lfc_col:
            deg_count = ((deg_df[p_value_col] < pvalue_thresh) &
                       (abs(deg_df[lfc_col]) > lfc_thresh)).sum()

        contrast_info.append({
            "Contrast Name": contrast_name,
            "Dataset": dataset_display,
            "Description": contrast["description"],
            "Justification": contrast["justification"],
            "Expression": contrast["expression"],
            "DEGs": int(deg_count)
        })

    return contrast_info


@log_streamlit_function
def _render_filtering_controls(df: pd.DataFrame) -> pd.DataFrame:
    """Render filtering controls and return filtered DataFrame."""
    st.subheader("Filter Contrasts")
    col1, col2, col3 = st.columns(3)

    with col1:
        dataset_filter = st.multiselect(
            "Filter by Dataset",
            options=sorted(df["Dataset"].unique()),
            default=[],
            key="contrasts_info_dataset_filter"
        )

    with col2:
        min_degs = st.number_input("Minimum DEGs", min_value=0, value=0, key="contrasts_info_min_degs")

    with col3:
        search_filter = st.text_input("Search Contrasts", "", key="contrasts_info_search")

    # Apply filters
    filtered_df = df
    if dataset_filter:
        filtered_df = filtered_df[filtered_df["Dataset"].isin(dataset_filter)]

    if min_degs > 0:
        filtered_df = filtered_df[filtered_df["DEGs"] >= min_degs]

    if search_filter:
        search_mask = filtered_df.apply(
            lambda row: any(search_filter.lower() in str(val).lower() for val in row),
            axis=1
        )
        filtered_df = filtered_df[search_mask]

    return filtered_df


@log_streamlit_function
def _render_contrast_table(filtered_df: pd.DataFrame, pvalue_thresh: float, lfc_thresh: float):
    """Render the contrast table for viewing."""
    # Sort by DEG count by default
    filtered_df = filtered_df.sort_values("DEGs", ascending=False)

    # Add parameter info for DEGs column
    deg_help_text = f"P-value < {pvalue_thresh}, |logFC| > {lfc_thresh}"

    # Display contrast information
    st.dataframe(
        filtered_df,
        hide_index=True,
        use_container_width=True,
        column_config={
            "Contrast Name": st.column_config.TextColumn("Contrast Name", width="small"),
            "Dataset": st.column_config.TextColumn("Dataset", width="medium"),
            "Description": st.column_config.TextColumn("Description", width="large"),
            "Justification": st.column_config.TextColumn("Justification", width="medium"),
            "Expression": st.column_config.TextColumn("Expression", width="small"),
            "DEGs": st.column_config.NumberColumn(
                "DEGs",
                format="%d",
                help=deg_help_text,
                width="small"
            )
        }
    )

    # Add small text about DEG parameters
    st.caption(f"DEG counts calculated using: {deg_help_text}")


@log_streamlit_function
def _render_selection_controls(filtered_df: pd.DataFrame):
    """Render information about contrast selection."""
    pass  # Removed info message about sidebar selection

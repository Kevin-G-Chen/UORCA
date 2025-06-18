"""
Contrasts Info Tab for UORCA Explorer.

This tab allows users to browse and filter contrast details.
"""

import pandas as pd
import streamlit as st
from typing import List, Dict, Any, Set, Tuple

from .helpers import check_ai_generating, setup_fragment_decorator, safe_rerun
from ResultsIntegration import ResultsIntegrator

# Set up fragment decorator
setup_fragment_decorator()


def render_contrasts_info_tab(ri: ResultsIntegrator, pvalue_thresh: float, lfc_thresh: float):
    """
    Render the contrasts info tab.

    Args:
        ri: ResultsIntegrator instance
        pvalue_thresh: P-value threshold for DEG counting
        lfc_thresh: Log fold change threshold for DEG counting
    """
    st.header("🔍 View Contrast Info")
    st.markdown("**🔍 Browse and filter contrast details.** View contrast descriptions, DEG counts, and filter by dataset or significance. Use this to understand what each comparison represents.")

    # Render the main contrasts interface
    _render_contrasts_interface(ri, pvalue_thresh, lfc_thresh)


@st.fragment
def _render_contrasts_interface(ri: ResultsIntegrator, pvalue_thresh: float, lfc_thresh: float):
    """Render the main contrasts interface using fragment isolation."""
    # Skip execution if AI is currently generating
    if check_ai_generating():
        return

    # Create contrast information DataFrame
    contrast_info = _create_contrast_info_dataframe(ri, pvalue_thresh, lfc_thresh)

    if not contrast_info:
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
        _render_contrast_table(filtered_df)
        _render_selection_controls(filtered_df)
    else:
        st.info("No contrasts match the current filters.")


def _create_contrast_info_dataframe(ri: ResultsIntegrator, pvalue_thresh: float, lfc_thresh: float) -> List[Dict[str, Any]]:
    """Create a list of contrast information dictionaries."""
    contrast_info = []
    for aid, contrasts in ri.deg_data.items():
        for contrast_id in contrasts.keys():
            # Get original contrast name and description
            original_name = contrast_id
            description = "No description available"

            # Try to get information from contrast_info
            if hasattr(ri, "contrast_info") and contrast_id in ri.contrast_info:
                description = ri.contrast_info[contrast_id].get('description', "No description available")
                # Use original name if available
                if 'name' in ri.contrast_info[contrast_id]:
                    original_name = ri.contrast_info[contrast_id]['name']

            # Count DEGs for this contrast
            deg_count = 0
            if aid in ri.deg_data and contrast_id in ri.deg_data[aid]:
                df = ri.deg_data[aid][contrast_id]
                # Use exact column names from DEG.csv file - prefer adjusted p-value
                p_value_col = None
                if 'adj.P.Val' in df.columns:
                    p_value_col = 'adj.P.Val'  # Adjusted p-value (preferred)
                elif 'P.Value' in df.columns:
                    p_value_col = 'P.Value'  # Unadjusted p-value (fallback)

                lfc_col = 'logFC' if 'logFC' in df.columns else None

                if p_value_col and lfc_col:
                    deg_count = ((df[p_value_col] < pvalue_thresh) &
                                 (abs(df[lfc_col]) > lfc_thresh)).sum()

            contrast_info.append({
                "Dataset": aid,
                "Accession": ri.analysis_info.get(aid, {}).get("accession", "Unknown"),
                "Contrast": original_name,
                "Original ID": contrast_id,
                "Description": description,
                "DEGs": deg_count
            })

    return contrast_info


def _render_filtering_controls(df: pd.DataFrame) -> pd.DataFrame:
    """Render filtering controls and return filtered DataFrame."""
    st.subheader("Filter Contrasts")
    col1, col2, col3 = st.columns(3)

    with col1:
        dataset_column = "Dataset" if "Dataset" in df.columns else "Accession"
        dataset_filter = st.multiselect(
            "Filter by Dataset",
            options=sorted(df[dataset_column].unique()),
            default=[]
        )

    with col2:
        min_degs = st.number_input("Minimum DEGs", min_value=0, value=0)

    with col3:
        search_filter = st.text_input("Search Contrasts", "")

    # Apply filters
    filtered_df = df
    if dataset_filter:
        # Get the column to filter by
        dataset_column = "Dataset" if "Dataset" in filtered_df.columns else "Accession"
        filtered_df = filtered_df[filtered_df[dataset_column].isin(dataset_filter)]

    if min_degs > 0:
        filtered_df = filtered_df[filtered_df["DEGs"] >= min_degs]

    if search_filter:
        search_mask = filtered_df.apply(
            lambda row: any(search_filter.lower() in str(val).lower() for val in row),
            axis=1
        )
        filtered_df = filtered_df[search_mask]

    return filtered_df


def _render_contrast_table(filtered_df: pd.DataFrame):
    """Render the interactive contrast table with selection checkboxes."""
    # Sort by DEG count by default
    filtered_df = filtered_df.sort_values("DEGs", ascending=False)

    # Add checkbox column for selection
    display_df = filtered_df.copy()
    display_df["✔"] = display_df.apply(
        lambda row: (row['Dataset'], row.get('Original ID', row['Contrast'])) in st.session_state.get('selected_contrasts', set()),
        axis=1
    )

    # Display contrast information with dataframe for interactivity
    edited_df = st.data_editor(
        display_df,
        hide_index=True,
        use_container_width=True,
        column_config={
            "✔": st.column_config.CheckboxColumn("Select", default=False),
            "Dataset": st.column_config.TextColumn("Dataset", width="medium"),
            "Accession": st.column_config.TextColumn("Accession", width="medium"),
            "Contrast": st.column_config.TextColumn("Contrast", width="medium"),
            "Original ID": st.column_config.TextColumn("Original ID", width="medium"),
            "Description": st.column_config.TextColumn("Description", width="large"),
            "DEGs": st.column_config.NumberColumn("DEGs", format="%d")
        },
        key="contrast_info_editor"
    )

    # Update selections based on checkboxes
    if not edited_df.empty:
        selected_from_info = set()
        for _, row in edited_df.iterrows():
            if row["✔"]:
                dataset = row['Dataset']
                contrast = row.get('Original ID', row['Contrast'])
                selected_from_info.add((dataset, contrast))
        st.session_state['selected_contrasts'] = selected_from_info


def _render_selection_controls(display_df: pd.DataFrame):
    """Render controls for selecting all visible contrasts."""
    # Add quick selection button
    if st.button("Select all visible contrasts", key="select_all_visible_contrasts"):
        visible_contrasts = set()
        for _, row in display_df.iterrows():
            dataset = row['Dataset']
            contrast = row.get('Original ID', row['Contrast'])
            visible_contrasts.add((dataset, contrast))
        st.session_state['selected_contrasts'] = visible_contrasts
        # Reset page number when changing contrasts
        st.session_state.page_num = 1
        st.success(f"Selected {len(visible_contrasts)} contrasts for analysis!")
        st.info("Switch to the Heat-map tab to view updated visualizations.")

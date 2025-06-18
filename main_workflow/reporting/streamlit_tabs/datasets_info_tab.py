"""
Datasets Info Tab for UORCA Explorer.

This tab allows users to browse and filter dataset metadata.
"""

import pandas as pd
import streamlit as st
from typing import List, Dict, Any, Set

from .helpers import check_ai_generating, setup_fragment_decorator, safe_rerun
from ..ResultsIntegration import ResultsIntegrator

# Set up fragment decorator
setup_fragment_decorator()


def render_datasets_info_tab(ri: ResultsIntegrator):
    """
    Render the datasets info tab.

    Args:
        ri: ResultsIntegrator instance
    """
    st.header("📋 View Dataset Info")
    st.markdown("**📋 Browse and filter dataset metadata.** View study details, organism information, sample counts, and experimental descriptions. Use filters to find specific datasets of interest.")

    # Render the main datasets interface
    _render_datasets_interface(ri)


@st.fragment
def _render_datasets_interface(ri: ResultsIntegrator):
    """Render the main datasets interface using fragment isolation."""
    # Skip execution if AI is currently generating
    if check_ai_generating():
        return

    # Create dataset information DataFrame
    dataset_info = _create_dataset_info_dataframe(ri)

    if not dataset_info:
        st.info("No dataset information available.")
        return

    df = pd.DataFrame(dataset_info)

    # Add filtering options
    filtered_df = _render_filtering_controls(df)

    # Display the filtered dataset information
    if not filtered_df.empty:
        _render_dataset_table(filtered_df)
        _render_dataset_details(filtered_df)
        _render_selection_controls(filtered_df)
    else:
        st.info("No datasets match the current filters.")


def _create_dataset_info_dataframe(ri: ResultsIntegrator) -> List[Dict[str, Any]]:
    """Create a list of dataset information dictionaries."""
    dataset_info = []
    for analysis_id, info in ri.analysis_info.items():
        # Build a dataset info dictionary
        dataset_dict = {
            "Accession": info.get("accession", "Unknown"),
            "Organism": info.get("organism", "Unknown"),
            "Number of Samples": info.get("number_of_samples", 0),
            "Number of Contrasts": info.get("number_of_contrasts", 0),
            "Dataset ID": analysis_id  # Keep dataset ID but place it last
        }

        # Add dataset metadata from analysis_info.json if available
        if hasattr(ri, "dataset_info") and analysis_id in getattr(ri, "dataset_info", {}):
            # Remove any "Title:" prefix from the title field
            title = ri.dataset_info[analysis_id].get("title", "")
            if isinstance(title, str) and title.startswith("Title:"):
                title = title[6:].strip()
            dataset_dict["Title"] = title

            dataset_dict["Summary"] = ri.dataset_info[analysis_id].get("summary", "")
            dataset_dict["Design"] = ri.dataset_info[analysis_id].get("design", "")

        dataset_info.append(dataset_dict)

    return dataset_info


def _render_filtering_controls(df: pd.DataFrame) -> pd.DataFrame:
    """Render filtering controls and return filtered DataFrame."""
    st.subheader("Filter Datasets")
    col1, col2 = st.columns(2)

    with col1:
        organism_filter = st.multiselect(
            "Filter by Organism",
            options=sorted(df["Organism"].unique()),
            default=[]
        )

    with col2:
        search_filter = st.text_input("Search Datasets", "")

    # Apply filters
    filtered_df = df
    if organism_filter:
        filtered_df = filtered_df[filtered_df["Organism"].isin(organism_filter)]

    if search_filter:
        search_mask = filtered_df.apply(
            lambda row: any(search_filter.lower() in str(val).lower() for val in row),
            axis=1
        )
        filtered_df = filtered_df[search_mask]

    return filtered_df


def _render_dataset_table(filtered_df: pd.DataFrame):
    """Render the interactive dataset table with selection checkboxes."""
    # Add checkbox column for selection
    display_df = filtered_df.copy()
    display_df["✔"] = display_df["Dataset ID"].isin(st.session_state.get('selected_datasets', set()))

    # Display dataset information with dataframe for interactivity
    edited_df = st.data_editor(
        display_df,
        hide_index=True,
        use_container_width=True,
        column_config={
            "✔": st.column_config.CheckboxColumn("Select", default=False),
            "Dataset ID": st.column_config.TextColumn("Dataset ID", width="medium"),
            "Title": st.column_config.TextColumn("Title", width="large"),
            "Accession": st.column_config.TextColumn("Accession", width="medium"),
            "Organism": st.column_config.TextColumn("Organism", width="medium"),
            "Number of Samples": st.column_config.NumberColumn("Samples", format="%d"),
            "Number of Contrasts": st.column_config.NumberColumn("Contrasts", format="%d")
        },
        key="dataset_info_editor"
    )

    # Update selections based on checkboxes
    if not edited_df.empty:
        selected_from_info = set(edited_df.loc[edited_df["✔"], "Dataset ID"].tolist())
        st.session_state['selected_datasets'] = selected_from_info


def _render_dataset_details(filtered_df: pd.DataFrame):
    """Render detailed information for each dataset."""
    # Show dataset details if available (always show by default)
    if any(col in filtered_df.columns for col in ["Title", "Summary", "Design"]):
        for _, row in filtered_df.iterrows():
            dataset_id = row.get("Dataset ID", row.get("Accession", "Unknown"))
            accession = row.get("Accession", "")

            with st.expander(f"Details for {dataset_id} {f'({accession})' if accession else ''}", expanded=True):
                if "Title" in filtered_df.columns and pd.notna(row.get("Title")):
                    st.subheader("Title")
                    st.markdown(str(row["Title"]))

                if "Summary" in filtered_df.columns and pd.notna(row.get("Summary")):
                    st.subheader("Summary")
                    st.markdown(str(row["Summary"]))

                if "Design" in filtered_df.columns and pd.notna(row.get("Design")):
                    st.subheader("Overall Design")
                    st.markdown(str(row["Design"]))


def _render_selection_controls(display_df: pd.DataFrame):
    """Render controls for selecting all visible datasets."""
    # Add quick selection button
    if st.button("Select all visible datasets", key="select_all_visible_datasets"):
        visible_datasets = set(display_df["Dataset ID"].tolist())
        st.session_state['selected_datasets'] = visible_datasets
        # Reset page number when changing datasets
        st.session_state.page_num = 1
        st.success(f"Selected {len(visible_datasets)} datasets for analysis!")
        st.info("Switch to the Heat-map or Expression tab to view updated visualizations.")

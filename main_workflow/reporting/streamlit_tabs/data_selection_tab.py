"""
Data Selection Tab for UORCA Explorer.

This tab allows users to select datasets and contrasts for analysis.
"""

import pandas as pd
import streamlit as st
from typing import Dict, Any, Set, Tuple

from .helpers import check_ai_generating, setup_fragment_decorator
from ResultsIntegration import ResultsIntegrator

# Set up fragment decorator
setup_fragment_decorator()


def render_data_selection_tab(ri: ResultsIntegrator, pvalue_thresh: float, lfc_thresh: float):
    """
    Render the data selection tab.

    Args:
        ri: ResultsIntegrator instance
        pvalue_thresh: P-value threshold for DEG counting
        lfc_thresh: Log fold change threshold for DEG counting
    """
    st.header("☑️ Select Data & Contrasts")
    st.markdown("**Select datasets and contrasts for analysis.** Use the tables below to choose your data, then click 'Regenerate Plots' to update visualizations.")

    # Dataset selection table
    st.subheader("Choose Datasets")
    _render_dataset_selection(ri)

    # Contrast selection table
    st.subheader("Choose Contrasts")
    _render_contrast_selection(ri, pvalue_thresh, lfc_thresh)

    # Selection summary
    _render_selection_summary()


@st.fragment
def _render_dataset_selection(ri: ResultsIntegrator):
    """Render the dataset selection interface."""
    # Skip execution if AI is currently generating
    if check_ai_generating():
        return

    dataset_rows = []
    for analysis_id, info in ri.analysis_info.items():
        title = getattr(ri, "dataset_info", {}).get(analysis_id, {}).get("title", "")
        title_display = title[:100] + ("..." if len(title) > 100 else "") if title else "No title available"

        dataset_rows.append({
            "✔": analysis_id in st.session_state.get('selected_datasets', set()),
            "Accession": info.get("accession", "Unknown"),
            "Dataset ID": analysis_id,
            "Title": title_display,
            "Organism": info.get("organism", "Unknown"),
            "Samples": info.get("number_of_samples", 0),
            "Contrasts": info.get("number_of_contrasts", 0)
        })

    if dataset_rows:
        ds_df = pd.DataFrame(dataset_rows)
        edited_ds = st.data_editor(
            ds_df,
            hide_index=True,
            use_container_width=True,
            column_config={
                "✔": st.column_config.CheckboxColumn("Select", default=False),
                "Title": st.column_config.TextColumn("Title", width="large"),
                "Accession": st.column_config.TextColumn("Accession", width="medium"),
                "Dataset ID": st.column_config.TextColumn("Dataset ID", width="medium")
            },
            key="selections_dataset_editor"
        )

        # Update selected datasets based on edited data
        if not edited_ds.empty:
            selected_datasets = set(edited_ds.loc[edited_ds["✔"], "Dataset ID"].tolist())
            st.session_state['selected_datasets'] = selected_datasets
            st.caption(f"📊 {len(selected_datasets)} datasets selected")


@st.fragment
def _render_contrast_selection(ri: ResultsIntegrator, pvalue_thresh: float, lfc_thresh: float):
    """Render the contrast selection interface."""
    # Skip execution if AI is currently generating
    if check_ai_generating():
        return

    if st.session_state['selected_datasets']:
        contrast_rows = []
        for analysis_id, contrasts in ri.deg_data.items():
            if analysis_id in st.session_state['selected_datasets']:
                for contrast_id in contrasts.keys():
                    # Get description
                    description = ri._get_contrast_description(analysis_id, contrast_id)
                    if description.startswith("Contrast: "):
                        description = description[10:]  # Remove "Contrast: " prefix

                    # Count DEGs
                    deg_count = 0
                    if analysis_id in ri.deg_data and contrast_id in ri.deg_data[analysis_id]:
                        df = ri.deg_data[analysis_id][contrast_id]
                        if 'adj.P.Val' in df.columns and 'logFC' in df.columns:
                            deg_count = ((df['adj.P.Val'] < pvalue_thresh) & (abs(df['logFC']) > lfc_thresh)).sum()

                    contrast_rows.append({
                        "✔": (analysis_id, contrast_id) in st.session_state.get('selected_contrasts', set()),
                        "Dataset": analysis_id,
                        "Accession": ri.analysis_info.get(analysis_id, {}).get("accession", "Unknown"),
                        "Contrast": contrast_id,
                        "Description": description[:150] + ("..." if len(description) > 150 else ""),
                        "DEGs": deg_count
                    })

        if contrast_rows:
            ctr_df = pd.DataFrame(contrast_rows)
            edited_ctr = st.data_editor(
                ctr_df,
                hide_index=True,
                use_container_width=True,
                column_config={
                    "✔": st.column_config.CheckboxColumn("Select", default=False),
                    "Description": st.column_config.TextColumn("Description", width="large"),
                    "Dataset": st.column_config.TextColumn("Dataset", width="medium"),
                    "Contrast": st.column_config.TextColumn("Contrast", width="medium"),
                    "DEGs": st.column_config.NumberColumn("DEGs", format="%d")
                },
                key="selections_contrast_editor"
            )

            # Update selected contrasts based on edited data
            if not edited_ctr.empty:
                selected_contrasts = set()
                for _, row in edited_ctr.iterrows():
                    if row["✔"]:
                        selected_contrasts.add((row["Dataset"], row["Contrast"]))
                st.session_state['selected_contrasts'] = selected_contrasts
                st.caption(f"🔍 {len(selected_contrasts)} contrasts selected")
        else:
            st.info("No contrasts available for selected datasets.")
    else:
        st.info("Please select at least one dataset to see available contrasts.")


def _render_selection_summary():
    """Render the selection summary section."""
    # Selection summary
    st.markdown("---")
    selected_datasets_count = len(st.session_state.get('selected_datasets', set()))
    selected_contrasts_count = len(st.session_state.get('selected_contrasts', set()))

    # Get gene count from session state or sidebar (this will need to be passed from main app)
    genes_count = len(st.session_state.get('current_gene_selection', []))

    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Datasets", selected_datasets_count)
    with col2:
        st.metric("Contrasts", selected_contrasts_count)
    with col3:
        st.metric("Genes", genes_count)

    if selected_datasets_count > 0 and selected_contrasts_count > 0 and genes_count > 0:
        st.success("✅ Ready for analysis! Switch to Heat-map or Expression tabs to view results.")
    else:
        missing = []
        if selected_datasets_count == 0:
            missing.append("datasets")
        if selected_contrasts_count == 0:
            missing.append("contrasts")
        if genes_count == 0:
            missing.append("genes")
        st.info(f"ℹ️ Please select {', '.join(missing)} to enable plot generation.")

"""
Sidebar Controls for UORCA Explorer - Form-Based Version with Dataset Selection.

This module handles all sidebar controls using Streamlit forms with expandable sections
for dataset selection, heatmap and expression plots configuration.
"""

import os
import logging
import streamlit as st
import pandas as pd
from typing import List, Dict, Any, Tuple, Optional

from .helpers import (
    cached_identify_important_genes,
    load_environment,
    log_streamlit_function,
    log_streamlit_event,
    log_streamlit_user_action,
    group_contrasts_by_organism,
    filter_genes_by_organism,
    get_organism_display_name
)
from ResultsIntegration import ResultsIntegrator

logger = logging.getLogger(__name__)

# Load environment variables
load_environment()


@log_streamlit_function
def render_sidebar_controls(ri: ResultsIntegrator, results_dir: str) -> Dict[str, Any]:
    """
    Render all sidebar controls using forms and return the selected parameters.

    Args:
        ri: ResultsIntegrator instance
        results_dir: Path to results directory

    Returns:
        Dictionary containing all selected parameters and gene selections
    """
    st.sidebar.title("UORCA Explorer")



    # Initialize return parameters with defaults
    params = {
        'gene_sel': [],
        'selected_contrasts': [],
        'selected_datasets': [],
        'heatmap_params': {
            'lfc_thresh': 1.0,
            'pvalue_thresh': 0.05,
            'gene_count': 50,
            'gene_selection_method': 'Frequent DEGs',
            'custom_genes': [],
            'contrasts': [],
            'datasets': []
        },

    }

    # Dataset Selection Forms
    dataset_params = _render_dataset_selection_section(ri, results_dir)
    if dataset_params:
        params['selected_datasets'] = dataset_params.get('datasets', [])
        params['selected_contrasts'] = dataset_params.get('contrasts', [])
        params['heatmap_params']['datasets'] = dataset_params.get('datasets', [])
        params['heatmap_params']['contrasts'] = dataset_params.get('contrasts', [])

    # Heatmap Configuration Form
    heatmap_params = _render_heatmap_form(ri, results_dir)
    if heatmap_params:
        params['heatmap_params'].update(heatmap_params)
        # Use contrasts from dataset selection if available, otherwise from heatmap
        if not params['selected_contrasts']:
            params['selected_contrasts'] = heatmap_params.get('contrasts', [])



    # Auto-select genes if we have contrasts
    if params['selected_contrasts']:
        gene_sel = _auto_select_genes(ri, results_dir, params['heatmap_params'])
        params['gene_sel'] = gene_sel
        st.session_state['current_gene_selection'] = gene_sel





    # Help section
    _render_help_section()

    return params





@log_streamlit_function
def _render_dataset_selection_section(ri: ResultsIntegrator, results_dir: str) -> Optional[Dict[str, Any]]:
    """Render the dataset selection section with two coordinated forms."""
    with st.sidebar.expander("Dataset & Contrast Selection", expanded=True):
        st.markdown("**Select your datasets and contrasts for analysis**")
        st.info("Selecting datasets will automatically update the contrast table below.")

        # Initialize session state for coordination
        if 'selected_datasets_from_form' not in st.session_state:
            st.session_state['selected_datasets_from_form'] = set()
        if 'selected_contrasts_from_form' not in st.session_state:
            st.session_state['selected_contrasts_from_form'] = set()

        # Form 1: Dataset Selection
        st.subheader("1. Select Datasets")
        with st.form("dataset_selection"):
            dataset_data = _create_dataset_table_data(ri)

            if dataset_data:
                df = pd.DataFrame(dataset_data)
                # Pre-select based on session state
                df['Select'] = df['Accession'].isin(st.session_state['selected_datasets_from_form'])

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
                    key="dataset_selection_table"
                )

                dataset_submitted = st.form_submit_button("Apply", type="primary")

                if dataset_submitted:
                    # Get selected datasets
                    selected_datasets = set()
                    if not edited_df.empty:
                        selected_rows = edited_df[edited_df["Select"] == True]
                        selected_datasets = set(selected_rows["Accession"].tolist())

                    # Update session state
                    st.session_state['selected_datasets_from_form'] = selected_datasets

                    # Auto-update contrasts based on selected datasets
                    if selected_datasets:
                        auto_contrasts = set()
                        for dataset in selected_datasets:
                            for contrast in ri.analysis_info.get(dataset, {}).get("contrasts", []):
                                auto_contrasts.add((dataset, contrast["name"]))
                        st.session_state['selected_contrasts_from_form'] = auto_contrasts
                    else:
                        st.session_state['selected_contrasts_from_form'] = set()

                    log_streamlit_user_action(f"Dataset form submitted: {len(selected_datasets)} datasets selected")
                    st.rerun()
            else:
                st.info("No datasets available")
                st.form_submit_button("Update Datasets", disabled=True)

        # Form 2: Contrast Selection
        st.subheader("2. Select Contrasts")
        with st.form("contrast_selection"):
            # Filter contrasts based on selected datasets
            if st.session_state['selected_datasets_from_form']:
                contrast_data = _create_contrast_table_data_filtered(
                    ri, st.session_state['selected_datasets_from_form']
                )
            else:
                contrast_data = []

            if contrast_data:
                df = pd.DataFrame(contrast_data)
                # Pre-select based on session state
                df['Select'] = df.apply(
                    lambda row: (row['Accession'], row['Contrast']) in st.session_state['selected_contrasts_from_form'],
                    axis=1
                )

                edited_df = st.data_editor(
                    df,
                    use_container_width=True,
                    hide_index=True,
                    column_config={
                        "Select": st.column_config.CheckboxColumn(
                            "",
                            help="Check to include this contrast",
                            default=False
                        ),
                        "Accession": st.column_config.TextColumn(
                            "Accession",
                            help="GEO accession number"
                        ),
                        "Contrast": st.column_config.TextColumn(
                            "Contrast",
                            help="Contrast identifier"
                        )
                    },
                    key="contrast_selection_table"
                )

                contrast_submitted = st.form_submit_button("Apply", type="secondary")

                if contrast_submitted:
                    # Get selected contrasts
                    selected_contrasts = set()
                    if not edited_df.empty:
                        selected_rows = edited_df[edited_df["Select"] == True]
                        selected_contrasts = set([
                            (row["Accession"], row["Contrast"])
                            for _, row in selected_rows.iterrows()
                        ])

                    # Update session state (but don't change dataset selection)
                    st.session_state['selected_contrasts_from_form'] = selected_contrasts

                    log_streamlit_user_action(f"Contrast form submitted: {len(selected_contrasts)} contrasts selected")
                    st.rerun()
            else:
                if st.session_state['selected_datasets_from_form']:
                    st.info("No contrasts available for selected datasets")
                else:
                    st.info("Please select datasets first")
                st.form_submit_button("Update Contrasts", disabled=True)

        # Return results if we have selections
        if st.session_state['selected_datasets_from_form'] or st.session_state['selected_contrasts_from_form']:
            return {
                'datasets': list(st.session_state['selected_datasets_from_form']),
                'contrasts': list(st.session_state['selected_contrasts_from_form'])
            }

    return None


@log_streamlit_function
def _render_heatmap_form(ri: ResultsIntegrator, results_dir: str) -> Optional[Dict[str, Any]]:
    """Render the heatmap configuration form."""
    with st.sidebar.expander("Heatmap Parameters", expanded=False):
        with st.form("heatmap_config"):
            st.subheader("Significance Thresholds")

            # Parameter inputs
            col1, col2 = st.columns(2)
            with col1:
                lfc_thresh = st.text_input(
                    "Log2FC Threshold",
                    value="1.0",
                    help="Absolute log2 fold change threshold (applies to both Frequent DEGs and Custom gene filtering)",
                    key="heatmap_lfc_threshold"
                )
            with col2:
                pvalue_thresh = st.text_input(
                    "P-value Threshold",
                    value="0.05",
                    help="Adjusted p-value threshold (applies to both Frequent DEGs and Custom gene filtering)",
                    key="heatmap_pvalue_threshold"
                )

            # Gene selection method
            st.subheader("Gene Selection")
            gene_selection_method = st.radio(
                "Gene Selection Method",
                options=["Frequent DEGs", "Custom"],
                index=0,
                key="heatmap_gene_selection_method",
                help="Choose between automatically selecting frequent DEGs or providing your own custom gene list"
            )

            if gene_selection_method == "Frequent DEGs":
                gene_count_input = st.text_input(
                    "Number of genes to display",
                    value="50",
                    help="Maximum number of top frequently differentially expressed genes to include",
                    key="heatmap_gene_count"
                )
                custom_genes_input = ""
            else:  # Custom
                st.info("Number of genes threshold does not apply to custom gene selection")
                custom_genes_input = st.text_area(
                    "Custom Gene List",
                    value="",
                    height=150,
                    placeholder="Enter one gene per line, e.g.:\nTP53\nEGFR\nMYC\nBRCA1",
                    help="Enter gene symbols, one per line. These genes will be filtered using the p-value and log2FC thresholds above.",
                    key="heatmap_custom_genes"
                )
                gene_count_input = "50"  # Default fallback

                # Show preview of custom genes
                if custom_genes_input.strip():
                    custom_genes_list = [gene.strip() for gene in custom_genes_input.strip().split('\n') if gene.strip()]
                    if custom_genes_list:
                        st.write(f"**Preview:** {len(custom_genes_list)} genes entered")
                        preview_text = ", ".join(custom_genes_list[:10])
                        if len(custom_genes_list) > 10:
                            preview_text += f", ... (+{len(custom_genes_list) - 10} more)"
                        st.caption(preview_text)

            # Validate parameters
            validation_error = False
            try:
                lfc_val = float(lfc_thresh)
                pval_val = float(pvalue_thresh)
            except ValueError:
                st.error("Please enter valid numeric values for thresholds")
                lfc_val, pval_val = 1.0, 0.05
                validation_error = True

            # Validate gene count for Frequent DEGs
            if gene_selection_method == "Frequent DEGs":
                try:
                    gene_count = int(gene_count_input)
                    if gene_count <= 0:
                        raise ValueError("Gene count must be positive")
                except ValueError:
                    st.error("Please enter a valid positive number for gene count")
                    gene_count = 50
                    validation_error = True
            else:
                gene_count = 50  # Not used for custom genes

            # Validate custom genes
            custom_genes_list = []
            if gene_selection_method == "Custom":
                if not custom_genes_input.strip():
                    st.error("Please enter at least one gene for custom selection")
                    validation_error = True
                else:
                    custom_genes_list = [gene.strip() for gene in custom_genes_input.strip().split('\n') if gene.strip()]
                    if not custom_genes_list:
                        st.error("Please enter valid gene names")
                        validation_error = True

            # Submit button
            submitted = st.form_submit_button("Apply", type="primary")

            if submitted and not validation_error:
                if gene_selection_method == "Frequent DEGs":
                    log_streamlit_user_action(f"Heatmap parameters updated: LFC={lfc_val}, P={pval_val}, genes={gene_count}, method=Frequent DEGs")
                else:
                    log_streamlit_user_action(f"Heatmap parameters updated: LFC={lfc_val}, P={pval_val}, method=Custom, custom_genes={len(custom_genes_list)}")

                return {
                    'lfc_thresh': lfc_val,
                    'pvalue_thresh': pval_val,
                    'gene_count': gene_count,
                    'gene_selection_method': gene_selection_method,
                    'custom_genes': custom_genes_list
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


@log_streamlit_function
def _create_contrast_table_data_filtered(ri: ResultsIntegrator, selected_datasets: set) -> List[Dict[str, Any]]:
    """Create data for the contrast selection table, filtered by selected datasets."""
    contrast_data = []

    for analysis_id in selected_datasets:
        if analysis_id in ri.analysis_info:
            info = ri.analysis_info[analysis_id]
            for contrast in info.get("contrasts", []):
                contrast_id = contrast["name"]

                contrast_data.append({
                    "Select": False,  # Default to unselected
                    "Accession": ri.analysis_info[analysis_id].get("accession", analysis_id),
                    "Contrast": contrast_id
                })

    return contrast_data


@log_streamlit_function
def _auto_select_genes(
    ri: ResultsIntegrator,
    results_dir: str,
    heatmap_params: Dict[str, Any]
) -> List[str]:
    """Select genes based on heatmap parameters (frequent DEGs or custom), considering species separation."""
    if not heatmap_params.get('contrasts'):
        return []

    # Check gene selection method
    gene_selection_method = heatmap_params.get('gene_selection_method', 'Frequent DEGs')

    try:
        # Group contrasts by organism first
        organism_groups = group_contrasts_by_organism(ri, heatmap_params['contrasts'])

        if gene_selection_method == 'Custom':
            # Custom gene selection - use provided gene list
            custom_genes = heatmap_params.get('custom_genes', [])
            if not custom_genes:
                st.sidebar.warning("No custom genes provided")
                return []

            # Filter custom genes to only include those present in the selected datasets
            available_genes = set()
            for organism_contrasts in organism_groups.values():
                for analysis_id, contrast_id in organism_contrasts:
                    if analysis_id in ri.deg_data and contrast_id in ri.deg_data[analysis_id]:
                        deg_df = ri.deg_data[analysis_id][contrast_id]
                        if 'Gene' in deg_df.columns:
                            available_genes.update(deg_df['Gene'].tolist())
                    if analysis_id in ri.cpm_data:
                        cpm_df = ri.cpm_data[analysis_id]
                        if 'Gene' in cpm_df.columns:
                            available_genes.update(cpm_df['Gene'].tolist())

            # Filter custom genes to only include available ones
            filtered_custom_genes = [gene for gene in custom_genes if gene in available_genes]
            missing_genes = [gene for gene in custom_genes if gene not in available_genes]

            # Display summary
            with st.sidebar.expander("Custom Gene Selection Summary", expanded=False):
                st.write(f"**Custom genes:** {len(custom_genes)} provided")
                st.write(f"**Available genes:** {len(filtered_custom_genes)} found in data")
                if missing_genes:
                    st.write(f"**Missing genes:** {len(missing_genes)} not found")
                    if len(missing_genes) <= 10:
                        st.write(f"Missing: {', '.join(missing_genes)}")
                    else:
                        st.write(f"Missing: {', '.join(missing_genes[:10])}, ... (+{len(missing_genes)-10} more)")

            log_streamlit_event(f"Custom gene selection: {len(filtered_custom_genes)} of {len(custom_genes)} genes available")
            return filtered_custom_genes

        else:
            # Frequent DEGs selection - use existing logic
            if len(organism_groups) == 1:
                # Single organism - use selected contrasts only
                organism = list(organism_groups.keys())[0]
                organism_contrasts = organism_groups[organism]

                top_genes = cached_identify_important_genes(
                    results_dir=results_dir,
                    selected_contrasts=organism_contrasts,
                    top_frequent=heatmap_params.get('gene_count', 50),
                    p_value_threshold=heatmap_params['pvalue_thresh'],
                    lfc_threshold=heatmap_params['lfc_thresh']
                )

                organism_display = get_organism_display_name(organism)

                # Display which contrasts were used for gene selection
                with st.sidebar.expander("Contrasts Used for Gene Selection", expanded=False):
                    st.write(f"**{organism_display}** - {len(organism_contrasts)} contrasts:")
                    for analysis_id, contrast_id in organism_contrasts:
                        # Get dataset accession for display
                        accession = ri.analysis_info.get(analysis_id, {}).get('accession', analysis_id)
                        st.write(f"• {accession}: {contrast_id}")

                log_streamlit_event(f"Auto-selected {len(top_genes)} frequent DEGs for {organism}")
                return top_genes

            else:
                # Multiple organisms - select genes for each and combine
                all_selected_genes = []
                gene_count_per_organism = heatmap_params.get('gene_count', 50)
                organism_summaries = []

                for organism, organism_contrasts in organism_groups.items():
                    # Get genes for this organism's contrasts only
                    organism_genes = cached_identify_important_genes(
                        results_dir=results_dir,
                        selected_contrasts=organism_contrasts,
                        top_frequent=gene_count_per_organism,
                        p_value_threshold=heatmap_params['pvalue_thresh'],
                        lfc_threshold=heatmap_params['lfc_thresh']
                    )

                    all_selected_genes.extend(organism_genes)
                    organism_display = get_organism_display_name(organism)
                    organism_summaries.append(f"{len(organism_genes)} genes from {organism_display}")

                # Display which contrasts were used for gene selection
                with st.sidebar.expander("Contrasts Used for Gene Selection", expanded=False):
                    for organism, organism_contrasts in organism_groups.items():
                        organism_display = get_organism_display_name(organism)
                        st.write(f"**{organism_display}** - {len(organism_contrasts)} contrasts:")
                        for analysis_id, contrast_id in organism_contrasts:
                            # Get dataset accession for display
                            accession = ri.analysis_info.get(analysis_id, {}).get('accession', analysis_id)
                            st.write(f"• {accession}: {contrast_id}")

                log_streamlit_event(f"Auto-selected {len(all_selected_genes)} frequent DEGs across {len(organism_groups)} organisms")
                return all_selected_genes

    except Exception as e:
        logger.error(f"Error in gene selection: {e}")
        st.sidebar.error(f"Error selecting genes: {str(e)}")
        return []





@log_streamlit_function
def _render_help_section():
    """Render the help section at the bottom of the sidebar."""
    st.sidebar.divider()
    with st.sidebar.expander("Help", expanded=False):
        st.markdown(
            """
            ### How to Use This Interface

            1. **Select Datasets**: Choose datasets in the first form and click "Apply"
            2. **Select Contrasts**: Fine-tune contrast selection in the second form and click "Apply"
            3. **Set Parameters**: Configure significance thresholds and gene selection method, then click "Apply"
            4. **View Results**: Switch to the heatmap tab to see results

            ### Gene Selection Methods

            **Frequent DEGs**: Automatically selects the most frequently differentially expressed genes across your selected contrasts
            - Uses the "Number of genes to display" setting
            - Genes are ranked by how often they appear as significant across contrasts

            **Custom**: Use your own gene list for targeted analysis
            - Enter gene symbols one per line in the text area
            - Number of genes setting does not apply
            - All entered genes will be used (if found in data)

            ### Threshold Application
            - **P-value and Log2FC thresholds** apply to both gene selection methods
            - **Number of genes** threshold only applies to Frequent DEGs method
            - Custom genes are filtered using the significance thresholds

            ### Important Notes
            - You must click "Apply" after making selections for changes to take effect
            - Selecting datasets automatically updates available contrasts
            - Each form works independently - apply changes as needed

            ### Form Coordination
            - Dataset selection populates available contrasts
            - Parameter changes affect gene selection and visualization
            - All changes require clicking the respective "Apply" button

            ---
            *Click the arrow above to close this help panel*
            """
        )

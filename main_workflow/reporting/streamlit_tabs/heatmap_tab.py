"""
Heatmap Tab for UORCA Explorer.

This tab displays interactive heatmaps showing log2 fold changes for selected genes across contrasts.
Features combined form for contrast and gene selection with proper select/clear all functionality.
"""

import logging
import streamlit as st
import pandas as pd
import traceback
from typing import List, Tuple, Optional, Dict, Any

from .helpers import (
    check_ai_generating,
    setup_fragment_decorator,
    cached_figure_creation,
    cached_identify_important_genes,
    log_streamlit_tab,
    log_streamlit_function,
    log_streamlit_event,
    log_streamlit_user_action,
    group_contrasts_by_organism,
    filter_genes_by_organism,
    get_organism_display_name
)
from ResultsIntegration import ResultsIntegrator

# Set up fragment decorator
setup_fragment_decorator()

logger = logging.getLogger(__name__)


@log_streamlit_tab("Heatmap")
def render_heatmap_tab(ri: ResultsIntegrator, selected_datasets: List[str], **kwargs):
    """
    Render the heatmap tab with combined contrast and gene selection form.

    Args:
        ri: ResultsIntegrator instance
        selected_datasets: List of selected dataset IDs from sidebar
        **kwargs: Additional arguments (maintained for compatibility)
    """
    st.header("Explore DEG Heatmap")
    st.markdown("**Interactive heatmap showing log2 fold changes for selected genes across contrasts.** Configure your analysis parameters below, then view the heatmap.")

    # Get selected datasets from sidebar
    if not selected_datasets:
        selected_datasets = st.session_state.get('selected_datasets_from_sidebar', [])

    if not selected_datasets:
        st.info("**Getting Started:**")
        st.markdown("""
        1. **Select Datasets** in the sidebar first
        2. **Choose Contrasts** using the form below (populated from your datasets)
        3. **Configure Gene Selection** using the gene selection form
        4. **View Results** - heatmap will be generated automatically

        **Tip:** Start with dataset selection in the sidebar to populate the forms below.
        """)
        return

    # Initialize session state for selected contrasts
    if 'selected_contrasts_heatmap' not in st.session_state:
        st.session_state['selected_contrasts_heatmap'] = set()

    # Check if datasets have changed and auto-select all contrasts by default
    current_datasets = set(selected_datasets)
    if 'previous_datasets_heatmap' not in st.session_state:
        st.session_state['previous_datasets_heatmap'] = set()

    if current_datasets != st.session_state['previous_datasets_heatmap']:
        # Datasets changed - select all contrasts by default
        contrast_data = _create_contrast_table_data_filtered(ri, selected_datasets)
        all_contrasts = set()
        for item in contrast_data:
            contrast_tuple = (item['Accession'], item['Contrast'])
            all_contrasts.add(contrast_tuple)
        st.session_state['selected_contrasts_heatmap'] = all_contrasts
        st.session_state['previous_datasets_heatmap'] = current_datasets

    # Quick selection buttons (outside form)
    st.subheader("1. Select Contrasts")
    contrast_data = _create_contrast_table_data_filtered(ri, selected_datasets)

    if contrast_data:
        col1, col2 = st.columns(2)
        with col1:
            if st.button("Select All Contrasts", key="select_all_contrasts_heatmap"):
                # Select all contrasts
                all_contrasts = set()
                for item in contrast_data:
                    contrast_tuple = (item['Accession'], item['Contrast'])
                    all_contrasts.add(contrast_tuple)
                st.session_state['selected_contrasts_heatmap'] = all_contrasts
                log_streamlit_user_action(f"Selected all {len(all_contrasts)} contrasts for heatmap")
                st.rerun()
        with col2:
            if st.button("Clear All Contrasts", key="clear_all_contrasts_heatmap"):
                # Clear all selections
                st.session_state['selected_contrasts_heatmap'] = set()
                log_streamlit_user_action("Cleared all contrast selections for heatmap")
                st.rerun()

    # Render combined form
    form_results = _render_combined_heatmap_form(ri, selected_datasets, contrast_data)

    if form_results:
        selected_contrasts = form_results['contrasts']
        gene_params = form_results['gene_params']

        # Convert accession-based contrasts to analysis_id-based contrasts
        analysis_contrasts = []
        if selected_contrasts:
            for accession, contrast_id in selected_contrasts:
                # Find the analysis_id for this accession
                analysis_id = None
                for aid in ri.analysis_info:
                    if ri.analysis_info[aid].get('accession') == accession:
                        analysis_id = aid
                        break

                if analysis_id:
                    analysis_contrasts.append((analysis_id, contrast_id))

        # Auto-select genes if we have contrasts and gene parameters
        gene_sel = []
        if analysis_contrasts and gene_params:
            gene_sel = _auto_select_genes(ri, ri.results_dir, analysis_contrasts, gene_params)

        # Display current selections summary
        if analysis_contrasts and gene_sel:
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Datasets", len(selected_datasets))
            with col2:
                st.metric("Contrasts", len(analysis_contrasts))
            with col3:
                st.metric("Genes", len(gene_sel))

            # Display gene selection method information
            if gene_params.get('gene_selection_method') == "Custom":
                custom_genes_count = len(gene_params.get('custom_genes', []))
                success_rate = (len(gene_sel) / custom_genes_count * 100) if custom_genes_count > 0 else 0
                if success_rate == 100:
                    st.success(f"🎯 **Custom Gene Selection**: All {len(gene_sel)} custom genes found in datasets")
                else:
                    st.info(f"🎯 **Custom Gene Selection**: Displaying {len(gene_sel)} of {custom_genes_count} custom genes ({success_rate:.0f}% found in datasets)")
            else:
                st.info(f"📊 **Frequent DEGs**: Displaying {len(gene_sel)} most frequently differentially expressed genes")

            # Group contrasts by organism and render heatmaps
            organism_groups = group_contrasts_by_organism(ri, analysis_contrasts)

            if len(organism_groups) == 1:
                # Single organism - no sub-tabs needed
                organism = list(organism_groups.keys())[0]
                organism_contrasts = organism_groups[organism]
                organism_genes = filter_genes_by_organism(ri, gene_sel, organism, organism_contrasts)

                _draw_heatmap(
                    ri,
                    organism_genes,
                    organism_contrasts,
                    gene_params.get('pvalue_thresh', 0.05),
                    gene_params.get('lfc_thresh', 1.0),
                    use_dynamic_filtering=True,
                    hide_empty_rows_cols=True,
                    font_size=12,
                    show_grid_lines=True,
                    grid_opacity=0.3
                )
            else:
                # Multiple organisms - create sub-tabs
                organism_names = list(organism_groups.keys())
                tab_names = [get_organism_display_name(org) for org in organism_names]
                tabs = st.tabs(tab_names)

                for i, organism in enumerate(organism_names):
                    with tabs[i]:
                        organism_contrasts = organism_groups[organism]
                        organism_genes = filter_genes_by_organism(ri, gene_sel, organism, organism_contrasts)

                        if organism_genes:
                            _draw_heatmap(
                                ri,
                                organism_genes,
                                organism_contrasts,
                                gene_params.get('pvalue_thresh', 0.05),
                                gene_params.get('lfc_thresh', 1.0),
                                use_dynamic_filtering=True,
                                hide_empty_rows_cols=True,
                                font_size=12,
                                show_grid_lines=True,
                                grid_opacity=0.3
                            )
                        else:
                            st.warning(f"No genes found for {get_organism_display_name(organism)} with current parameters.")
                            st.info("Try adjusting the significance thresholds in the gene selection form or selecting more datasets for this species.")

        elif analysis_contrasts and not gene_sel:
            st.warning("No genes selected. Please configure gene selection parameters above.")
        elif not analysis_contrasts:
            st.info("Please select contrasts using the form above.")


@log_streamlit_function
def _render_combined_heatmap_form(ri: ResultsIntegrator, selected_datasets: List[str], contrast_data: List[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    """Render the combined form for contrast and gene selection."""

    with st.form("heatmap_combined_form"):
        st.subheader("Analysis Configuration")

        # Contrast Selection Section
        st.markdown("**Contrasts**")

        if contrast_data:
            df = pd.DataFrame(contrast_data)
            # Pre-select based on session state
            df['Select'] = df.apply(
                lambda row: (row['Accession'], row['Contrast']) in st.session_state['selected_contrasts_heatmap'],
                axis=1
            )

            st.info(f"Available contrasts from {len(selected_datasets)} selected datasets")

            edited_df = st.data_editor(
                df,
                use_container_width=True,
                hide_index=True,
                column_config={
                    "Select": st.column_config.CheckboxColumn(
                        "",
                        help="Check to include this contrast in the heatmap",
                        default=False
                    ),
                    "Accession": st.column_config.TextColumn(
                        "Dataset",
                        help="GEO accession number",
                        width=100
                    ),
                    "Contrast": st.column_config.TextColumn(
                        "Contrast",
                        help="Contrast identifier",
                        width=200
                    ),
                    "Description": st.column_config.TextColumn(
                        "Description",
                        help="Contrast description",
                        width=400
                    )
                },
                key="heatmap_combined_contrast_table"
            )
        else:
            edited_df = pd.DataFrame()
            st.info("No contrasts available for selected datasets")

        # Gene Selection Section
        st.markdown("---")
        st.markdown("**Gene Selection**")

        # Gene selection method
        gene_selection_method = st.radio(
            "Choose gene selection method",
            options=["Frequent DEGs", "Custom"],
            key="heatmap_combined_gene_method",
            help="**Frequent DEGs**: Automatically finds genes that are consistently differentially expressed across multiple contrasts. **Custom**: Use your own gene list for targeted analysis.",
            horizontal=True
        )

        # Custom gene input (show immediately when Custom is selected)
        custom_genes_input = ""
        if gene_selection_method == "Custom":
            custom_genes_input = st.text_area(
                "Custom Gene List",
                height=150,
                placeholder="Enter one gene per line, e.g.:\nTP53\nEGFR\nMYC\nBRCA1",
                help="Enter gene symbols, one per line. Genes will be filtered by significance thresholds.",
                key="heatmap_combined_custom_genes"
            )

            # Show preview for custom genes
            if custom_genes_input.strip():
                custom_genes_list = [gene.strip() for gene in custom_genes_input.strip().split('\n') if gene.strip()]
                if custom_genes_list:
                    st.write(f"**Preview:** {len(custom_genes_list)} genes entered")
                    preview_text = ", ".join(custom_genes_list[:10])
                    if len(custom_genes_list) > 10:
                        preview_text += f", ... (+{len(custom_genes_list) - 10} more)"
                    st.caption(preview_text)

        # Gene count (for Frequent DEGs)
        if gene_selection_method == "Frequent DEGs":
            gene_count_input = st.text_input(
                "Number of genes to display",
                value="50",
                help="Number of top frequently differentially expressed genes to include",
                key="heatmap_combined_gene_count"
            )
        else:
            gene_count_input = "50"  # Default, not used for custom

        # Significance Thresholds
        st.markdown("**Significance Thresholds**")
        col1, col2 = st.columns(2)
        with col1:
            lfc_thresh = st.text_input(
                "Log2FC Threshold",
                value="1.0",
                help="Absolute log2 fold change threshold for significance",
                key="heatmap_combined_lfc_threshold"
            )
        with col2:
            pvalue_thresh = st.text_input(
                "P-value Threshold",
                value="0.05",
                help="Adjusted p-value threshold for significance",
                key="heatmap_combined_pvalue_threshold"
            )

        # Validation
        validation_error = False
        try:
            lfc_val = float(lfc_thresh)
            pval_val = float(pvalue_thresh)
        except ValueError:
            st.error("Please enter valid numeric values for thresholds")
            lfc_val, pval_val = 1.0, 0.05
            validation_error = True

        # Validate gene count
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
            gene_count = 50

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

        # Single submit button
        submitted = st.form_submit_button("Generate Heatmap Analysis", type="primary")

        if submitted and not validation_error:
            # Get selected contrasts from current form state
            selected_contrasts = []
            if not edited_df.empty:
                selected_rows = edited_df[edited_df["Select"] == True]
                selected_contrasts = [
                    (row["Accession"], row["Contrast"])
                    for _, row in selected_rows.iterrows()
                ]

            # Update session state with current selections
            st.session_state['selected_contrasts_heatmap'] = set(selected_contrasts)

            # Validation
            if not selected_contrasts:
                st.error("Please select at least one contrast")
                return None

            # Gene validation with detailed checking
            if gene_selection_method == "Custom":
                with st.expander("Gene Validation", expanded=False):
                    st.write(f"**Total genes entered:** {len(custom_genes_list)}")

                    # Check which genes are actually found in the selected contrasts
                    available_genes = set()
                    for accession, contrast_id in selected_contrasts:
                        # Find the analysis_id for this accession
                        analysis_id = None
                        for aid in ri.analysis_info:
                            if ri.analysis_info[aid].get('accession') == accession:
                                analysis_id = aid
                                break

                        if analysis_id and analysis_id in ri.deg_data and contrast_id in ri.deg_data[analysis_id]:
                            deg_df = ri.deg_data[analysis_id][contrast_id]
                            if 'Gene' in deg_df.columns:
                                available_genes.update(deg_df['Gene'].tolist())

                    found_genes = [gene for gene in custom_genes_list if gene in available_genes]
                    missing_genes = [gene for gene in custom_genes_list if gene not in available_genes]

                    if missing_genes:
                        st.warning(f"**{len(missing_genes)} genes not found in selected contrasts:**")
                        if len(missing_genes) <= 10:
                            st.write(", ".join(missing_genes))
                        else:
                            st.write(f"{', '.join(missing_genes[:10])}, ... (+{len(missing_genes)-10} more)")

                    if found_genes:
                        st.success(f"**{len(found_genes)} genes found and will be displayed**")
                    else:
                        st.error("No genes found in selected contrasts!")
                        return None

            if gene_selection_method == "Frequent DEGs":
                log_streamlit_user_action(f"Heatmap analysis: {len(selected_contrasts)} contrasts, LFC={lfc_val}, P={pval_val}, genes={gene_count}, method=Frequent DEGs")
            else:
                log_streamlit_user_action(f"Heatmap analysis: {len(selected_contrasts)} contrasts, LFC={lfc_val}, P={pval_val}, method=Custom, custom_genes={len(custom_genes_list)}")

            return {
                'contrasts': selected_contrasts,
                'gene_params': {
                    'lfc_thresh': lfc_val,
                    'pvalue_thresh': pval_val,
                    'gene_count': gene_count,
                    'gene_selection_method': gene_selection_method,
                    'custom_genes': custom_genes_list
                }
            }

    return None


@log_streamlit_function
def _create_contrast_table_data_filtered(ri: ResultsIntegrator, selected_datasets: List[str]) -> List[Dict[str, Any]]:
    """Create data for the contrast selection table, filtered by selected datasets."""
    contrast_data = []

    for analysis_id in selected_datasets:
        if analysis_id in ri.analysis_info:
            info = ri.analysis_info[analysis_id]
            accession = info.get("accession", analysis_id)

            for contrast in info.get("contrasts", []):
                contrast_id = contrast["name"]
                description = contrast.get("description", "")

                contrast_data.append({
                    "Select": False,
                    "Accession": accession,
                    "Contrast": contrast_id,
                    "Description": description
                })

    return contrast_data


@log_streamlit_function
def _auto_select_genes(
    ri: ResultsIntegrator,
    results_dir: str,
    selected_contrasts: List[Tuple[str, str]],
    gene_params: Dict[str, Any]
) -> List[str]:
    """Select genes based on gene parameters and selected contrasts (analysis_id format)."""
    if not selected_contrasts:
        return []

    gene_selection_method = gene_params.get('gene_selection_method', 'Frequent DEGs')

    try:
        # Group contrasts by organism (selected_contrasts should already be in analysis_id format)
        organism_groups = group_contrasts_by_organism(ri, selected_contrasts)

        if gene_selection_method == 'Custom':
            # Custom gene selection
            custom_genes = gene_params.get('custom_genes', [])
            if not custom_genes:
                return []

            # Collect all available genes from selected contrasts
            available_genes = set()
            for organism_contrasts in organism_groups.values():
                for analysis_id, contrast_id in organism_contrasts:
                    if analysis_id in ri.deg_data and contrast_id in ri.deg_data[analysis_id]:
                        deg_df = ri.deg_data[analysis_id][contrast_id]
                        if 'Gene' in deg_df.columns:
                            available_genes.update(deg_df['Gene'].tolist())

            # Filter custom genes to only those available
            available_custom_genes = [gene for gene in custom_genes if gene in available_genes]
            missing_genes = [gene for gene in custom_genes if gene not in available_genes]

            # Show validation in an expander
            with st.expander("Custom Gene Validation", expanded=bool(missing_genes)):
                st.write(f"**Total custom genes:** {len(custom_genes)}")
                st.write(f"**Found in selected contrasts:** {len(available_custom_genes)}")

                if missing_genes:
                    st.warning(f"**{len(missing_genes)} genes not found in selected contrasts:**")
                    if len(missing_genes) <= 15:
                        st.write(", ".join(missing_genes))
                    else:
                        st.write(f"{', '.join(missing_genes[:15])}, ... (+{len(missing_genes)-15} more)")

            log_streamlit_event(f"Custom gene selection: {len(available_custom_genes)} of {len(custom_genes)} genes available")
            return available_custom_genes

        else:
            # Frequent DEGs selection
            if len(organism_groups) == 1:
                # Single organism
                organism = list(organism_groups.keys())[0]
                organism_contrasts = organism_groups[organism]

                top_genes = cached_identify_important_genes(
                    results_dir=results_dir,
                    selected_contrasts=organism_contrasts,
                    top_frequent=gene_params.get('gene_count', 50),
                    p_value_threshold=gene_params['pvalue_thresh'],
                    lfc_threshold=gene_params['lfc_thresh']
                )

                log_streamlit_event(f"Auto-selected {len(top_genes)} frequent DEGs for {organism}")
                return top_genes

            else:
                # Multiple organisms - combine genes from each
                all_selected_genes = []
                gene_count_per_organism = gene_params.get('gene_count', 50)

                for organism, organism_contrasts in organism_groups.items():
                    organism_genes = cached_identify_important_genes(
                        results_dir=results_dir,
                        selected_contrasts=organism_contrasts,
                        top_frequent=gene_count_per_organism,
                        p_value_threshold=gene_params['pvalue_thresh'],
                        lfc_threshold=gene_params['lfc_thresh']
                    )
                    all_selected_genes.extend(organism_genes)

                log_streamlit_event(f"Auto-selected {len(all_selected_genes)} frequent DEGs across {len(organism_groups)} organisms")
                return all_selected_genes

    except Exception as e:
        logger.error(f"Error in gene selection: {e}")
        st.error(f"Error selecting genes: {str(e)}")
        return []


@st.fragment
@log_streamlit_function
def _draw_heatmap(
    ri: ResultsIntegrator,
    gene_selection: List[str],
    contrast_pairs: List[Tuple[str, str]],
    p_thresh: Optional[float],
    lfc_thresh_val: Optional[float],
    use_dynamic_filtering: bool,
    hide_empty_rows_cols: bool,
    font_size: int,
    show_grid_lines: bool,
    grid_opacity: float
):
    """
    Create and display the heatmap using fragment isolation.
    """
    # Skip execution if AI is currently generating
    if check_ai_generating():
        return

    with st.spinner("Generating heatmap..."):
        try:
            # Use cached figure creation for better performance
            p_threshold = p_thresh if use_dynamic_filtering else None
            lfc_threshold = lfc_thresh_val if use_dynamic_filtering else None

            fig = cached_figure_creation(
                "create_lfc_heatmap",
                ri.results_dir,
                gene_selection,
                contrast_pairs,
                None,  # output_file
                p_threshold,
                lfc_threshold,
                hide_empty_rows_cols,
                font_size,
                show_grid_lines,
                grid_opacity,
            )

            if fig:
                log_streamlit_event("Heatmap generated successfully")
                st.info("**Heatmap Tips:** Hover over cells to see contrast descriptions and gene information. Use the forms above to modify your analysis.")
                st.plotly_chart(fig, use_container_width=True)

                # Display information about filtered genes/contrasts if hide_empty_rows_cols is True
                if hide_empty_rows_cols:
                    _display_filtered_elements_info(ri)
            else:
                log_streamlit_event("Failed to generate heatmap")
                st.error("Could not generate heatmap. Please check your selections and try adjusting the parameters above.")

        except Exception as e:
            logger.error(f"Error generating heatmap: {str(e)}", exc_info=True)
            st.error(f"Error generating heatmap: {str(e)}")
            _display_heatmap_error_details(e)


@log_streamlit_function
def _display_filtered_elements_info(ri: ResultsIntegrator):
    """Display information about filtered genes and contrasts."""
    try:
        filtered_info = ri.get_last_heatmap_filtered_info()
        filtered_genes = filtered_info.get('genes', [])
        filtered_contrasts = filtered_info.get('contrasts', [])

        if filtered_genes or filtered_contrasts:
            with st.expander("Hidden Elements (No Significant Values)", expanded=False):
                if filtered_genes:
                    st.subheader("Hidden Genes")
                    st.write(f"**{len(filtered_genes)} genes** hidden because they have no significant values in any selected contrasts:")
                    genes_text = ", ".join(filtered_genes)
                    st.text_area("Genes:", value=genes_text, height=100, disabled=True, key="heatmap_filtered_genes_display")

                if filtered_contrasts:
                    st.subheader("Hidden Contrasts")
                    st.write(f"**{len(filtered_contrasts)} contrasts** hidden because they have no significant values for any selected genes:")
                    for contrast in filtered_contrasts:
                        st.write(f"• {contrast}")

    except Exception as e:
        logger.warning(f"Could not display filtered elements info: {e}")


@log_streamlit_function
def _display_heatmap_error_details(error: Exception):
    """Display detailed error information in an expandable section."""
    with st.expander("Heatmap Error Details", expanded=False):
        st.code(traceback.format_exc())

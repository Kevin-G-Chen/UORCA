"""
Heatmap Tab for UORCA Explorer.

This tab displays interactive heatmaps showing log2 fold changes for selected genes across contrasts.
Features combined form for contrast and gene selection with proper select/clear all functionality.
"""

import logging
import re
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
    get_organism_display_name,
    get_valid_contrasts_with_data,
    plotly_fig_to_pdf_bytes,
    generate_plot_filename
)
from ResultsIntegration import ResultsIntegrator
import io
import zipfile
import json
from datetime import datetime
import textwrap

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
    st.markdown("Construct a heatmap of selected genes and contrasts. The generated plots can be downloaded. The data and script to produce the heatmap can also be downloaded for further customisation.")

    # Get selected datasets from sidebar
    if not selected_datasets:
        selected_datasets = st.session_state.get('selected_datasets_from_sidebar', [])

    if not selected_datasets:
        st.markdown("To get started, select datasets from the sidebar on the left.")
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
    st.subheader("Select Contrasts")
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

        # Display current selections summary (removed metrics and gene selection info)
        if analysis_contrasts and gene_sel:
            pass  # Metrics and info messages removed

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

        # Always show both gene selection options
        # Gene count (for Frequent DEGs)
        gene_count_input = st.text_input(
            "Number of genes to display (Frequent DEGs)",
            value="50",
            help="Number of top frequently differentially expressed genes to include",
            key="heatmap_combined_gene_count"
        )

        # Custom gene input (always visible)
        custom_genes_input = st.text_area(
            "Custom Gene List",
            height=150,
            placeholder="Enter one gene per line, e.g.:\nTP53\nEGFR\nMYC\nBRCA1",
            help="Enter gene symbols, one per line. Genes will be filtered by significance thresholds.",
            key="heatmap_combined_custom_genes"
        )

        # Parse custom genes without preview display
        custom_genes_list = []
        if custom_genes_input.strip():
            custom_genes_list = [gene.strip() for gene in custom_genes_input.strip().split('\n') if gene.strip()]

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

        # Validate gene count (always validate since it's always shown)
        try:
            gene_count = int(gene_count_input)
            if gene_count <= 0:
                raise ValueError("Gene count must be positive")
        except ValueError:
            st.error("Please enter a valid positive number for gene count")
            gene_count = 50
            validation_error = True

        # Method-specific validation
        if gene_selection_method == "Custom":
            if not custom_genes_list:
                st.error("Please enter at least one gene for custom selection")
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

            # Gene validation with detailed checking (like expression plots)
            if gene_selection_method == "Custom":
                if not custom_genes_list:
                    st.error("Please enter at least one gene for custom selection")
                    return None

                # Show gene and contrast validation with detailed checking
                with st.expander("Gene and contrast validation", expanded=False):
                    st.write(f"**Total genes entered:** {len(custom_genes_list)}")
                    st.write(f"**Total contrasts selected:** {len(selected_contrasts)}")

                    st.write("")  # Add spacing

                    # Check which genes are actually found in the selected contrasts
                    available_genes = set()
                    significant_genes = set()
                    contrasts_with_data = []
                    contrasts_without_data = []
                    contrasts_with_no_significant_custom_genes = []

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
                                all_contrast_genes = set(deg_df['Gene'].tolist())
                                available_genes.update(all_contrast_genes)
                                contrasts_with_data.append(f"{accession}:{contrast_id}")

                                # Check if this contrast has any significant custom genes
                                contrast_has_significant_custom_genes = False
                                if 'adj.P.Val' in deg_df.columns and 'logFC' in deg_df.columns:
                                    significant_contrast_genes = deg_df[
                                        (deg_df['adj.P.Val'] < pval_val) &
                                        (abs(deg_df['logFC']) > lfc_val)
                                    ]['Gene'].tolist()
                                    significant_genes.update(significant_contrast_genes)

                                    # Check if any of our custom genes are significant in this contrast
                                    custom_genes_in_contrast = [gene for gene in custom_genes_list if gene in significant_contrast_genes]
                                    if custom_genes_in_contrast:
                                        contrast_has_significant_custom_genes = True

                                if not contrast_has_significant_custom_genes:
                                    contrasts_with_no_significant_custom_genes.append(f"{accession}:{contrast_id}")
                        else:
                            contrasts_without_data.append(f"{accession}:{contrast_id}")

                    # Categorise genes
                    found_genes = [gene for gene in custom_genes_list if gene in available_genes]
                    missing_genes = [gene for gene in custom_genes_list if gene not in available_genes]
                    found_but_filtered = [gene for gene in found_genes if gene not in significant_genes]
                    will_display = [gene for gene in found_genes if gene in significant_genes]

                    # Gene validation results
                    st.write("**Gene Analysis:**")

                    if will_display:
                        st.write(f"**{len(will_display)} genes will be displayed** in the heatmap")

                    if missing_genes:
                        st.write(f"\n**{len(missing_genes)} genes not found in selected contrasts:**")
                        if len(missing_genes) <= 10:
                            st.write(", ".join(missing_genes))
                        else:
                            st.write(f"{', '.join(missing_genes[:10])}, ... (+{len(missing_genes)-10} more)")

                    if found_but_filtered:
                        st.write(f"\n**{len(found_but_filtered)} genes found but will be filtered out** (never significant in selected contrasts):")
                        if len(found_but_filtered) <= 10:
                            st.write(", ".join(found_but_filtered))
                        else:
                            st.write(f"{', '.join(found_but_filtered[:10])}, ... (+{len(found_but_filtered)-10} more)")

                    # Contrast validation results
                    st.write(f"\n**Contrast Analysis:**")
                    st.write(f"**{len(contrasts_with_data)} contrasts have data available**")

                    if contrasts_without_data:
                        st.write(f"\n**{len(contrasts_without_data)} contrasts missing data:**")
                        for contrast in contrasts_without_data:
                            st.write(f"• {contrast}")

                    if contrasts_with_no_significant_custom_genes:
                        st.write(f"\n**{len(contrasts_with_no_significant_custom_genes)} contrasts have no significant custom genes:**")
                        for contrast in contrasts_with_no_significant_custom_genes:
                            st.write(f"• {contrast}")

                    if not will_display:
                        st.write(f"\n**ERROR:** No genes will be displayed! All custom genes are either missing or never significant.")
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
    """Create data for the contrast selection table using standardized validation logic."""
    # Use the centralized validation function
    valid_contrasts = get_valid_contrasts_with_data(ri, selected_datasets)

    # Convert to format expected by heatmap table
    contrast_data = []
    for contrast in valid_contrasts:
        contrast_data.append({
            "Select": False,
            "Accession": contrast["accession"],
            "Contrast": contrast["contrast_name"],
            "Description": contrast["description"]
        })
    # Sort by GEO accession numeric portion, then by contrast name
    def _acc_num(acc: str) -> int:
        m = re.search(r"(\d+)", str(acc) or "")
        return int(m.group(1)) if m else float('inf')

    contrast_data.sort(key=lambda x: (_acc_num(x["Accession"]), x["Accession"], x["Contrast"]))
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

            # Filter custom genes to only those available (removed validation expandable)
            available_custom_genes = [gene for gene in custom_genes if gene in available_genes]

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
                st.plotly_chart(fig, use_container_width=True)

                # Download buttons section
                st.markdown("---")
                st.subheader("Download Options")

                # Prepare reproducible download package immediately if heatmap data available
                try:
                    last_df = getattr(ri, 'last_heatmap_df', None)
                    last_ctx = getattr(ri, 'last_heatmap_context', None)

                    # Create two columns for download buttons
                    col1, col2 = st.columns(2)

                    # PDF download button in first column
                    with col1:
                        pdf_bytes = plotly_fig_to_pdf_bytes(fig, width=1600, height=1200, scale=2)
                        if pdf_bytes:
                            st.download_button(
                                label="Download as PDF",
                                data=pdf_bytes,
                                file_name=generate_plot_filename("heatmap", "pdf"),
                                mime="application/pdf",
                                help="Download the heatmap as a high-resolution PDF suitable for publications"
                            )

                    # Reproducible package button in second column
                    with col2:
                        if last_df is not None and last_ctx is not None:
                            # Build in-memory ZIP immediately so users can download right away
                            zip_buffer = io.BytesIO()
                            with zipfile.ZipFile(zip_buffer, mode='w', compression=zipfile.ZIP_DEFLATED) as zf:
                                # heatmap data CSV
                                csv_bytes = last_df.reset_index().to_csv(index=False).encode('utf-8')
                                zf.writestr('heatmap_data.csv', csv_bytes)

                                # metadata JSON - include x-axis default labels mapping
                                try:
                                    clustered = last_ctx.get('clustered_contrasts') or []
                                    contrast_labels = last_ctx.get('contrast_labels') or []
                                    simplified = last_ctx.get('simplified_labels') or []
                                    # Build x-axis label mapping: contrast_label -> simplified_label (if available)
                                    x_axis_labels = {c: s for c, s in zip(contrast_labels, simplified) if c in clustered}
                                except Exception:
                                    x_axis_labels = {}

                                meta = {
                                    'generated_at': datetime.utcnow().isoformat() + 'Z',
                                    'p_value_threshold': last_ctx.get('p_value_threshold'),
                                    'lfc_threshold': last_ctx.get('lfc_threshold'),
                                    'genes': last_ctx.get('genes'),
                                    'contrast_labels': last_ctx.get('contrast_labels'),
                                    'clustered_genes': last_ctx.get('clustered_genes'),
                                    'clustered_contrasts': last_ctx.get('clustered_contrasts'),
                                    'x_axis_labels': x_axis_labels,
                                    # Default plotting parameters for reproduce script
                                    'output_format': 'png',
                                    'dpi': 300,
                                    'width': 12.0,
                                    'height_per_gene': 0.15,
                                    'font_size': 10.0,
                                    'left_margin': 2.5,
                                    'bottom_margin': 1.5
                                }
                                zf.writestr('metadata.json', json.dumps(meta, indent=2))

                                # reproducible script (static output: PNG/PDF/SVG, supports label mapping and layout params)
                                repro_script = _build_repro_script_static()
                                zf.writestr('reproduce_heatmap.py', repro_script)

                                # README describing files and usage
                                readme_text = _build_readme_text()
                                zf.writestr('README.txt', readme_text)

                            zip_buffer.seek(0)

                            file_name = f"uorca_heatmap_repro_{datetime.utcnow().strftime('%Y%m%dT%H%M%SZ')}.zip"
                            # Rename button and add descriptive text
                            st.download_button(
                                label="Reproduce and edit heatmap",
                                data=zip_buffer.getvalue(),
                                file_name=file_name,
                                mime='application/zip',
                                help="Download a reproducible package (CSV, metadata, editable script, and README). Adjust labels, fonts, margins, or export a high-resolution PDF for publication. Run: python reproduce_heatmap.py --input heatmap_data.csv --output heatmap.pdf --output-format pdf --dpi 300 --font-size 10"
                            )
                except Exception as e:
                    logger.debug(f"Could not prepare reproducible package: {e}")

                # Filtering information is now included in the form validation
            else:
                log_streamlit_event("Failed to generate heatmap")
                st.error("Could not generate heatmap. Please check your selections and try adjusting the parameters above.")

        except Exception as e:
            logger.error(f"Error generating heatmap: {str(e)}", exc_info=True)
            st.error(f"Error generating heatmap: {str(e)}")
            _display_heatmap_error_details(e)


@log_streamlit_function
def _display_comprehensive_filtering_info(ri: ResultsIntegrator, requested_genes: List[str], selected_contrasts: List[Tuple[str, str]]):
    """Display comprehensive information about genes and contrasts removed from the heatmap."""
    # This function is now integrated into the form validation and no longer displays separately
    pass

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


def _build_repro_script() -> str:
    """Return contents of a standalone script that reproduces the heatmap from CSV."""
    # kept for backward compatibility but not used
    return _build_repro_script_static()


def _build_repro_script_static() -> str:
    """Return contents of a standalone script that reproduces the heatmap as a static PNG.

    The script supports an optional JSON file to remap x-axis labels:
      python reproduce_heatmap.py --input heatmap_data.csv --output heatmap.png --xlabel-mapping my_labels.json
    where my_labels.json is a dict mapping original column names to desired labels.
    """
    script_lines = [
    '#!/usr/bin/env python3',
    '"""',
    'Reproduce a static heatmap PNG from `heatmap_data.csv`.',
    '',
    'Usage:',
    '  python reproduce_heatmap.py --input heatmap_data.csv --output heatmap.png --xlabel-mapping labels.json --dpi 300',
    '',
    'The optional `labels.json` should be a JSON object mapping existing column names to new labels,',
    'e.g. {"analysis1_contrastA": "Control vs Treated", "analysis2_contrastB": "Condition X"}',
    '"""',
    'import argparse',
    'import json',
    'import pandas as pd',
    'import numpy as np',
    '',
    '',
    'def main():',
    '    parser = argparse.ArgumentParser()',
    "    parser.add_argument('--input', default='heatmap_data.csv')",
    "    parser.add_argument('--output', default='heatmap.png')",
    "    parser.add_argument('--output-format', choices=['png','pdf','svg'], default=None, help='Output file format; if omitted, taken from metadata.json or defaults to png')",
    "    parser.add_argument('--xlabel-mapping', default=None, help='JSON file mapping original column names to new labels')",
    "    parser.add_argument('--dpi', type=int, default=None)",
    "    parser.add_argument('--width', type=float, default=None, help='Figure width in inches')",
    "    parser.add_argument('--height-per-gene', type=float, default=None, help='Height per gene in inches')",
    "    parser.add_argument('--font-size', type=float, default=None, help='Base font size for labels')",
    "    parser.add_argument('--left-margin', type=float, default=None, help='Left margin in inches')",
    "    parser.add_argument('--bottom-margin', type=float, default=None, help='Bottom margin in inches')",
    '    args = parser.parse_args()',
    '',
    '    # Load metadata.json if present and use it as defaults for parameters',
    '    meta = {}',
    '    try:',
    "        with open('metadata.json', 'r') as mfh:",
    '            meta = json.load(mfh)',
    '    except Exception:',
    '        meta = {}',
    '',
    '    # Resolve parameters: CLI overrides metadata; metadata overrides hardcoded defaults',
    "    out_fmt = args.output_format or meta.get('output_format', 'png')",
    "    dpi = args.dpi or meta.get('dpi', 300)",
    "    width = args.width or meta.get('width', 12.0)",
    "    height_per_gene = args.height_per_gene or meta.get('height_per_gene', 0.15)",
    "    font_size = args.font_size or meta.get('font_size', 10.0)",
    "    left_margin = args.left_margin or meta.get('left_margin', 2.5)",
    "    bottom_margin = args.bottom_margin or meta.get('bottom_margin', 1.5)",
    '',
    '    # Read input data',
    '    df = pd.read_csv(args.input)',
    "    if 'Gene' in df.columns:",
    "        df = df.set_index('Gene')",
    '',
    '    # Apply xlabel mapping: CLI file if provided, otherwise metadata x_axis_labels',
    '    xlabel_map = {}',
    '    if args.xlabel_mapping:',
    "        with open(args.xlabel_mapping, 'r') as fh:",
    '            xlabel_map = json.load(fh)',
    '    else:',
    "        xlabel_map = meta.get('x_axis_labels', {})",
    '',
    '    # Reorder columns if metadata provides clustered_contrasts',
    "    clustered = meta.get('clustered_contrasts')",
    '    if clustered and set(clustered).issubset(set(df.columns)):',
    '        df = df[clustered]',
    '',
    '    # Replace column labels according to mapping, preserving order',
    '    new_cols = [xlabel_map.get(c, c) for c in df.columns]',
    '',
    '    # Create figure size based on number of genes and requested width',
    '    n_genes = df.shape[0]',
    '    fig_height = max(4, n_genes * height_per_gene)',
    '',
    '    # Perform plotting inside a guarded block so errors are written to an error report',
    '    try:',
    '        # Now import plotting libraries (delayed to allow --help without heavy deps)',
    '        import matplotlib',
    "        matplotlib.use('Agg')",
    '        import matplotlib.pyplot as plt',
    '        import seaborn as sns',
    '',
    '        fig, ax = plt.subplots(figsize=(width, fig_height))',
    '',
    '        # Use diverging palette centered at zero',
    '        max_abs = max(np.nanmax(np.abs(df.values)), 1.0)',
    '        # Adjust fonts and margins',
    "        plt.rcParams.update({'font.size': font_size})",
    '',
    '        sns.heatmap(',
    '            df,',
    "            cmap='RdBu_r',",
    '            vmin=-max_abs,',
    '            vmax=max_abs,',
    "            cbar_kws={'label': 'Log2FC'},",
    '            xticklabels=new_cols,',
    '            yticklabels=True,',
    '            ax=ax',
    '        )',
    '',
    "        plt.xticks(rotation=45, ha='right', fontsize=font_size)",
    '        plt.yticks(fontsize=font_size)',
    '',
    '        # Convert margins in inches to figure fraction',
    '        left = left_margin / width',
    '        bottom = bottom_margin / fig_height',
    '        plt.subplots_adjust(left=left, right=0.98, top=0.98, bottom=bottom)',
    '',
    '        # Save in requested format',
    '        ofmt = out_fmt.lower()',
    "        if ofmt == 'pdf':",
    "            fig.savefig(args.output, dpi=dpi, format='pdf')",
    "        elif ofmt == 'svg':",
    "            fig.savefig(args.output, dpi=dpi, format='svg')",
    '        else:',
    '            fig.savefig(args.output, dpi=dpi)',
    '',
    '    except Exception as exc:',
    '        # Write an error report so users can inspect failures when running outside the app',
    '        import traceback',
    '        report = {',
    "            'error': str(exc),",
    "            'traceback': traceback.format_exc(),",
    "            'note': 'Ensure matplotlib and seaborn are installed or run via the UORCA uv environment (e.g. `uv run python reproduce_heatmap.py ...`).'",
    '        }',
    '        try:',
    "            with open('error_report.txt', 'w') as ef:",
    "                ef.write('Error reproducing heatmap\\n')",
    "                ef.write('Error: ' + report['error'] + '\\n\\n')",
    "                ef.write(report['traceback'])",
    "                ef.write('\\n\\n')",
    "                ef.write(report['note'])",
    '        except Exception:',
    '            pass',
    '        # Re-raise so CLI returns non-zero status',
    '        raise',
    '',
    '',
    "if __name__ == '__main__':",
    "    main()",
    "",
    ]
    return '\n'.join(script_lines) + '\n'


def _build_readme_text() -> str:
    return (
        "UORCA Heatmap Reproducible Package\n"
        "--------------------------------\n\n"
        "Included files:\n"
        "- heatmap_data.csv: CSV table of genes (rows) and contrasts (columns) containing log2 fold change values (non-significant values set to 0).\n"
        "- metadata.json: JSON metadata including thresholds, clustered contrast order, and default x-axis label mapping.\n"
        "- reproduce_heatmap.py: Standalone script to generate a static heatmap (PNG/PDF/SVG). Supports options to remap x-axis labels, adjust font sizes, margins, output format and DPI.\n"
        "- README.txt: This file.\n\n"
        "Quick start:\n"
        "1. Unzip the package.\n"
        "2. (Optional) Edit or create a JSON file mapping original contrast column names to publication-friendly labels, e.g.:\n"
        "   {\"GSE111143_contrast1\": \"Control vs Treated\", \"GSE222222_contrastA\": \"Condition X\"}\n"
        "3. Run the script to produce a high-resolution PDF for publication: \n"
        "   python reproduce_heatmap.py --input heatmap_data.csv --output heatmap.pdf --output-format pdf --dpi 300 --xlabel-mapping my_labels.json --font-size 10 --left-margin 2.5\n\n"
        "Why download?\n"
        "- Edit axis labels for clarity in figures.\n"
        "- Produce vector PDF/SVG outputs suitable for publication.\n"
        "- Adjust fonts and margins for journal layouts.\n\n"
        "If you need the original DEG/CPM source files included for provenance, ask to include them in the package."
    )

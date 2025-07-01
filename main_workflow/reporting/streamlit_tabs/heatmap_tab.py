"""
Heatmap Tab for UORCA Explorer.

This tab displays interactive heatmaps showing log2 fold changes for selected genes across contrasts.
"""

import logging
import streamlit as st
import traceback
from typing import List, Tuple, Optional

from .helpers import (
    check_ai_generating,
    setup_fragment_decorator,
    cached_figure_creation,
    log_streamlit_tab,
    log_streamlit_function,
    log_streamlit_event,
    group_contrasts_by_organism,
    filter_genes_by_organism,
    get_organism_display_name
)
from ResultsIntegration import ResultsIntegrator

# Set up fragment decorator
setup_fragment_decorator()

logger = logging.getLogger(__name__)


@log_streamlit_tab("Heatmap")
def render_heatmap_tab(
    ri: ResultsIntegrator,
    gene_sel: List[str],
    selected_contrasts: List[Tuple[str, str]],
    effective_pvalue_thresh: Optional[float],
    effective_lfc_thresh: Optional[float],
    use_dynamic_filtering: bool,
    hide_empty_rows_cols: bool,
    gene_selection_method: str = "Frequent DEGs",
    custom_genes_count: int = 0
):
    """
    Render the heatmap tab.

    Args:
        ri: ResultsIntegrator instance
        gene_sel: List of selected genes
        selected_contrasts: List of (analysis_id, contrast_id) tuples
        effective_pvalue_thresh: P-value threshold for filtering
        effective_lfc_thresh: Log fold change threshold for filtering
        use_dynamic_filtering: Whether to use dynamic filtering
        hide_empty_rows_cols: Whether to hide empty rows/columns
        gene_selection_method: Method used for gene selection ("Frequent DEGs" or "Custom")
        custom_genes_count: Number of custom genes provided (if applicable)
    """
    st.header("Explore DEG Heatmap")
    st.markdown("**Interactive heatmap showing log2 fold changes for selected genes across contrasts.** Hover over cells for details. Use the Dataset & Contrast Selection and Heatmap Parameters forms in the sidebar.")

    # Display settings for heatmap
    display_settings = _render_display_settings()

    if not selected_contrasts:
        log_streamlit_event("No contrasts selected for heatmap")
        st.info("**Getting Started:**")
        st.markdown("""
        1. **Select Datasets** in the "Dataset & Contrast Selection" section in the sidebar
        2. **Select Contrasts** using the second form (auto-populated from your datasets)
        3. **Configure Parameters** in the "Heatmap Parameters" section (thresholds and gene count)
        4. **View Results** - genes are automatically selected and heatmap generated

        **Tip:** Start with dataset selection - this will automatically populate available contrasts.
        """)
    elif not gene_sel:
        log_streamlit_event("No genes selected for heatmap")
        if gene_selection_method == "Custom":
            st.warning("No custom genes found in the selected datasets. Please check your gene list and ensure the genes are present in your data.")
        else:
            st.warning("No genes were automatically selected with the current parameters. Try adjusting the significance thresholds in the sidebar form.")
    else:
        # Display gene selection method information
        if gene_selection_method == "Custom":
            st.info(f"🎯 **Custom Gene Selection**: Displaying {len(gene_sel)} custom genes (from {custom_genes_count} provided)")
        else:
            st.info(f"📊 **Frequent DEGs**: Displaying {len(gene_sel)} most frequently differentially expressed genes")
        log_streamlit_event(f"Heatmap: {len(gene_sel)} genes, {len(selected_contrasts)} contrasts")

        # Group contrasts by organism
        organism_groups = group_contrasts_by_organism(ri, selected_contrasts)

        if len(organism_groups) == 1:
            # Single organism - no sub-tabs needed
            organism = list(organism_groups.keys())[0]
            organism_contrasts = organism_groups[organism]
            organism_genes = filter_genes_by_organism(ri, gene_sel, organism, organism_contrasts)

            _draw_heatmap(
                ri,
                organism_genes,
                organism_contrasts,
                effective_pvalue_thresh,
                effective_lfc_thresh,
                use_dynamic_filtering,
                hide_empty_rows_cols,
                **display_settings
            )
        else:
            # Multiple organisms - create sub-tabs
            # Create organism-specific tabs
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
                            effective_pvalue_thresh,
                            effective_lfc_thresh,
                            use_dynamic_filtering,
                            hide_empty_rows_cols,
                            **display_settings
                        )
                    else:
                        st.warning(f"No genes found for {get_organism_display_name(organism)} with current parameters.")
                        st.info("Try adjusting the significance thresholds in the sidebar or selecting more datasets for this species.")


@log_streamlit_function
def _render_display_settings() -> dict:
    """Render display settings controls in the sidebar and return the settings."""
    with st.sidebar.expander("Display Settings", expanded=False):
        heatmap_font_size = st.slider("Font size", 8, 16, 12, key="heatmap_tab_font")
        show_grid_lines = st.checkbox("Show grid lines", value=True, key="heatmap_tab_grid")
        grid_opacity = st.slider("Grid opacity", 0.1, 1.0, 0.3, key="heatmap_tab_grid_opacity")

    return {
        "font_size": heatmap_font_size,
        "show_grid_lines": show_grid_lines,
        "grid_opacity": grid_opacity,
    }


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

    Args:
        ri: ResultsIntegrator instance
        gene_selection: List of selected genes
        contrast_pairs: List of (analysis_id, contrast_id) tuples
        p_thresh: P-value threshold for filtering
        lfc_thresh_val: Log fold change threshold for filtering
        use_dynamic_filtering: Whether to use dynamic filtering
        hide_empty_rows_cols: Whether to hide empty rows/columns
        font_size: Font size for heatmap
        show_grid_lines: Whether to show grid lines
        grid_opacity: Opacity of grid lines
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
                _display_heatmap_info()
                st.plotly_chart(fig, use_container_width=True)

                # Display information about filtered genes/contrasts if hide_empty_rows_cols is True
                if hide_empty_rows_cols:
                    _display_filtered_elements_info(ri)
            else:
                log_streamlit_event("Failed to generate heatmap")
                st.error("Could not generate heatmap. Please check your selections and try adjusting the parameters in the sidebar form.")

        except Exception as e:
            logger.error(f"Error generating heatmap: {str(e)}", exc_info=True)
            st.error(f"Error generating heatmap: {str(e)}")
            _display_heatmap_error_details(e)


@log_streamlit_function
def _display_heatmap_info():
    """Display informational messages about the heatmap."""
    st.info("**Heatmap Tips:** Hover over cells to see contrast descriptions and gene information. Use the Dataset & Contrast Selection and Heatmap Parameters forms in the sidebar to modify your analysis.")


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
                    # Display genes in a more compact format
                    genes_text = ", ".join(filtered_genes)
                    st.text_area("Genes:", value=genes_text, height=100, disabled=True, key="heatmap_filtered_genes_display")

                if filtered_contrasts:
                    st.subheader("Hidden Contrasts")
                    st.write(f"**{len(filtered_contrasts)} contrasts** hidden because they have no significant values for any selected genes:")
                    for contrast in filtered_contrasts:
                        st.write(f"• {contrast}")

                if not filtered_genes and not filtered_contrasts:
                    st.info("No elements were hidden - all genes and contrasts have at least one significant value.")
    except Exception as e:
        logger.warning(f"Could not display filtered elements info: {e}")


@log_streamlit_function
def _display_heatmap_error_details(error: Exception):
    """Display detailed error information in an expandable section."""
    with st.expander("Heatmap Error Details", expanded=False):
        st.code(traceback.format_exc())

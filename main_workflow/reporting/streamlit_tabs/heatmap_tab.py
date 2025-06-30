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
    hide_empty_rows_cols: bool
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
        st.warning("No genes were automatically selected with the current parameters. Try adjusting the significance thresholds in the sidebar form.")
    else:
        log_streamlit_event(f"Heatmap: {len(gene_sel)} genes, {len(selected_contrasts)} contrasts")

        # Group contrasts by organism
        organism_groups = group_contrasts_by_organism(ri, selected_contrasts)

        if len(organism_groups) == 1:
            # Single organism - no sub-tabs needed
            organism = list(organism_groups.keys())[0]
            organism_contrasts = organism_groups[organism]
            organism_genes = filter_genes_by_organism(ri, gene_sel, organism, organism_contrasts)

            st.success(f"Ready to display heatmap with {len(organism_genes)} genes across {len(organism_contrasts)} contrasts ({get_organism_display_name(organism)})")

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
            st.success(f"Found {len(organism_groups)} species in selected data. Creating species-specific heatmaps to prevent gene name conflicts.")

            # Show organism breakdown
            with st.expander("Species Breakdown", expanded=True):
                for organism, contrasts in organism_groups.items():
                    datasets = set(contrast[0] for contrast in contrasts)
                    st.write(f"**{get_organism_display_name(organism)}**: {len(datasets)} datasets, {len(contrasts)} contrasts")

            st.info("**Important**: Gene expression data is separated by species to ensure biological accuracy and prevent cross-species gene name conflicts.")

            # Create organism-specific tabs
            organism_names = list(organism_groups.keys())
            tab_names = [get_organism_display_name(org) for org in organism_names]
            tabs = st.tabs(tab_names)

            for i, organism in enumerate(organism_names):
                with tabs[i]:
                    organism_contrasts = organism_groups[organism]
                    organism_genes = filter_genes_by_organism(ri, gene_sel, organism, organism_contrasts)

                    if organism_genes:
                        st.success(f"**{get_organism_display_name(organism)} Analysis**")
                        st.write(f"Displaying {len(organism_genes)} genes across {len(organism_contrasts)} contrasts")

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
        heatmap_font_size = st.slider("Font size", 8, 16, 12, key="heatmap_font")
        show_grid_lines = st.checkbox("Show grid lines", value=True, key="heatmap_grid")
        grid_opacity = st.slider("Grid opacity", 0.1, 1.0, 0.3, key="heatmap_grid_opacity")

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

                # Show current configuration
                st.caption(f"Showing {len(gene_selection)} genes across {len(contrast_pairs)} contrasts with P-value ≤ {p_thresh:.3f} and |Log2FC| ≥ {lfc_thresh_val:.1f}")
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
def _display_heatmap_error_details(error: Exception):
    """Display detailed error information in an expandable section."""
    with st.expander("Heatmap Error Details", expanded=False):
        st.code(traceback.format_exc())

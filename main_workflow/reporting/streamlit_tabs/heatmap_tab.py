"""
Heatmap Tab for UORCA Explorer.

This tab displays interactive heatmaps showing log2 fold changes for selected genes across contrasts.
"""

import logging
import streamlit as st
import traceback
from typing import List, Tuple, Optional

from .helpers import check_ai_generating, setup_fragment_decorator, cached_figure_creation
from ResultsIntegration import ResultsIntegrator

# Set up fragment decorator
setup_fragment_decorator()

logger = logging.getLogger(__name__)


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
    st.header("🌡️ Explore DEG Heatmap")
    st.markdown("**📊 Interactive heatmap showing log2 fold changes for selected genes across contrasts.** Hover over cells for details. Use the sidebar to adjust significance thresholds and filtering options.")

    # Display settings for heatmap
    display_settings = _render_display_settings()

    if not gene_sel:
        st.info("Please select genes from the sidebar.")
    elif not selected_contrasts:
        st.info("Please select contrasts in the 'Selections' tab.")
    else:
        # Create and display the heatmap
        _draw_heatmap(
            ri,
            gene_sel,
            selected_contrasts,
            effective_pvalue_thresh,
            effective_lfc_thresh,
            use_dynamic_filtering,
            hide_empty_rows_cols,
            **display_settings
        )


def _render_display_settings() -> dict:
    """Render display settings controls in the sidebar and return the settings."""
    with st.sidebar.expander("🎨 Display Settings", expanded=False):
        heatmap_font_size = st.slider("Font size", 8, 16, 12, key="heatmap_font")
        show_grid_lines = st.checkbox("Show grid lines", value=True, key="heatmap_grid")
        grid_opacity = st.slider("Grid opacity", 0.1, 1.0, 0.3, key="heatmap_grid_opacity")

    return {
        "font_size": heatmap_font_size,
        "show_grid_lines": show_grid_lines,
        "grid_opacity": grid_opacity,
    }


@st.fragment
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
            # Try to use cached figure creation if available
            if 'cached_figure_creation' in globals():
                fig = cached_figure_creation(
                    "create_lfc_heatmap",
                    ri,
                    gene_selection,
                    contrast_pairs,
                    None,
                )
            else:
                # Use direct creation with all parameters
                p_threshold = p_thresh if use_dynamic_filtering else None
                lfc_threshold = lfc_thresh_val if use_dynamic_filtering else None

                fig = ri.create_lfc_heatmap(
                    genes=gene_selection,
                    contrasts=contrast_pairs,
                    output_file=None,
                    p_value_threshold=p_threshold,
                    lfc_threshold=lfc_threshold,
                    hide_empty_rows_cols=hide_empty_rows_cols,
                    font_size=font_size,
                    show_grid_lines=show_grid_lines,
                    grid_opacity=grid_opacity,
                )

            if fig:
                _display_heatmap_info()
                st.plotly_chart(fig, use_container_width=True)
            else:
                st.error("Could not generate heatmap. Please check your selections.")

        except Exception as e:
            logger.error(f"Error generating heatmap: {str(e)}", exc_info=True)
            st.error(f"Error generating heatmap: {str(e)}")
            _display_heatmap_error_details(e)


def _display_heatmap_info():
    """Display informational messages about the heatmap."""
    info_messages = [
        "💡 Hover over cells in the heatmap to see contrast descriptions and gene information."
    ]
    for msg in info_messages:
        st.info(msg)


def _display_heatmap_error_details(error: Exception):
    """Display detailed error information in an expandable section."""
    with st.expander("🔍 Heatmap Error Details", expanded=False):
        st.code(traceback.format_exc())

"""
Expression Plots Tab for UORCA Explorer.

This tab displays violin plots showing gene expression distributions across sample groups.
"""

import logging
import streamlit as st
import traceback
from typing import List, Tuple, Dict, Any

from .helpers import (
    check_ai_generating,
    setup_fragment_decorator,
    cached_figure_creation,
    calculate_pagination_info,
    safe_rerun
)
from ResultsIntegration import ResultsIntegrator

# Set up fragment decorator
setup_fragment_decorator()

logger = logging.getLogger(__name__)


def render_expression_plots_tab(
    ri: ResultsIntegrator,
    gene_sel: List[str],
    selected_datasets: List[str],
    hide_x_labels: bool
):
    """
    Render the expression plots tab.

    Args:
        ri: ResultsIntegrator instance
        gene_sel: List of selected genes
        selected_datasets: List of selected dataset IDs
        hide_x_labels: Whether to hide x-axis labels
    """
    st.header("📈 Plot Gene Expression")
    st.markdown("**🎻 Violin plots showing gene expression distributions across sample groups.** Each panel represents one gene, with samples grouped by experimental conditions.")

    # Display settings for expression plots
    display_settings = _render_display_settings()

    if not gene_sel:
        st.info("Please select genes from the sidebar.")
    elif not selected_datasets:
        st.info("Please select datasets in the 'Selections' tab.")
    else:
        # Calculate pagination information
        total_pages, current_page, genes_per_page, current_genes = calculate_pagination_info(gene_sel)

        # Add pagination controls at the top if needed
        if total_pages > 1:
            _render_top_pagination_controls(current_page, total_pages)

        # Create and display the expression plots
        _draw_expression_plots(
            ri,
            gene_sel,
            selected_datasets,
            hide_x_labels,
            current_page,
            total_pages,
            **display_settings
        )


def _render_display_settings() -> Dict[str, Any]:
    """Render display settings controls in the sidebar and return the settings."""
    with st.sidebar.expander("🎨 Display Settings", expanded=False):
        facet_font_size = st.slider("Facet title size", 8, 16, 10, key="violin_font")
        lock_y_axis = st.checkbox("Lock y-axis across genes", value=False, key="violin_lock_y")
        show_raw_points = st.checkbox("Show raw points", value=True, key="violin_points")
        legend_position = st.selectbox("Legend position", ["Bottom", "Right", "Top"], index=0, key="violin_legend")

    return {
        "facet_font_size": facet_font_size,
        "lock_y_axis": lock_y_axis,
        "show_raw_points": show_raw_points,
        "legend_position": legend_position,
    }


def _render_top_pagination_controls(current_page: int, total_pages: int):
    """Render pagination controls at the top of the plots for convenience."""
    cols = st.columns([2, 1, 1, 1, 2])

    with cols[1]:
        prev_disabled = current_page <= 1
        if st.button("◀ Previous", disabled=prev_disabled, key="prev_main"):
            st.session_state.page_num = max(1, current_page - 1)
            safe_rerun()

    with cols[2]:
        st.markdown(f"**Page {current_page}/{total_pages}**")

    with cols[3]:
        next_disabled = current_page >= total_pages
        if st.button("Next ▶", disabled=next_disabled, key="next_main"):
            st.session_state.page_num = min(total_pages, current_page + 1)
            safe_rerun()


@st.fragment
def _draw_expression_plots(
    ri: ResultsIntegrator,
    gene_selection: List[str],
    dataset_selection: List[str],
    hide_labels: bool,
    page_num: int,
    total_pgs: int,
    facet_font_size: int,
    lock_y_axis: bool,
    show_raw_points: bool,
    legend_position: str
):
    """
    Create and display the expression plots using fragment isolation.

    Args:
        ri: ResultsIntegrator instance
        gene_selection: List of selected genes
        dataset_selection: List of selected datasets
        hide_labels: Whether to hide x-axis labels
        page_num: Current page number
        total_pgs: Total number of pages
        facet_font_size: Font size for facet titles
        lock_y_axis: Whether to lock y-axis across genes
        show_raw_points: Whether to show raw data points
        legend_position: Position of the legend
    """
    # Skip execution if AI is currently generating
    if check_ai_generating():
        return

    with st.spinner("Generating expression plots..."):
        try:
            # Calculate gene slice for the current page
            genes_per_page = 30
            start_idx = (page_num - 1) * genes_per_page
            end_idx = min(start_idx + genes_per_page, len(gene_selection))
            current_genes = gene_selection[start_idx:end_idx]

            # Try to use cached figure creation if available
            if 'cached_figure_creation' in globals():
                fig2 = cached_figure_creation(
                    "create_expression_plots",
                    ri,
                    current_genes,
                    dataset_selection,
                    "violin",
                    None,
                    hide_labels,
                    page_num,
                    facet_font_size,
                    lock_y_axis,
                    show_raw_points,
                    legend_position,
                    True,  # show_grid_lines
                    0.3    # grid_opacity
                )
            else:
                # Use direct creation
                fig2 = ri.create_expression_plots(
                    genes=current_genes,
                    analyses=dataset_selection,
                    plot_type="violin",
                    output_file=None,
                    hide_x_labels=hide_labels,
                    page_number=page_num,
                    facet_font_size=facet_font_size,
                    lock_y_axis=lock_y_axis,
                    show_raw_points=show_raw_points,
                    legend_position=legend_position,
                    show_grid_lines=True,
                    grid_opacity=0.3
                )

            if fig2:
                st.plotly_chart(fig2, use_container_width=True)
            else:
                st.error("Could not generate expression plots. Please check your selections.")

        except Exception as e:
            logger.error(f"Error generating expression plots: {str(e)}", exc_info=True)
            st.error(f"Error generating expression plots: {str(e)}")
            _display_expression_plot_error_details(e)


def _display_expression_plot_error_details(error: Exception):
    """Display detailed error information in an expandable section."""
    with st.expander("🔍 Expression Plot Error Details", expanded=False):
        st.code(traceback.format_exc())

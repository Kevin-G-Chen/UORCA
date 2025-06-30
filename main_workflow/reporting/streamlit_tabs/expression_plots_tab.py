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
    safe_rerun,
    log_streamlit_tab,
    log_streamlit_function,
    log_streamlit_event,
    group_datasets_by_organism,
    filter_genes_by_organism_datasets,
    get_organism_display_name
)
from ResultsIntegration import ResultsIntegrator

# Set up fragment decorator
setup_fragment_decorator()

logger = logging.getLogger(__name__)


@log_streamlit_tab("Expression Plots")
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
    st.header("Plot Gene Expression")
    st.markdown("**Violin plots showing gene expression distributions across sample groups.** Each panel represents one gene, with samples grouped by experimental conditions. Select datasets using the Dataset & Contrast Selection form in the sidebar.")

    # Display settings for expression plots
    display_settings = _render_display_settings()

    if not selected_datasets:
        log_streamlit_event("No datasets selected for expression plots")
        st.info("**Getting Started with Expression Plots:**")
        st.markdown("""
        1. **Select Datasets** in the "Dataset & Contrast Selection" section in the sidebar
        2. **Configure genes** using the "Heatmap Parameters" section (genes auto-selected based on your contrasts)
        3. **View expression plots** here once datasets and genes are selected

        **Alternative options:**
        - View individual dataset plots in the **Analyze Experiments** tab
        - Browse dataset information in the **View Dataset Info** tab
        """)
    elif not gene_sel:
        log_streamlit_event("No genes selected for expression plots")
        st.info("No genes available. Please select datasets and configure parameters in the sidebar to auto-select genes first.")
    else:
        # Group datasets by organism
        organism_groups = group_datasets_by_organism(ri, selected_datasets)

        if len(organism_groups) == 1:
            # Single organism - no sub-tabs needed
            organism = list(organism_groups.keys())[0]
            organism_datasets = organism_groups[organism]
            organism_genes = filter_genes_by_organism_datasets(ri, gene_sel, organism, organism_datasets)

            # Calculate pagination information for this organism
            total_pages, current_page, genes_per_page, current_genes = calculate_pagination_info(organism_genes)
            log_streamlit_event(f"Expression plots pagination: page {current_page}/{total_pages}, showing {len(current_genes)} genes for {organism}")

            # Add pagination controls at the top if needed
            if total_pages > 1:
                _render_top_pagination_controls(current_page, total_pages, "expression_single")

            # Create and display the expression plots
            _draw_expression_plots(
                ri,
                organism_genes,
                organism_datasets,
                hide_x_labels,
                current_page,
                total_pages,
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
                    organism_datasets = organism_groups[organism]
                    organism_genes = filter_genes_by_organism_datasets(ri, gene_sel, organism, organism_datasets)

                    if organism_genes:

                        # Calculate pagination information for this organism
                        total_pages, current_page, genes_per_page, current_genes = calculate_pagination_info(organism_genes)

                        # Add pagination controls at the top if needed
                        if total_pages > 1:
                            _render_top_pagination_controls(current_page, total_pages, f"expression_{organism}")

                        # Create and display the expression plots
                        _draw_expression_plots(
                            ri,
                            organism_genes,
                            organism_datasets,
                            hide_x_labels,
                            current_page,
                            total_pages,
                            **display_settings
                        )
                    else:
                        st.warning(f"No genes found for {get_organism_display_name(organism)} with current parameters.")
                        st.info("Try adjusting the significance thresholds in the sidebar or selecting more datasets for this species.")



@log_streamlit_function
def _render_display_settings() -> Dict[str, Any]:
    """Render display settings controls in the sidebar and return the settings."""
    with st.sidebar.expander("Display Settings", expanded=False):
        facet_font_size = st.slider("Facet title size", 8, 16, 10, key="expression_plots_font")
        lock_y_axis = st.checkbox("Lock y-axis across genes", value=False, key="expression_plots_lock_y")
        show_raw_points = st.checkbox("Show raw points", value=True, key="expression_plots_points")
        legend_position = st.selectbox("Legend position", ["Bottom", "Right", "Top"], index=0, key="expression_plots_legend")

    return {
        "facet_font_size": facet_font_size,
        "lock_y_axis": lock_y_axis,
        "show_raw_points": show_raw_points,
        "legend_position": legend_position,
    }


@log_streamlit_function
def _render_top_pagination_controls(current_page: int, total_pages: int, prefix: str = "expression"):
    """Render pagination controls at the top of the plots for convenience."""
    cols = st.columns([2, 1, 1, 1, 2])

    with cols[1]:
        prev_disabled = current_page <= 1
        if st.button("Previous", disabled=prev_disabled, key=f"{prefix}_prev"):
            st.session_state.page_num = max(1, current_page - 1)
            safe_rerun()

    with cols[2]:
        st.markdown(f"**Page {current_page}/{total_pages}**")

    with cols[3]:
        next_disabled = current_page >= total_pages
        if st.button("Next", disabled=next_disabled, key=f"{prefix}_next"):
            st.session_state.page_num = min(total_pages, current_page + 1)
            safe_rerun()


@st.fragment
@log_streamlit_function
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

            # Use cached figure creation for better performance
            fig2 = cached_figure_creation(
                "create_expression_plots",
                ri.results_dir,
                current_genes,
                "violin",
                dataset_selection,
                None,  # output_file
                hide_labels,
                page_num,
                facet_font_size,
                lock_y_axis,
                show_raw_points,
                legend_position,
                True,  # show_grid_lines
                0.3    # grid_opacity
            )

            if fig2:
                log_streamlit_event("Expression plots generated successfully")
                st.plotly_chart(fig2, use_container_width=True)
            else:
                log_streamlit_event("Failed to generate expression plots")
                st.error("Could not generate expression plots. Please check your selections.")

        except Exception as e:
            logger.error(f"Error generating expression plots: {str(e)}", exc_info=True)
            st.error(f"Error generating expression plots: {str(e)}")
            _display_expression_plot_error_details(e)


@log_streamlit_function
def _display_expression_plot_error_details(error: Exception):
    """Display detailed error information in an expandable section."""
    with st.expander("Expression Plot Error Details", expanded=False):
        st.code(traceback.format_exc())

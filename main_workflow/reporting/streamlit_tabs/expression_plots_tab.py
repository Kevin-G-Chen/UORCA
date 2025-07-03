"""
Expression Plots Tab for UORCA Explorer.

This tab displays violin plots showing gene expression distributions across sample groups.
Includes group selection and custom gene selection forms at the top.
"""

import logging
import streamlit as st
import pandas as pd
import traceback
from typing import List, Dict, Any, Optional

from .helpers import (
    check_ai_generating,
    setup_fragment_decorator,
    cached_figure_creation,
    calculate_pagination_info,
    safe_rerun,
    log_streamlit_tab,
    log_streamlit_function,
    log_streamlit_event,
    log_streamlit_user_action,
    group_datasets_by_organism,
    filter_genes_by_organism_datasets,
    get_organism_display_name
)
from ResultsIntegration import ResultsIntegrator

# Set up fragment decorator
setup_fragment_decorator()

logger = logging.getLogger(__name__)


@log_streamlit_tab("Expression Plots")
def render_expression_plots_tab(ri: ResultsIntegrator, selected_datasets: List[str], **kwargs):
    """
    Render the expression plots tab with group selection and gene selection forms.

    Args:
        ri: ResultsIntegrator instance
        selected_datasets: List of selected dataset IDs from sidebar
        **kwargs: Additional arguments (maintained for compatibility)
    """
    st.header("Plot Gene Expression")
    st.markdown("**Violin plots showing gene expression distributions across sample groups.** Configure your analysis parameters below, then view the expression plots.")

    # Get selected datasets from sidebar
    if not selected_datasets:
        selected_datasets = st.session_state.get('selected_datasets_from_sidebar', [])

    if not selected_datasets:
        st.info("**Getting Started with Expression Plots:**")
        st.markdown("""
        1. **Select Datasets** in the sidebar first
        2. **Choose Sample Groups** using the form below (populated from your datasets)
        3. **Enter Custom Genes** using the gene list form
        4. **View Results** - expression plots will be generated automatically

        **Note:** This tab shows expression levels across sample groups within datasets.
        """)
        return

    # Render group selection form
    selected_groups = _render_group_selection_form(ri, selected_datasets)

    # Render gene selection form
    selected_genes = _render_gene_selection_form(ri)

    # Display current selections summary
    if selected_groups and selected_genes:
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Datasets", len(selected_datasets))
        with col2:
            st.metric("Sample Groups", len(selected_groups))
        with col3:
            st.metric("Genes", len(selected_genes))

        st.info(f"🎯 **Custom Gene Selection**: Displaying {len(selected_genes)} genes across sample groups")

        # Group datasets by organism and render expression plots
        organism_groups = group_datasets_by_organism(ri, selected_datasets)

        if len(organism_groups) == 1:
            # Single organism - no sub-tabs needed
            organism = list(organism_groups.keys())[0]
            organism_datasets = organism_groups[organism]
            organism_genes = filter_genes_by_organism_datasets(ri, selected_genes, organism, organism_datasets)

            # Calculate pagination information for this organism
            total_pages, current_page, genes_per_page, current_genes = calculate_pagination_info(organism_genes, genes_per_page=20)  # 20 genes per page for 2x10 layout
            log_streamlit_event(f"Expression plots pagination: page {current_page}/{total_pages}, showing {len(current_genes)} genes for {organism}")

            # Add pagination controls at the top if needed
            if total_pages > 1:
                _render_top_pagination_controls(current_page, total_pages, "expression_single")

            # Create and display the expression plots
            _draw_expression_plots(
                ri,
                organism_genes,
                organism_datasets,
                selected_groups,
                current_page,
                total_pages
            )
        else:
            # Multiple organisms - create sub-tabs
            organism_names = list(organism_groups.keys())
            tab_names = [get_organism_display_name(org) for org in organism_names]
            tabs = st.tabs(tab_names)

            for i, organism in enumerate(organism_names):
                with tabs[i]:
                    organism_datasets = organism_groups[organism]
                    organism_genes = filter_genes_by_organism_datasets(ri, selected_genes, organism, organism_datasets)

                    if organism_genes:
                        # Calculate pagination information for this organism
                        total_pages, current_page, genes_per_page, current_genes = calculate_pagination_info(organism_genes, genes_per_page=20)

                        # Add pagination controls at the top if needed
                        if total_pages > 1:
                            _render_top_pagination_controls(current_page, total_pages, f"expression_{organism}")

                        # Create and display the expression plots
                        _draw_expression_plots(
                            ri,
                            organism_genes,
                            organism_datasets,
                            selected_groups,
                            current_page,
                            total_pages
                        )
                    else:
                        st.warning(f"No genes found for {get_organism_display_name(organism)} with current parameters.")
                        st.info("Try entering different genes or selecting more datasets for this species.")

    elif selected_groups and not selected_genes:
        st.warning("No genes selected. Please enter genes using the form above.")
    elif not selected_groups:
        st.info("Please select sample groups using the form above.")


@log_streamlit_function
def _render_group_selection_form(ri: ResultsIntegrator, selected_datasets: List[str]) -> List[str]:
    """Render the sample group selection form."""
    st.subheader("1. Select Sample Groups")

    # Initialize session state
    if 'selected_groups_expression' not in st.session_state:
        st.session_state['selected_groups_expression'] = set()

    with st.form("expression_group_selection"):
        # Get all unique groups from selected datasets
        all_groups = set()
        dataset_group_mapping = {}

        for analysis_id in selected_datasets:
            if analysis_id in ri.analysis_info:
                info = ri.analysis_info[analysis_id]
                accession = info.get("accession", analysis_id)
                unique_groups = info.get("unique_groups", [])

                dataset_group_mapping[accession] = unique_groups
                all_groups.update(unique_groups)

        if all_groups:
            # Create group data for display
            group_data = []
            for group in sorted(all_groups):
                # Find which datasets contain this group
                containing_datasets = []
                for dataset, groups in dataset_group_mapping.items():
                    if group in groups:
                        containing_datasets.append(dataset)

                group_data.append({
                    "Select": group in st.session_state['selected_groups_expression'],
                    "Group": group,
                    "Datasets": ", ".join(containing_datasets),
                    "Count": len(containing_datasets)
                })

            df = pd.DataFrame(group_data)

            st.info(f"Available sample groups from {len(selected_datasets)} selected datasets")

            edited_df = st.data_editor(
                df,
                use_container_width=True,
                hide_index=True,
                column_config={
                    "Select": st.column_config.CheckboxColumn(
                        "",
                        help="Check to include this sample group in expression plots",
                        default=False
                    ),
                    "Group": st.column_config.TextColumn(
                        "Sample Group",
                        help="Sample group name",
                        width=300
                    ),
                    "Datasets": st.column_config.TextColumn(
                        "Found in Datasets",
                        help="Datasets containing this group",
                        width=200
                    ),
                    "Count": st.column_config.NumberColumn(
                        "# Datasets",
                        help="Number of datasets containing this group",
                        width=100
                    )
                },
                key="expression_group_selection_table"
            )

            # Quick selection buttons
            col1, col2, col3, col4 = st.columns(4)
            with col1:
                select_all = st.form_submit_button("Select All", type="secondary")
            with col2:
                clear_all = st.form_submit_button("Clear All", type="secondary")
            with col3:
                apply_selection = st.form_submit_button("Apply Selection", type="primary")
            with col4:
                st.write("")  # Spacer

            if select_all:
                # Select all groups
                all_group_names = set([row["Group"] for _, row in df.iterrows()])
                st.session_state['selected_groups_expression'] = all_group_names
                log_streamlit_user_action(f"Selected all {len(all_group_names)} groups for expression plots")
                st.rerun()

            elif clear_all:
                # Clear all selections
                st.session_state['selected_groups_expression'] = set()
                log_streamlit_user_action("Cleared all group selections for expression plots")
                st.rerun()

            elif apply_selection:
                # Apply current selections
                selected_groups = set()
                if not edited_df.empty:
                    selected_rows = edited_df[edited_df["Select"] == True]
                    selected_groups = set(selected_rows["Group"].tolist())

                st.session_state['selected_groups_expression'] = selected_groups
                log_streamlit_user_action(f"Applied {len(selected_groups)} group selections for expression plots")
                st.rerun()

        else:
            st.info("No sample groups found in selected datasets")
            st.form_submit_button("Apply Selection", disabled=True)

    # Return current selections as list
    return list(st.session_state['selected_groups_expression'])


@log_streamlit_function
def _render_gene_selection_form(ri: ResultsIntegrator) -> List[str]:
    """Render the custom gene selection form."""
    st.subheader("2. Enter Custom Genes")

    # Initialize session state
    if 'custom_genes_expression' not in st.session_state:
        st.session_state['custom_genes_expression'] = []

    with st.form("expression_gene_selection"):
        st.markdown("**Enter genes to display in expression plots**")
        st.info("Expression plots show all available genes without significance filtering.")

        # Custom gene input
        custom_genes_input = st.text_area(
            "Gene List",
            height=150,
            placeholder="Enter one gene per line, e.g.:\nTP53\nEGFR\nMYC\nBRCA1",
            help="Enter gene symbols, one per line. All available genes will be shown regardless of significance.",
            key="expression_custom_genes",
            value="\n".join(st.session_state['custom_genes_expression'])
        )

        # Show preview for genes
        if custom_genes_input.strip():
            custom_genes_list = [gene.strip() for gene in custom_genes_input.strip().split('\n') if gene.strip()]
            if custom_genes_list:
                st.write(f"**Preview:** {len(custom_genes_list)} genes entered")
                preview_text = ", ".join(custom_genes_list[:10])
                if len(custom_genes_list) > 10:
                    preview_text += f", ... (+{len(custom_genes_list) - 10} more)"
                st.caption(preview_text)

        # Submit button
        submitted = st.form_submit_button("Apply Gene Selection", type="primary")

        if submitted:
            # Validate custom genes
            if not custom_genes_input.strip():
                st.error("Please enter at least one gene")
                return st.session_state['custom_genes_expression']

            custom_genes_list = [gene.strip() for gene in custom_genes_input.strip().split('\n') if gene.strip()]
            if not custom_genes_list:
                st.error("Please enter valid gene names")
                return st.session_state['custom_genes_expression']

            # Update session state
            st.session_state['custom_genes_expression'] = custom_genes_list
            log_streamlit_user_action(f"Expression plot genes updated: {len(custom_genes_list)} custom genes")
            st.rerun()

    # Show current gene validation
    if st.session_state['custom_genes_expression']:
        with st.expander("Gene Validation", expanded=False):
            genes = st.session_state['custom_genes_expression']
            st.write(f"**Total genes:** {len(genes)}")
            st.info("All entered genes will be displayed if found in the selected datasets, regardless of significance levels.")

    return st.session_state['custom_genes_expression']


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
    selected_groups: List[str],
    page_num: int,
    total_pgs: int
):
    """
    Create and display the expression plots using fragment isolation.
    """
    # Skip execution if AI is currently generating
    if check_ai_generating():
        return

    with st.spinner("Generating expression plots..."):
        try:
            # Calculate gene slice for the current page (20 genes per page for 2 plots per row)
            genes_per_page = 20
            start_idx = (page_num - 1) * genes_per_page
            end_idx = min(start_idx + genes_per_page, len(gene_selection))
            current_genes = gene_selection[start_idx:end_idx]

            # Fixed parameters for expression plots (2 plots per row)
            hide_x_labels = True
            facet_font_size = 10
            lock_y_axis = False
            show_raw_points = True
            legend_position = "Bottom"
            show_grid_lines = True
            grid_opacity = 0.3

            # Use the passed selected_groups parameter

            # Use cached figure creation for better performance
            fig2 = cached_figure_creation(
                "create_expression_plots",
                ri.results_dir,
                current_genes,
                "violin",  # plot_type
                dataset_selection,
                None,  # output_file
                hide_x_labels,
                page_num,
                facet_font_size,
                lock_y_axis,
                show_raw_points,
                legend_position,
                show_grid_lines,
                grid_opacity,
                selected_groups,  # Pass selected groups for filtering
                2  # plots_per_row - set to 2 for expression plots
            )

            if fig2:
                log_streamlit_event("Expression plots generated successfully")

                # Display information about the plots
                st.success(f"**Expression plots for {len(current_genes)} genes** - Page {page_num} of {total_pgs}")

                # Display the plot
                st.plotly_chart(fig2, use_container_width=True)

                # Add informational notes
                with st.expander("💡 About Expression Plots", expanded=False):
                    st.markdown("""
                    **How to interpret these plots:**
                    - Each panel shows one gene's expression across samples
                    - Samples are grouped by experimental conditions (sample groups)
                    - Violin plots show the distribution of expression values
                    - Points represent individual samples
                    - Values are in log₂CPM (counts per million, log-transformed)

                    **Sample grouping:**
                    - Groups are determined from the experimental design metadata
                    - Each group represents a different experimental condition
                    - Groups are shared across datasets when they have the same name

                    **Layout:**
                    - **2 plots per row** for optimal viewing
                    - Large gene sets are automatically split into pages
                    - Use pagination controls above to navigate between pages

                    **Data filtering:**
                    - **No significance filtering** - all genes are shown if present in data
                    - Unlike heatmaps, expression plots show raw expression levels
                    - Missing genes indicate they were not found in the selected datasets
                    """)
            else:
                log_streamlit_event("Failed to generate expression plots")
                st.error("Could not generate expression plots. Please check your selections.")

                # Provide troubleshooting information
                with st.expander("Troubleshooting", expanded=True):
                    st.markdown("""
                    **Common issues:**
                    - Selected genes not found in the datasets
                    - Sample grouping information missing
                    - Dataset compatibility issues

                    **Solutions:**
                    - Try selecting different genes or datasets
                    - Check that your datasets have completed analysis
                    - Verify that expression data is available
                    - Ensure sample groups are properly defined in the metadata
                    """)

        except Exception as e:
            logger.error(f"Error generating expression plots: {str(e)}", exc_info=True)
            st.error(f"Error generating expression plots: {str(e)}")
            _display_expression_plot_error_details(e)


@log_streamlit_function
def _display_expression_plot_error_details(error: Exception):
    """Display detailed error information in an expandable section."""
    with st.expander("Expression Plot Error Details", expanded=False):
        st.code(traceback.format_exc())
        st.markdown("""
        **Common issues:**
        - Selected genes not found in the datasets
        - Sample grouping information missing from analysis metadata
        - Dataset compatibility issues across different analyses
        - Expression data not available for selected genes
        """)

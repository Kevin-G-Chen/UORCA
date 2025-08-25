"""
Expression Plots Tab for UORCA Explorer.

This tab displays violin plots showing gene expression distributions across sample groups.
Simplified interface with combined form for group and gene selection.
"""

import os
import re
import logging
import streamlit as st
import pandas as pd
import traceback
from typing import List, Dict, Any, Optional

from .helpers import (
    check_ai_generating,
    setup_fragment_decorator,
    cached_figure_creation,
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
def render_expression_plots_tab(ri: ResultsIntegrator, selected_datasets: List[str]):
    """
    Render the expression plots tab with combined form for group and gene selection.

    Args:
        ri: ResultsIntegrator instance
        selected_datasets: List of selected dataset IDs from sidebar
    """
    st.header("Expression Plots")
    st.markdown("Visualise gene expression across sample groups.")

    # Get selected datasets from sidebar
    if not selected_datasets:
        selected_datasets = st.session_state.get('selected_datasets_from_sidebar', [])

    if not selected_datasets:
        st.markdown("To get started, select datasets from the sidebar on the left.")
        return

    # Render combined form for groups and genes
    form_results = _render_combined_form(ri, selected_datasets)

    if form_results:
        selected_groups = form_results['groups']
        selected_genes = form_results['genes']

        # Generate expression plots if we have both groups and genes
        if selected_groups and selected_genes:
            # Group datasets by organism
            organism_groups = group_datasets_by_organism(ri, selected_datasets)

            if len(organism_groups) == 1:
                # Single organism - no sub-tabs needed
                organism = list(organism_groups.keys())[0]
                organism_datasets = [ds for ds in selected_datasets if ds in organism_groups[organism]]
                organism_genes = filter_genes_by_organism_datasets(ri, selected_genes, organism, organism_datasets)

                if organism_genes:
                    log_streamlit_event(f"Expression plots: showing {len(organism_genes)} genes for {organism}")

                    # Create and display all expression plots at once
                    _draw_expression_plots(
                        ri,
                        organism_genes,
                        organism_datasets,
                        selected_groups
                    )
                else:
                    st.warning(f"No genes found for {get_organism_display_name(organism)} in selected datasets.")
            else:
                # Multiple organisms - create sub-tabs
                organism_names = list(organism_groups.keys())
                tab_names = [get_organism_display_name(org) for org in organism_names]
                tabs = st.tabs(tab_names)

                for i, organism in enumerate(organism_names):
                    with tabs[i]:
                        organism_datasets = [ds for ds in selected_datasets if ds in organism_groups[organism]]
                        organism_genes = filter_genes_by_organism_datasets(ri, selected_genes, organism, organism_datasets)

                        if organism_genes:
                            # Create and display all expression plots at once
                            _draw_expression_plots(
                                ri,
                                organism_genes,
                                organism_datasets,
                                selected_groups
                            )
                        else:
                            st.warning(f"No genes found for {get_organism_display_name(organism)} with current parameters.")


@log_streamlit_function
def _get_group_dataset_combinations(ri: ResultsIntegrator, selected_datasets: List[str]) -> List[Dict[str, Any]]:
    """Get all unique group-dataset combinations with duplicate handling."""
    group_dataset_combinations = []

    # Log a single message once per session to avoid verbosity
    global _logged_groups_once
    try:
        _logged_groups_once
    except NameError:
        _logged_groups_once = False
    if not _logged_groups_once:
        log_streamlit_event("Loading sample groups for all datasets")
        _logged_groups_once = True

    for dataset_id in selected_datasets:
        if dataset_id in ri.analysis_info:
            info = ri.analysis_info[dataset_id]
            accession = info.get("accession", dataset_id)
            unique_groups = info.get("unique_groups", [])

            # Do not log per-dataset to keep logs concise

            # Get dataset title
            dataset_title = ""
            if hasattr(ri, 'dataset_info') and dataset_id in ri.dataset_info:
                title = ri.dataset_info[dataset_id].get('title', '')
                if title and title.startswith('Title:'):
                    title = title[6:].strip()
                dataset_title = title

            for group in unique_groups:
                sample_count = _get_sample_count_for_group(ri, dataset_id, group)
                group_dataset_combinations.append({
                    "Select": True,  # Default to selected
                    "Sample Group": group,
                    "Dataset": accession,
                    "Title": dataset_title,
                    "Species": info.get("organism", "Unknown"),
                    "Samples": sample_count
                })

    # Check for duplicate group names and make them unique
    if group_dataset_combinations:
        # First pass: count occurrences of each group name
        group_counts = {}
        for item in group_dataset_combinations:
            original_group = item["Sample Group"]
            group_counts[original_group] = group_counts.get(original_group, 0) + 1

        # Second pass: rename items that have duplicates
        for item in group_dataset_combinations:
            original_group = item["Sample Group"]
            if group_counts[original_group] > 1:
                # Duplicate found - make unique by appending dataset
                unique_group = f"{original_group} ({item['Dataset']})"
                item["Sample Group"] = unique_group

    # Pre-select based on session state (similar to sidebar pattern)
    if 'selected_groups_from_expression' in st.session_state:
        selected_groups = st.session_state['selected_groups_from_expression']
        for item in group_dataset_combinations:
            group_id = f"{item['Sample Group']}_{item['Dataset']}"
            item['Select'] = group_id in selected_groups

    # Sort combinations by dataset GEO accession numeric portion, then by dataset and group
    def _acc_num(acc: str) -> int:
        m = re.search(r"(\d+)", str(acc) or "")
        return int(m.group(1)) if m else float('inf')

    group_dataset_combinations.sort(key=lambda x: (_acc_num(x["Dataset"]), x["Dataset"], x["Sample Group"]))

    return group_dataset_combinations


def _get_sample_count_for_group(ri: ResultsIntegrator, dataset_id: str, group_name: str) -> int:
    """Get the number of samples in a specific group for a dataset.

    Verbose logging is handled at the dataset level in _get_group_dataset_combinations
    to avoid logging per group.
    """
    try:
        # Look for the sample mapping file
        analysis_dir = os.path.join(ri.results_dir, dataset_id)
        sample_file_locations = [
            os.path.join(analysis_dir, "metadata", "edger_analysis_samples.csv"),
            os.path.join(analysis_dir, "edger_analysis_samples.csv"),
        ]

        edger_samples_file = None
        for path in sample_file_locations:
            if os.path.isfile(path):
                edger_samples_file = path
                break

        if not edger_samples_file:
            return 0

        # Load the sample mapping file
        sample_mapping_df = pd.read_csv(edger_samples_file, nrows=0)
        first_col = sample_mapping_df.columns[0]
        has_index = first_col == '' or first_col.startswith('Unnamed')
        sample_mapping_df = pd.read_csv(edger_samples_file, index_col=0 if has_index else None)

        # Get the group column from analysis_info
        group_col = None
        if ri.analysis_info and dataset_id in ri.analysis_info and 'merged_column' in ri.analysis_info[dataset_id]:
            group_col = ri.analysis_info[dataset_id]['merged_column']

        if group_col and group_col in sample_mapping_df.columns:
            # Count samples in this group
            count = (sample_mapping_df[group_col] == group_name).sum()
            return int(count)

        return 0
    except Exception as e:
        logger.warning(f"Error getting sample count for group {group_name} in {dataset_id}: {e}")
        return 0


@log_streamlit_function
def _render_combined_form(ri: ResultsIntegrator, selected_datasets: List[str]) -> Optional[Dict[str, Any]]:
    """Render the combined form for sample group and gene selection."""

    # Initialize session state for selected groups
    if 'selected_groups_from_expression' not in st.session_state:
        st.session_state['selected_groups_from_expression'] = set()

    # Check if datasets have changed and auto-select all groups by default
    current_datasets = set(selected_datasets)
    if 'previous_datasets_expression' not in st.session_state:
        st.session_state['previous_datasets_expression'] = set()

    if current_datasets != st.session_state['previous_datasets_expression']:
        # Datasets changed - select all groups by default
        temp_combinations = _get_group_dataset_combinations(ri, selected_datasets)
        all_groups = set()
        for item in temp_combinations:
            group_id = f"{item['Sample Group']}_{item['Dataset']}"
            all_groups.add(group_id)
        st.session_state['selected_groups_from_expression'] = all_groups
        st.session_state['previous_datasets_expression'] = current_datasets

    # Quick selection buttons (outside form)
    if group_dataset_combinations := _get_group_dataset_combinations(ri, selected_datasets):
        col1, col2 = st.columns(2)
        with col1:
            if st.button("Select All Groups", key="select_all_groups_expression"):
                # Select all group IDs (similar to sidebar pattern)
                all_groups = set()
                for item in group_dataset_combinations:
                    group_id = f"{item['Sample Group']}_{item['Dataset']}"
                    all_groups.add(group_id)
                st.session_state['selected_groups_from_expression'] = all_groups
                st.rerun()
        with col2:
            if st.button("Clear All Groups", key="clear_all_groups_expression"):
                st.session_state['selected_groups_from_expression'] = set()
                st.rerun()

    with st.form("expression_combined_form"):
        st.subheader("Configuration")

        # Sample Group Selection
        st.markdown("**Sample Groups**")

        # Get all unique group-dataset combinations
        group_dataset_combinations = _get_group_dataset_combinations(ri, selected_datasets)

        if group_dataset_combinations:
            df = pd.DataFrame(group_dataset_combinations)

            edited_df = st.data_editor(
                df,
                use_container_width=True,
                hide_index=True,
                column_config={
                    "Select": st.column_config.CheckboxColumn(
                        "",
                        help="Check to include this sample group",
                        default=True
                    ),
                    "Sample Group": st.column_config.TextColumn(
                        "Sample Group",
                        help="Sample group name",
                        width=160
                    ),
                    "Dataset": st.column_config.TextColumn(
                        "Dataset",
                        help="Dataset accession",
                        width=90
                    ),
                    "Species": st.column_config.TextColumn(
                        "Species",
                        help="Species/organism",
                        width=100
                    ),
                    "Samples": st.column_config.NumberColumn(
                        "Samples",
                        help="Number of samples in this group",
                        width=70
                    ),
                    "Title": st.column_config.TextColumn(
                        "Title",
                        help="Dataset title",
                        width=220
                    )
                },
                key="expression_combined_group_table"
            )
        else:
            edited_df = pd.DataFrame()

        # Gene Selection
        st.markdown("**Genes**")

        custom_genes_input = st.text_area(
            "Gene List",
            height=150,
            placeholder="Enter one gene per line, e.g.:\nTP53\nEGFR\nMYC\nBRCA1",
            help="Enter gene symbols, one per line",
            key="expression_combined_genes"
        )

        # Submit button
        submitted = st.form_submit_button("Generate Expression Plots", type="primary")

        if submitted:
            # Get selected groups from current form state
            selected_groups = []
            if not edited_df.empty:
                selected_rows = edited_df[edited_df["Select"] == True]
                # Create unique identifiers for each group-dataset combination
                for _, row in selected_rows.iterrows():
                    group_id = f"{row['Sample Group']}_{row['Dataset']}"
                    selected_groups.append(group_id)

                # Update session state with current selections
                st.session_state['selected_groups_from_expression'] = set(selected_groups)

            # Get genes from current form state
            selected_genes = []
            if custom_genes_input.strip():
                selected_genes = [gene.strip() for gene in custom_genes_input.strip().split('\n') if gene.strip()]

            # Validation
            if not selected_groups:
                st.error("Please select at least one sample group")
                return None

            if not selected_genes:
                st.error("Please enter at least one gene")
                return None

            # Show gene validation with detailed checking
            with st.expander("Gene Validation", expanded=False):
                st.write(f"**Total genes entered:** {len(selected_genes)}")

                # Check which genes are actually found in the selected datasets
                available_genes = set()
                for dataset_id in selected_datasets:
                    if dataset_id in ri.cpm_data:
                        cpm_df = ri.cpm_data[dataset_id]
                        if 'Gene' in cpm_df.columns:
                            available_genes.update(cpm_df['Gene'].tolist())

                found_genes = [gene for gene in selected_genes if gene in available_genes]
                missing_genes = [gene for gene in selected_genes if gene not in available_genes]

                if found_genes:
                    st.write(f"\n**{len(found_genes)} genes found and will be displayed**")

                if missing_genes:
                    st.write(f"\n**{len(missing_genes)} genes not found in ANY selected datasets:**")
                    if len(missing_genes) <= 10:
                        st.write(", ".join(missing_genes))
                    else:
                        st.write(f"{', '.join(missing_genes[:10])}, ... (+{len(missing_genes)-10} more)")

                if not found_genes:
                    st.write(f"\n**ERROR:** No genes found in selected datasets!")
                    return None

            log_streamlit_user_action(f"Expression plots: {len(selected_groups)} groups, {len(selected_genes)} genes")

            return {
                'groups': selected_groups,
                'genes': selected_genes
            }

    return None




@st.fragment
@log_streamlit_function
def _draw_expression_plots(
    ri: ResultsIntegrator,
    gene_selection: List[str],
    dataset_selection: List[str],
    selected_groups: List[str]
):
    """
    Create and display the expression plots using fragment isolation.
    """
    # Skip execution if AI is currently generating
    if check_ai_generating():
        return

    with st.spinner("Generating expression plots..."):
        try:
            # Use all genes at once (no pagination)
            current_genes = gene_selection

            # Fixed parameters for expression plots (2 plots per row)
            hide_x_labels = False  # Show x and y axis titles on each facet
            facet_font_size = 10
            lock_y_axis = False
            show_raw_points = True
            legend_position = "None"  # No legend
            show_grid_lines = True
            grid_opacity = 0.3

            # Convert group IDs back to actual group names for filtering
            # selected_groups contains items like "GroupName_DatasetAccession"
            actual_groups = []
            for group_id in selected_groups:
                # Extract just the group name part (before the underscore)
                if '_' in group_id:
                    group_name = group_id.rsplit('_', 1)[0]
                    actual_groups.append(group_name)
                else:
                    actual_groups.append(group_id)

            # Use cached figure creation for better performance
            fig2 = cached_figure_creation(
                "create_expression_plots",
                ri.results_dir,
                current_genes,
                "violin",  # plot_type
                dataset_selection,
                None,  # output_file
                hide_x_labels,
                1,  # page_num - always 1 since no pagination
                facet_font_size,
                lock_y_axis,
                show_raw_points,
                legend_position,
                show_grid_lines,
                grid_opacity,
                actual_groups,  # Pass actual group names for filtering
                2,  # plots_per_row - set to 2 for expression plots
                False  # show_box - remove boxplot from violin
            )

            if fig2:
                log_streamlit_event("Expression plots generated successfully")

                # Add shared y-axis title and formatting
                fig2.update_layout(
                    yaxis_title="log2(CPM)",
                    showlegend=False  # Remove legend
                )

                # Apply y-axis formatting to all subplots
                fig2.update_yaxes(
                    title_text="log2(CPM)",
                    tickformat=".1f"  # Limit y-axis to 1 decimal place
                )

                # Display the plot
                st.plotly_chart(fig2, use_container_width=True)

            else:
                log_streamlit_event("Failed to generate expression plots")
                st.error("Could not generate expression plots. Please check your selections.")

                # Provide troubleshooting information
                with st.expander("Troubleshooting", expanded=False):
                    st.markdown("""
                    **Common issues:**
                    - Selected genes not found in the datasets
                    - Sample grouping information missing
                    - No samples matching selected groups

                    **Solutions:**
                    - Try selecting different genes or datasets
                    - Verify that expression data is available
                    """)

        except Exception as e:
            logger.error(f"Error generating expression plots: {str(e)}", exc_info=True)
            st.error(f"Error generating expression plots: {str(e)}")
            _display_expression_plot_error_details(e)


@log_streamlit_function
def _display_expression_plot_error_details(error: Exception):
    """Display detailed error information in an expandable section."""
    with st.expander("Error Details", expanded=False):
        st.code(traceback.format_exc())

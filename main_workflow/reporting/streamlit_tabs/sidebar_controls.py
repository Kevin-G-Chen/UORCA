"""
Sidebar Controls for UORCA Explorer - Form-Based Version.

This module handles all sidebar controls using Streamlit forms with expandable sections
for heatmap and expression plots configuration.
"""

import os
import logging
import streamlit as st
import pandas as pd
from typing import List, Dict, Any, Tuple, Optional

from .helpers import (
    cached_identify_important_genes,
    get_integrator,
    cached_get_all_genes_from_integrator,
    load_environment,
    log_streamlit_function,
    log_streamlit_event,
    log_streamlit_user_action
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
    st.sidebar.title("🧬 UORCA Explorer")

    # Results directory input (outside of forms)
    _render_results_directory_input(results_dir)

    # Initialize return parameters with defaults
    params = {
        'gene_sel': [],
        'selected_contrasts': [],
        'heatmap_params': {
            'lfc_thresh': 1.0,
            'pvalue_thresh': 0.05,
            'contrasts': []
        },
        'expression_params': {
            'datasets': []
        }
    }

    # Heatmap Configuration Form
    heatmap_params = _render_heatmap_form(ri, results_dir)
    if heatmap_params:
        params['heatmap_params'] = heatmap_params
        params['selected_contrasts'] = heatmap_params['contrasts']

    # Expression Plots Configuration Form (placeholder for now)
    expression_params = _render_expression_form(ri)
    if expression_params:
        params['expression_params'] = expression_params

    # Auto-select genes based on heatmap parameters
    if heatmap_params and heatmap_params['contrasts']:
        gene_sel = _auto_select_genes(ri, results_dir, heatmap_params)
        params['gene_sel'] = gene_sel
        st.session_state['current_gene_selection'] = gene_sel

    # Legacy parameter mapping for compatibility
    params.update({
        'pvalue_thresh': heatmap_params['pvalue_thresh'] if heatmap_params else 0.05,
        'lfc_thresh': heatmap_params['lfc_thresh'] if heatmap_params else 1.0,
        'effective_pvalue_thresh': heatmap_params['pvalue_thresh'] if heatmap_params else 0.05,
        'effective_lfc_thresh': heatmap_params['lfc_thresh'] if heatmap_params else 1.0,
        'use_dynamic_filtering': True,
        'hide_empty_rows_cols': True,
        'hide_x_labels': True,
        'show_advanced': False
    })

    # Configuration status display
    _render_configuration_status(params)

    # Help section
    _render_help_section()

    return params


@log_streamlit_function
def _render_results_directory_input(current_results_dir: str):
    """Render the results directory input field."""
    # Get the default results directory
    default_dir = os.getenv("UORCA_DEFAULT_RESULTS_DIR")
    if not default_dir:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        default_dir = os.path.join(os.path.dirname(os.path.dirname(script_dir)), "UORCA_results")
        if not os.path.exists(default_dir):
            default_dir = os.path.dirname(os.path.dirname(script_dir))

    # Results directory input
    results_dir = st.sidebar.text_input(
        "Results directory",
        value="/UORCA_results" if os.path.exists('/UORCA_results') else current_results_dir or default_dir,
        help="Path to UORCA results. If running in container, use '/UORCA_results'. Otherwise, use the full path to your results directory."
    )


@log_streamlit_function
def _render_heatmap_form(ri: ResultsIntegrator, results_dir: str) -> Optional[Dict[str, Any]]:
    """Render the heatmap configuration form."""
    with st.sidebar.expander("🌡️ Heatmap Configuration", expanded=True):
        with st.form("heatmap_config"):
            st.subheader("Parameters")

            # Parameter inputs
            col1, col2 = st.columns(2)
            with col1:
                lfc_thresh = st.text_input(
                    "Log2FC Threshold",
                    value="1.0",
                    help="Absolute log2 fold change threshold"
                )
            with col2:
                pvalue_thresh = st.text_input(
                    "P-value Threshold",
                    value="0.05",
                    help="Adjusted p-value threshold"
                )

            # Validate parameters
            try:
                lfc_val = float(lfc_thresh)
                pval_val = float(pvalue_thresh)
            except ValueError:
                st.error("Please enter valid numeric values for thresholds")
                lfc_val, pval_val = 1.0, 0.05

            st.subheader("Contrast Selection")

            # Create contrast selection table
            contrast_data = _create_contrast_table_data(ri)

            if contrast_data:
                # Convert to DataFrame for easier handling
                df = pd.DataFrame(contrast_data)

                # Use data_editor for interactive selection
                edited_df = st.data_editor(
                    df,
                    use_container_width=True,
                    hide_index=True,
                    column_config={
                        "Select": st.column_config.CheckboxColumn(
                            "Select",
                            help="Check to include this contrast in the heatmap",
                            default=False
                        ),
                        "Dataset": st.column_config.TextColumn(
                            "Dataset",
                            help="Dataset identifier"
                        ),
                        "Contrast": st.column_config.TextColumn(
                            "Contrast",
                            help="Contrast identifier"
                        ),
                        "Description": st.column_config.TextColumn(
                            "Description",
                            help="Contrast description",
                            width="large"
                        )
                    },
                    key="heatmap_contrast_selection"
                )

                # Submit button
                submitted = st.form_submit_button("Update Heatmap", type="primary")

                if submitted:
                    # Get selected contrasts
                    selected_contrasts = []
                    if not edited_df.empty:
                        selected_rows = edited_df[edited_df["Select"] == True]
                        selected_contrasts = [
                            (row["Dataset"], row["Contrast"])
                            for _, row in selected_rows.iterrows()
                        ]

                    log_streamlit_user_action(f"Heatmap form submitted: {len(selected_contrasts)} contrasts selected")

                    return {
                        'lfc_thresh': lfc_val,
                        'pvalue_thresh': pval_val,
                        'contrasts': selected_contrasts
                    }
            else:
                st.info("No contrasts available")
                st.form_submit_button("Update Heatmap", disabled=True)

    return None


@log_streamlit_function
def _render_expression_form(ri: ResultsIntegrator) -> Optional[Dict[str, Any]]:
    """Render the expression plots configuration form (placeholder)."""
    with st.sidebar.expander("📈 Expression Plots Configuration", expanded=False):
        with st.form("expression_config"):
            st.info("Expression plot configuration will be implemented in future versions.")
            st.form_submit_button("Update Expression Plots", disabled=True)

    return None


@log_streamlit_function
def _create_contrast_table_data(ri: ResultsIntegrator) -> List[Dict[str, Any]]:
    """Create data for the contrast selection table."""
    contrast_data = []

    for analysis_id, info in ri.analysis_info.items():
        for contrast in info.get("contrasts", []):
            contrast_id = contrast["name"]
            description = contrast.get("description", "")

            # Truncate long descriptions
            if len(description) > 100:
                description = description[:97] + "..."

            contrast_data.append({
                "Select": False,  # Default to unselected
                "Dataset": analysis_id,
                "Contrast": contrast_id,
                "Description": description
            })

    return contrast_data


@log_streamlit_function
def _auto_select_genes(
    ri: ResultsIntegrator,
    results_dir: str,
    heatmap_params: Dict[str, Any]
) -> List[str]:
    """Auto-select genes based on heatmap parameters."""
    if not heatmap_params.get('contrasts'):
        return []

    try:
        # Use cached function to get important genes
        top_genes = cached_identify_important_genes(
            results_dir=results_dir,
            top_frequent=20,  # Default to 20 frequent genes
            top_unique=0,     # No unique genes for now
            max_contrasts_for_unique=0,
            min_unique_per_contrast=1,
            p_value_threshold=heatmap_params['pvalue_thresh'],
            lfc_threshold=heatmap_params['lfc_thresh']
        )

        # Limit to 50 genes for performance
        limited_genes = top_genes[:50] if len(top_genes) > 50 else top_genes

        if len(top_genes) > 50:
            st.sidebar.info(f"Auto-selected top 50 of {len(top_genes)} important genes")
        else:
            st.sidebar.success(f"Auto-selected {len(limited_genes)} important genes")

        log_streamlit_event(f"Auto-selected {len(limited_genes)} genes")
        return limited_genes

    except Exception as e:
        logger.error(f"Error in auto gene selection: {e}")
        st.sidebar.error(f"Error selecting genes: {str(e)}")
        return []


@log_streamlit_function
def _render_configuration_status(params: Dict[str, Any]):
    """Display current configuration status in the sidebar."""
    st.sidebar.divider()
    st.sidebar.subheader("📊 Current Configuration")

    # Show heatmap configuration
    heatmap_params = params.get('heatmap_params', {})
    gene_count = len(params.get('gene_sel', []))
    contrast_count = len(params.get('selected_contrasts', []))

    if heatmap_params:
        st.sidebar.success("✅ Heatmap Configured")
        st.sidebar.caption(
            f"🧬 {gene_count} genes | "
            f"🔍 {contrast_count} contrasts | "
            f"📊 LFC≥{heatmap_params.get('lfc_thresh', 'N/A')} | "
            f"P≤{heatmap_params.get('pvalue_thresh', 'N/A')}"
        )
    else:
        st.sidebar.info("ℹ️ Configure heatmap parameters above")

    # Show expression configuration status
    expression_params = params.get('expression_params', {})
    if expression_params and expression_params.get('datasets'):
        st.sidebar.success("✅ Expression Plots Configured")
        dataset_count = len(expression_params['datasets'])
        st.sidebar.caption(f"📈 {dataset_count} datasets selected")
    else:
        st.sidebar.info("📈 Expression plots: Not yet configured")


@log_streamlit_function
def _render_help_section():
    """Render the help section at the bottom of the sidebar."""
    st.sidebar.divider()
    with st.sidebar.expander("ℹ️ Help", expanded=False):
        st.markdown(
            """
            ### How to Use Forms

            1. **Configure Parameters**: Set your significance thresholds
            2. **Select Contrasts**: Check the contrasts you want to analyze
            3. **Submit Form**: Click "Update Heatmap" to apply changes
            4. **View Results**: Switch to the heatmap tab to see results

            ### Form-Based Interface
            - Changes only apply when you submit the form
            - This prevents constant recomputation as you adjust settings
            - Submit the form after making all your selections

            ### Auto Gene Selection
            - Genes are automatically selected based on your parameters
            - Uses frequency across selected contrasts
            - Limited to top 50 genes for performance

            ---
            *Powered by Streamlit & Plotly*
            """
        )

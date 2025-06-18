"""
Sidebar Controls for UORCA Explorer.

This module handles all sidebar controls including gene selection, parameters,
and visualization options.
"""

import os
import logging
import streamlit as st
import tempfile
import shutil
import pandas as pd
from datetime import datetime
from pathlib import Path
from typing import List, Dict, Any, Tuple, Optional

from .helpers import (
    cached_identify_important_genes,
    get_all_genes_from_integrator,
    get_integrator,
    cached_get_all_genes_from_integrator,
    safe_rerun,
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
    Render all sidebar controls and return the selected parameters.

    Args:
        ri: ResultsIntegrator instance
        results_dir: Path to results directory

    Returns:
        Dictionary containing all selected parameters and gene selections
    """
    st.sidebar.title("🧬 UORCA Explorer")

    # Results directory input
    _render_results_directory_input(results_dir)

    # Get all genes for selection (cached)
    all_genes = cached_get_all_genes_from_integrator(results_dir)

    # Advanced options toggle
    show_advanced = st.sidebar.checkbox("Show advanced options", value=False)

    # Parameter controls
    params = _render_parameter_controls(show_advanced)

    # Gene selection
    gene_sel = _render_gene_selection(ri, results_dir, all_genes, params, show_advanced)

    # Visualization options
    viz_options = _render_visualization_options()

    # Pagination controls
    pagination_info = _render_pagination_controls(gene_sel)

    # Export options
    _render_export_options(ri, params, gene_sel, viz_options)

    # Help section
    _render_help_section()

    # Store current gene selection in session state for other tabs
    st.session_state['current_gene_selection'] = gene_sel

    # Return all parameters
    return {
        'gene_sel': gene_sel,
        'show_advanced': show_advanced,
        'pagination_info': pagination_info,
        'viz_options': viz_options,
        **params
    }





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
def _render_parameter_controls(show_advanced: bool) -> Dict[str, Any]:
    """Render parameter control widgets and return selected values."""
    st.sidebar.subheader("Advanced Options")

    # Initialize default values
    params = {
        'pvalue_thresh': 0.05,
        'lfc_thresh': 1.0,
        'min_unique': 0,
        'max_contrasts_unique': 0,
        'top_frequent_genes': 20,
        'top_unique_genes': 0,
        'hide_x_labels': True,
        'use_separate_heatmap_filters': False,
        'heatmap_pvalue_thresh': 0.05,
        'heatmap_lfc_thresh': 1.0,
        'use_dynamic_filtering': True,
        'effective_pvalue_thresh': 0.05,
        'effective_lfc_thresh': 1.0,
        'hide_empty_rows_cols': True
    }

    if show_advanced:
        params.update(_render_advanced_parameters())
    else:
        params.update(_render_simple_parameters())

    return params


@log_streamlit_function
def _render_advanced_parameters() -> Dict[str, Any]:
    """Render advanced parameter controls."""
    params = {}

    # P-value threshold
    col1, col2 = st.sidebar.columns([3, 1])
    with col1:
        params['pvalue_thresh'] = st.slider(
            "Adjusted P-value threshold", 0.0, 0.1, 0.05, 0.001,
            help="Genes must have adj.P.Val below this threshold to be considered significant"
        )
    with col2:
        pvalue_text = st.text_input(
            "P-value threshold",
            value=f"{params['pvalue_thresh']:.3f}",
            key="pvalue_text",
            label_visibility="collapsed"
        )
        try:
            params['pvalue_thresh'] = float(pvalue_text)
        except ValueError:
            pass

    # Log fold change threshold
    col1, col2 = st.sidebar.columns([3, 1])
    with col1:
        params['lfc_thresh'] = st.slider(
            "Log2FC threshold", 0.0, 5.0, 1.0, 0.1,
            help="Absolute log2 fold change threshold"
        )
    with col2:
        lfc_text = st.text_input(
            "LFC threshold",
            value=f"{params['lfc_thresh']:.1f}",
            key="lfc_text",
            label_visibility="collapsed"
        )
        try:
            params['lfc_thresh'] = float(lfc_text)
        except ValueError:
            pass

    # Auto-selection parameters
    st.sidebar.subheader("Auto-selection Parameters")
    st.sidebar.markdown("**Auto-selection strategy:**")
    st.sidebar.markdown("**Frequent DEGs**: Genes significant across multiple contrasts")

    # Frequent DEGs
    col1, col2 = st.sidebar.columns([3, 1])
    with col1:
        params['top_frequent_genes'] = st.slider(
            "Frequent DEGs to include", 5, 50, 20,
            help="Number of genes that are differentially expressed across the most contrasts"
        )
    with col2:
        freq_text = st.text_input(
            "Frequent genes count",
            value=str(params['top_frequent_genes']),
            key="freq_text",
            label_visibility="collapsed"
        )
        try:
            params['top_frequent_genes'] = int(freq_text)
        except ValueError:
            pass



    # Visualization options
    st.sidebar.subheader("Visualization Options")
    params['hide_x_labels'] = st.sidebar.checkbox(
        "Hide x-axis labels in expression plots", value=True
    )

    # Heatmap options
    st.sidebar.subheader("Heatmap Options")
    params['use_separate_heatmap_filters'] = st.sidebar.checkbox(
        "Use different significance filters for heatmap", value=False,
        help="By default, heatmap uses same filters as Auto-selected DEGs. Enable this to set different thresholds."
    )

    if params['use_separate_heatmap_filters']:
        params.update(_render_separate_heatmap_filters())
        params['use_dynamic_filtering'] = True
        params['effective_pvalue_thresh'] = params['heatmap_pvalue_thresh']
        params['effective_lfc_thresh'] = params['heatmap_lfc_thresh']
    else:
        params['use_dynamic_filtering'] = True
        params['effective_pvalue_thresh'] = params['pvalue_thresh']
        params['effective_lfc_thresh'] = params['lfc_thresh']
        st.sidebar.info("🎯 Heatmap uses same significance filters as Auto-selected DEGs")

    params['hide_empty_rows_cols'] = True

    return params


@log_streamlit_function
def _render_separate_heatmap_filters() -> Dict[str, Any]:
    """Render separate heatmap filter controls."""
    st.sidebar.markdown("**Heatmap-specific filters:**")

    params = {}

    # Heatmap P-value threshold
    col1, col2 = st.sidebar.columns([3, 1])
    with col1:
        params['heatmap_pvalue_thresh'] = st.slider(
            "Heatmap P-value threshold", 0.0, 0.1, 0.05, 0.001,
            help="P-value threshold for heatmap coloring"
        )
    with col2:
        heatmap_pvalue_text = st.text_input(
            "Heatmap P-value",
            value=f"{params['heatmap_pvalue_thresh']:.3f}",
            key="heatmap_pvalue_text",
            label_visibility="collapsed"
        )
        try:
            params['heatmap_pvalue_thresh'] = float(heatmap_pvalue_text)
        except ValueError:
            pass

    # Heatmap LFC threshold
    col1, col2 = st.sidebar.columns([3, 1])
    with col1:
        params['heatmap_lfc_thresh'] = st.slider(
            "Heatmap LFC threshold", 0.0, 5.0, 1.0, 0.1,
            help="LFC threshold for heatmap coloring"
        )
    with col2:
        heatmap_lfc_text = st.text_input(
            "Heatmap LFC",
            value=f"{params['heatmap_lfc_thresh']:.1f}",
            key="heatmap_lfc_text",
            label_visibility="collapsed"
        )
        try:
            params['heatmap_lfc_thresh'] = float(heatmap_lfc_text)
        except ValueError:
            pass

    return params


@log_streamlit_function
def _render_simple_parameters() -> Dict[str, Any]:
    """Render simple parameter controls."""
    params = {}

    # P-value threshold
    col1, col2 = st.sidebar.columns([3, 1])
    with col1:
        params['pvalue_thresh'] = st.slider(
            "Adjusted P-value threshold", 0.0, 0.1, 0.05, 0.001,
            help="Use adjusted P-value (adj.P.Val) for filtering DEGs"
        )
    with col2:
        pvalue_text = st.text_input(
            "P-value threshold",
            value=f"{params['pvalue_thresh']:.3f}",
            key="pvalue_text_simple",
            label_visibility="collapsed"
        )
        try:
            params['pvalue_thresh'] = float(pvalue_text)
        except ValueError:
            pass

    # LFC threshold
    col1, col2 = st.sidebar.columns([3, 1])
    with col1:
        params['lfc_thresh'] = st.slider(
            "abs(Log2FC) threshold", 0.0, 5.0, 1.0, 0.1,
            help="Genes must have absolute log2 fold change above this threshold to be considered significant"
        )
    with col2:
        lfc_text = st.text_input(
            "LFC threshold",
            value=f"{params['lfc_thresh']:.1f}",
            key="lfc_text_simple",
            label_visibility="collapsed"
        )
        try:
            params['lfc_thresh'] = float(lfc_text)
        except ValueError:
            pass

    # Use same filters for heatmap
    params['use_dynamic_filtering'] = True
    params['effective_pvalue_thresh'] = params['pvalue_thresh']
    params['effective_lfc_thresh'] = params['lfc_thresh']
    params['hide_empty_rows_cols'] = True

    st.sidebar.info("🎯 Heatmap uses same significance filters as Auto-selected DEGs")

    # Auto-selection parameters (simplified)
    st.sidebar.markdown("**Auto-selection strategy:**")
    st.sidebar.markdown("**Frequent DEGs**: Genes significant across multiple contrasts")

    # Frequent DEGs
    col1, col2 = st.sidebar.columns([3, 1])
    with col1:
        params['top_frequent_genes'] = st.slider(
            "Frequent DEGs to include", 5, 50, 20,
            help="Number of genes that are differentially expressed across the most contrasts"
        )
    with col2:
        freq_text = st.text_input(
            "Frequent genes count",
            value=str(params['top_frequent_genes']),
            key="freq_text_simple",
            label_visibility="collapsed"
        )
        try:
            params['top_frequent_genes'] = int(freq_text)
        except ValueError:
            pass

    # Set defaults for simple mode (removed contrast-specific controls)
    params['top_unique_genes'] = 0
    params['max_contrasts_unique'] = 0
    params['min_unique'] = 0
    params['hide_x_labels'] = True

    return params


@log_streamlit_function
def _render_gene_selection(
    ri: ResultsIntegrator,
    results_dir: str,
    all_genes: List[str],
    params: Dict[str, Any],
    show_advanced: bool
) -> List[str]:
    """Render gene selection controls and return selected genes."""
    st.sidebar.subheader("Gene Selection")
    gene_select_method = st.sidebar.radio(
        "Selection method:",
        ["Auto-selected DEGs", "Custom"],
        index=0  # Default to "Auto-selected DEGs"
    )

    if gene_select_method == "Custom":
        log_streamlit_user_action("Selected custom gene input method")
        return _render_custom_gene_selection(all_genes)
    else:
        log_streamlit_user_action("Selected auto gene selection method")
        return _render_auto_gene_selection(ri, results_dir, params, show_advanced)


@log_streamlit_function
def _render_custom_gene_selection(all_genes: List[str]) -> List[str]:
    """Render custom gene selection interface."""
    st.sidebar.write("Enter genes (one per line):")
    gene_input = st.sidebar.text_area(
        "Gene list:",
        height=150,
        placeholder="MYCN\nALK\nPHOX2B\nATRX\nTP53\nTERT\nARID1A\nARID1B\nNF1\nBARD1",
        label_visibility="collapsed"
    )

    if gene_input.strip():
        # Parse genes from textarea (one per line)
        input_genes = [gene.strip() for gene in gene_input.strip().split('\n') if gene.strip()]
        # Filter to only genes that exist in the data
        gene_sel = [gene for gene in input_genes if gene in all_genes]

        # Show validation info
        if len(input_genes) != len(gene_sel):
            missing_genes = [gene for gene in input_genes if gene not in all_genes]
            st.sidebar.warning(
                f"⚠️ {len(missing_genes)} genes not found in data: "
                f"{', '.join(missing_genes[:5])}{' ...' if len(missing_genes) > 5 else ''}"
            )

        if gene_sel:
            log_streamlit_user_action(f"Custom gene selection: {len(gene_sel)} valid genes")
            st.sidebar.success(f"✅ {len(gene_sel)} genes selected")
            st.sidebar.caption(f"🧬 {len(gene_sel)} genes selected")
        else:
            log_streamlit_event("Custom gene selection: no valid genes found")
            st.sidebar.error("❌ No valid genes found")

        return gene_sel
    else:
        return []


@log_streamlit_function
def _render_auto_gene_selection(
    ri: ResultsIntegrator,
    results_dir: str,
    params: Dict[str, Any],
    show_advanced: bool
) -> List[str]:
    """Render auto gene selection interface with caching."""
    # Use cached function to get important genes
    top_genes = cached_identify_important_genes(
        results_dir=results_dir,
        top_frequent=params['top_frequent_genes'],
        top_unique=params['top_unique_genes'],
        max_contrasts_for_unique=params['max_contrasts_unique'],
        min_unique_per_contrast=params['min_unique'] if show_advanced else 1,
        p_value_threshold=params['pvalue_thresh'],
        lfc_threshold=params['lfc_thresh']
    )

    # Handle gene limiting and ranking
    limited_genes = _limit_and_rank_genes(ri, top_genes)

    # Render selection interface
    return _render_gene_multiselect(limited_genes, len(top_genes))


@st.cache_data(show_spinner=False)
@log_streamlit_function
def _get_max_lfc_per_gene(results_dir: str, top_genes_tuple: tuple) -> Dict[str, float]:
    """Cached vectorized computation of max absolute LFC per gene."""
    # Convert tuple back to set for faster lookup
    top_genes_set = set(top_genes_tuple)

    # Load the integrator to get DEG data
    ri, _ = get_integrator(results_dir)
    if not ri:
        return {}

    # Collect all DEG tables into one DataFrame
    pieces = []
    for contrasts in ri.deg_data.values():
        for df in contrasts.values():
            if {'Gene', 'logFC'}.issubset(df.columns):
                pieces.append(df[['Gene', 'logFC']])

    if not pieces:
        return {}

    # Concatenate all DEG data
    big_df = pd.concat(pieces, ignore_index=True)

    # Filter to only genes in top_genes for efficiency
    big_df = big_df[big_df['Gene'].isin(top_genes_set)]

    if big_df.empty:
        return {}

    # Compute max absolute logFC per gene using vectorized operations
    max_lfc = big_df.assign(absFC=big_df['logFC'].abs()) \
                   .groupby('Gene')['absFC'] \
                   .max()

    return max_lfc.to_dict()

@log_streamlit_function
def _limit_and_rank_genes(ri: ResultsIntegrator, top_genes: List[str]) -> List[str]:
    """Limit genes to 200 maximum, selecting by highest LFC if needed."""
    if len(top_genes) <= 200:
        return top_genes

    # Get cached max LFC data using vectorized operations
    # Convert to tuple for caching (lists aren't hashable)
    gene_lfc_map = _get_max_lfc_per_gene(ri.results_dir, tuple(top_genes))

    # Sort by highest LFC and take top 200
    sorted_genes = sorted(top_genes, key=lambda g: gene_lfc_map.get(g, 0), reverse=True)
    limited_genes = sorted_genes[:200]

    st.sidebar.warning(
        f"⚠️ Showing top 200 of {len(top_genes)} important genes (ranked by highest LFC). "
        f"{len(top_genes)-200} genes excluded."
    )

    return limited_genes


@log_streamlit_function
def _render_gene_multiselect(limited_genes: List[str], total_genes: int) -> List[str]:
    """Render the gene multiselect widget with select all button."""
    col1, col2 = st.sidebar.columns([3, 1])
    with col2:
        st.markdown("")  # Add some spacing
        if st.button("Select All", key="select_all_genes_important") and limited_genes:
            log_streamlit_user_action(f"Selected all {len(limited_genes)} auto-selected genes")
            st.session_state['gene_sel_important'] = limited_genes
            safe_rerun()

    # Use session state if "Select All" was clicked or use default selection
    default_genes = st.session_state.get('gene_sel_important', limited_genes if limited_genes else [])
    gene_sel = st.sidebar.multiselect(
        "Auto-selected DEGs:",
        limited_genes if limited_genes else [],
        format_func=lambda x: x,
        default=default_genes,
        key="gene_sel_important",
    )

    # Add count badge
    if gene_sel:
        st.sidebar.caption(f"🧬 {len(gene_sel)} genes selected")

    # Clear the session state after use
    if 'gene_sel_important' in st.session_state:
        del st.session_state['gene_sel_important']

    return gene_sel


@log_streamlit_function
def _render_visualization_options() -> Dict[str, Any]:
    """Render visualization options and return selected values."""
    st.sidebar.subheader("Plot Options")
    return {
        'plot_type': "violin",  # Always use violin plots
        'hide_x_labels': True  # This is set in parameters section
    }


@log_streamlit_function
def _render_pagination_controls(gene_sel: List[str]) -> Dict[str, Any]:
    """Render pagination controls if needed."""
    genes_per_page = 30
    total_pages = (len(gene_sel) + genes_per_page - 1) // genes_per_page if gene_sel else 1

    if total_pages > 1:
        st.sidebar.markdown("---")
        st.sidebar.subheader("Gene Pagination")
        st.sidebar.info(f"Your selection contains {len(gene_sel)} genes, showing {genes_per_page} per page.")

        # Create three columns for pagination controls
        col1, col2, col3 = st.sidebar.columns([1, 2, 1])

        with col1:
            prev_disabled = st.session_state.get('page_num', 1) <= 1
            if st.button("◀", disabled=prev_disabled, key="prev_page"):
                log_streamlit_user_action("Navigated to previous page")
                st.session_state.page_num = max(1, st.session_state.get('page_num', 1) - 1)
                safe_rerun()

        with col2:
            # Initialize page_num in session state if not present
            if 'page_num' not in st.session_state:
                st.session_state.page_num = 1

            page_num = st.number_input(
                "Page",
                min_value=1,
                max_value=total_pages,
                value=st.session_state.page_num,
                step=1,
                key="page_input"
            )
            st.session_state.page_num = page_num

        with col3:
            next_disabled = st.session_state.get('page_num', 1) >= total_pages
            if st.button("▶", disabled=next_disabled, key="next_page"):
                log_streamlit_user_action("Navigated to next page")
                st.session_state.page_num = min(total_pages, st.session_state.get('page_num', 1) + 1)
                safe_rerun()

        # Show page indicator
        st.sidebar.caption(f"Page {st.session_state.page_num} of {total_pages}")
    else:
        if 'page_num' not in st.session_state:
            st.session_state.page_num = 1

    return {
        'total_pages': total_pages,
        'current_page': st.session_state.get('page_num', 1),
        'genes_per_page': genes_per_page
    }


@log_streamlit_function
def _render_export_options(
    ri: ResultsIntegrator,
    params: Dict[str, Any],
    gene_sel: List[str],
    viz_options: Dict[str, Any]
):
    """Render export options section."""
    st.sidebar.divider()
    st.sidebar.subheader("Export Options")
    export_format = st.sidebar.selectbox("Export format:", ["HTML", "CSV"])

    if st.sidebar.button("Export Current View"):
        log_streamlit_user_action(f"Started export in {export_format} format")
        if export_format == "HTML":
            _export_html(ri, params, gene_sel, viz_options)
        else:
            _export_csv(ri, gene_sel)


@log_streamlit_function
def _export_html(
    ri: ResultsIntegrator,
    params: Dict[str, Any],
    gene_sel: List[str],
    viz_options: Dict[str, Any]
):
    """Export current view as HTML."""
    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            # Generate the interactive report
            output_dir = ri.create_integrated_report(
                top_frequent=params['top_frequent_genes'],
                top_unique=params['top_unique_genes'],
                plot_type=viz_options['plot_type'],
                gene_list=gene_sel,
                max_genes=100,
                min_unique_per_contrast=params['min_unique'],
                p_value_threshold=params['pvalue_thresh'],
                lfc_threshold=params['lfc_thresh'],
                max_contrasts_for_unique=params['max_contrasts_unique'],
                hide_x_labels=params['hide_x_labels'],
                output_dir=tmpdir
            )

            # Create a zip file
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            zip_filename = f"uorca_explorer_export_{timestamp}"
            zip_path = os.path.join(tmpdir, f"{zip_filename}.zip")
            shutil.make_archive(os.path.join(tmpdir, zip_filename), 'zip', output_dir)

            # Provide download link
            with open(zip_path, "rb") as f:
                st.sidebar.download_button(
                    label="Download HTML Report",
                    data=f,
                    file_name="uorca_explorer_export.zip",
                    mime="application/zip"
                )
    except Exception as e:
        logger.error(f"Error exporting HTML: {str(e)}", exc_info=True)
        st.sidebar.error(f"Error exporting HTML: {str(e)}")


@log_streamlit_function
def _export_csv(ri: ResultsIntegrator, gene_sel: List[str]):
    """Export current view as CSV."""
    try:
        import io
        import zipfile

        # Create a buffer for the zip file
        buffer = io.BytesIO()

        with zipfile.ZipFile(buffer, 'w') as zf:
            # Export gene expression data
            selected_datasets = st.session_state.get('selected_datasets', set())
            for dataset_id in selected_datasets:
                if dataset_id in ri.cpm_data:
                    df = ri.cpm_data[dataset_id]
                    if 'Gene' in df.columns:
                        # Filter to selected genes if any
                        if gene_sel:
                            df = df[df['Gene'].isin(gene_sel)]
                        csv_data = df.to_csv(index=False)
                        zf.writestr(f"{dataset_id}_expression.csv", csv_data)

            # Export DEG data
            selected_contrasts = st.session_state.get('selected_contrasts', set())
            for aid, cid in selected_contrasts:
                if aid in ri.deg_data and cid in ri.deg_data[aid]:
                    df = ri.deg_data[aid][cid]
                    if 'Gene' in df.columns:
                        # Filter to selected genes if any
                        if gene_sel:
                            df = df[df['Gene'].isin(gene_sel)]
                        csv_data = df.to_csv(index=False)
                        zf.writestr(f"{aid}_{cid}_DEG.csv", csv_data)

        # Reset buffer position
        buffer.seek(0)

        # Provide download button
        st.sidebar.download_button(
            label="Download CSV Data",
            data=buffer,
            file_name="uorca_explorer_data.zip",
            mime="application/zip"
        )
    except Exception as e:
        logger.error(f"Error exporting CSV: {str(e)}", exc_info=True)
        st.sidebar.error(f"Error exporting CSV: {str(e)}")


@log_streamlit_function
def _render_help_section():
    """Render the help section at the bottom of the sidebar."""
    st.sidebar.divider()
    st.sidebar.markdown(
        """
        *Powered by Streamlit & Plotly*

        **How to access this app remotely:**
        1. Run this app on your server: `streamlit run uorca_explorer.py`
        2. Create an SSH tunnel from your local machine:
           ```
           ssh -L 8000:127.0.0.1:8501 user@server
           ```
        3. Open `http://127.0.0.1:8000` in your browser
        """
    )

#!/usr/bin/env python3
"""
UORCA Explorer - A Streamlit app for interactive exploration of UORCA RNA-seq results

This app allows you to:
- Select any subset of genes, datasets, or contrasts
- Instantly see updated heatmaps and expression plots
- Explore your RNA-seq results interactively

Usage:
  streamlit run uorca_explorer.py

Or with custom port:
  streamlit run uorca_explorer.py --server.port 8501
"""

import os
import streamlit as st
import pandas as pd
import plotly.express as px
from pathlib import Path
import sys
from datetime import datetime

# Use current directory to import ResultsIntegration
script_dir = os.path.dirname(os.path.abspath(__file__))
# ResultsIntegration.py is in the same directory
from ResultsIntegration import ResultsIntegrator

# Set page configuration
st.set_page_config(
    page_title="UORCA Explorer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Add CSS for dataframe text wrapping
st.markdown("""
<style>
  /* Fix for text wrapping in dataframes */
  .stDataFrame tbody tr td { 
    white-space: normal !important; 
    word-wrap: break-word !important;
    max-width: 300px;
  }
  .stDataFrame th { 
    white-space: normal !important; 
    word-wrap: break-word !important;
    max-width: 300px;
  }
  /* Ensure Description column has more space */
  .stDataFrame td:nth-child(3) {
    min-width: 250px;
    max-width: 500px;
  }
</style>
""", unsafe_allow_html=True)

# ---------- 1. sidebar: data & user parameters -----------------
st.sidebar.title("🧬 UORCA Explorer")

# Get the default results directory
default_dir = os.path.join(os.path.dirname(os.path.dirname(script_dir)), "UORCA_results")
if not os.path.exists(default_dir):
    default_dir = os.path.dirname(os.path.dirname(script_dir))

# absolute path on the server
results_dir = st.sidebar.text_input(
    "Results directory",
    value=default_dir
)

# expensive → cache
@st.cache_data
def load_integrator(path):
    try:
        ri = ResultsIntegrator(results_dir=path)
        ri.load_data()
        return ri, None
    except Exception as e:
        return None, str(e)

# Load the integrator (cached)
with st.sidebar.status("Loading data...", expanded=True) as status:
    ri, error = load_integrator(results_dir)

    if error:
        status.update(label=f"Error loading data: {error}", state="error")
    elif not ri or not ri.cpm_data:
        status.update(label="No data found. Please check the directory path.", state="error")
    else:
        status.update(label=f"✅ Loaded {len(ri.cpm_data)} datasets", state="complete")

# Only show the rest of the UI if we successfully loaded data
if ri and ri.cpm_data:
    # Get all genes across all datasets
    all_genes = set()
    for cpm_df in ri.cpm_data.values():
        if 'Gene' in cpm_df.columns:
            all_genes.update(cpm_df['Gene'].tolist())
    all_genes = sorted(all_genes)

    # Advanced options (moved up so variables are available throughout)
    st.sidebar.subheader("Advanced Options")
    show_advanced = st.sidebar.checkbox("Show advanced options", value=False)

    # Initialize default values
    pvalue_thresh = 0.05
    lfc_thresh = 1.0
    min_unique = 1
    max_contrasts_unique = 2
    top_frequent_genes = 20
    top_unique_genes = 10
    hide_x_labels = True
    debug_mode = False

    if show_advanced:
        # Add slider and number input side by side for precision
        col1, col2 = st.sidebar.columns([3, 1])
        with col1:
            pvalue_thresh = st.slider("Adjusted P-value threshold", 0.0, 0.1, 0.05, 0.001, help="Use adjusted P-value (adj.P.Val) for filtering DEGs")
        with col2:
            pvalue_text = st.text_input("", value=f"{pvalue_thresh:.3f}", key="pvalue_text")
            try:
                pvalue_thresh = float(pvalue_text)
            except ValueError:
                pass

        col1, col2 = st.sidebar.columns([3, 1])
        with col1:
            lfc_thresh = st.slider("Log2FC threshold", 0.0, 5.0, 1.0, 0.1, help="Absolute log2 fold change threshold")
        with col2:
            lfc_text = st.text_input("", value=f"{lfc_thresh:.1f}", key="lfc_text")
            try:
                lfc_thresh = float(lfc_text)
            except ValueError:
                pass

        st.sidebar.subheader("Gene Selection Parameters")
        col1, col2 = st.sidebar.columns([3, 1])
        with col1:
            top_frequent_genes = st.slider("Top frequently occurring genes", 5, 50, 20, help="Number of genes that appear across multiple contrasts")
        with col2:
            freq_text = st.text_input("", value=str(top_frequent_genes), key="freq_text")
            try:
                top_frequent_genes = int(freq_text)
            except ValueError:
                pass

        col1, col2 = st.sidebar.columns([3, 1])
        with col1:
            top_unique_genes = st.slider("Top unique genes per contrast", 1, 20, 10, help="Number of genes unique to each contrast to include")
        with col2:
            unique_text = st.text_input("", value=str(top_unique_genes), key="unique_text")
            try:
                top_unique_genes = int(unique_text)
            except ValueError:
                pass

        col1, col2 = st.sidebar.columns([3, 1])
        with col1:
            max_contrasts_unique = st.slider("Max contrasts for 'unique' genes", 1, 10, 2, help="Maximum number of contrasts a gene can appear in to be considered 'unique'")
        with col2:
            max_text = st.text_input("", value=str(max_contrasts_unique), key="max_text")
            try:
                max_contrasts_unique = int(max_text)
            except ValueError:
                pass

        col1, col2 = st.sidebar.columns([3, 1])
        with col1:
            min_unique = st.slider("Min unique genes per contrast", 0, 10, 1, help="Minimum number of unique genes to select from each contrast")
        with col2:
            min_text = st.text_input("", value=str(min_unique), key="min_text")
            try:
                min_unique = int(min_text)
            except ValueError:
                pass

        st.sidebar.subheader("Visualization Options")
        hide_x_labels = st.sidebar.checkbox("Hide x-axis labels in expression plots", value=True)
        debug_mode = st.sidebar.checkbox("Debug mode", value=False)

    # Initialize session state for selections
    if 'datasets_selected' not in st.session_state:
        st.session_state['datasets_selected'] = None

    if 'contrasts_selected' not in st.session_state:
        st.session_state['contrasts_selected'] = None

    # Gene selection
    st.sidebar.subheader("Gene Selection")
    gene_select_method = st.sidebar.radio(
        "Selection method:",
        ["From list", "Search", "Important genes"]
    )

    if gene_select_method == "From list":
        gene_sel = st.sidebar.multiselect(
            "Select genes:",
            all_genes,
            max_selections=30
        )
    elif gene_select_method == "Search":
        gene_search = st.sidebar.text_input("Search genes:")
        matching_genes = [g for g in all_genes if gene_search.upper() in g.upper()]
        if len(matching_genes) > 100:
            st.sidebar.warning(f"Found {len(matching_genes)} matches. Showing first 100.")
            matching_genes = matching_genes[:100]
        gene_sel = st.sidebar.multiselect(
            "Matching genes:",
            matching_genes,
            max_selections=30
        )
    else:  # Important genes
        # Use p-value and lfc thresholds (these will be defined either here or in advanced options)
        if not show_advanced:
            # Add slider and number input side by side for precision
            col1, col2 = st.sidebar.columns([3, 1])
            with col1:
                pvalue_thresh = st.slider("Adjusted P-value threshold", 0.0, 0.1, 0.05, 0.001, help="Use adjusted P-value (adj.P.Val) for filtering DEGs")
            with col2:
                pvalue_text = st.text_input("", value=f"{pvalue_thresh:.3f}", key="pvalue_text_simple")
                try:
                    pvalue_thresh = float(pvalue_text)
                except ValueError:
                    pass

            col1, col2 = st.sidebar.columns([3, 1])
            with col1:
                lfc_thresh = st.slider("abs(Log2FC) threshold", 0.0, 5.0, 1.0, 0.1, help="Absolute log2 fold change threshold")
            with col2:
                lfc_text = st.text_input("", value=f"{lfc_thresh:.1f}", key="lfc_text_simple")
                try:
                    lfc_thresh = float(lfc_text)
                except ValueError:
                    pass

            # Set default values for parameters with number inputs
            col1, col2 = st.sidebar.columns([3, 1])
            with col1:
                top_frequent_genes = st.slider("Top frequently occurring genes", 5, 50, 20, help="Number of genes that appear across multiple contrasts")
            with col2:
                freq_text = st.text_input("", value=str(top_frequent_genes), key="freq_text_simple")
                try:
                    top_frequent_genes = int(freq_text)
                except ValueError:
                    pass

            col1, col2 = st.sidebar.columns([3, 1])
            with col1:
                top_unique_genes = st.slider("Top unique genes per contrast", 1, 20, 10, help="Number of genes unique to each contrast to include")
            with col2:
                unique_text = st.text_input("", value=str(top_unique_genes), key="unique_text_simple")
                try:
                    top_unique_genes = int(unique_text)
                except ValueError:
                    pass

            col1, col2 = st.sidebar.columns([3, 1])
            with col1:
                max_contrasts_unique = st.slider("Max contrasts for 'unique' genes", 1, 10, 2, help="Maximum number of contrasts a gene can appear in to be considered 'unique'")
            with col2:
                max_text = st.text_input("", value=str(max_contrasts_unique), key="max_text_simple")
                try:
                    max_contrasts_unique = int(max_text)
                except ValueError:
                    pass

        # Get all significant genes
        top_genes = ri.identify_important_genes(
            top_frequent=top_frequent_genes,
            top_unique=top_unique_genes,
            max_contrasts_for_unique=max_contrasts_unique,
            min_unique_per_contrast=min_unique if show_advanced else 1,
            p_value_threshold=pvalue_thresh,
            lfc_threshold=lfc_thresh
        )
        top_genes = sorted(top_genes)
        gene_sel = st.sidebar.multiselect(
            "Important genes:",
            top_genes,
            default=top_genes[:min(50, len(top_genes))] if top_genes else None,
            max_selections=100
        )

    # Dataset selection
    st.sidebar.subheader("Dataset Selection")
    all_dsets = list(ri.cpm_data.keys())

    # Use session state if available
    default_datasets = st.session_state.get('datasets_selected', all_dsets)
    if default_datasets is None:
        default_datasets = all_dsets

    ds_sel = st.sidebar.multiselect(
        "Select datasets:",
        all_dsets,
        default=default_datasets
    )

    # Contrast selection
    st.sidebar.subheader("Contrast Selection")
    all_contr = []
    contrast_map = {}  # Map display string to (analysis_id, contrast_id)
    for aid, contrasts in ri.deg_data.items():
        for cid in contrasts.keys():
            # Get original contrast name and description if available
            original_name = cid
            desc = ""
                
            if hasattr(ri, "contrast_info") and cid in ri.contrast_info:
                # Use name from contrast_info if available
                if 'name' in ri.contrast_info[cid]:
                    original_name = ri.contrast_info[cid]['name']
                    
                # Add description snippet
                if 'description' in ri.contrast_info[cid]:
                    desc_text = ri.contrast_info[cid]['description']
                    if len(desc_text) > 30:
                        desc = f" - {desc_text[:30]}..."
                    else:
                        desc = f" - {desc_text}"
                
            # Create display string with original name
            display_string = f"{aid}:{original_name}{desc}"
            all_contr.append(display_string)
            contrast_map[display_string] = (aid, cid)

    # Use session state if available
    default_contrasts = st.session_state.get('contrasts_selected',
                                           all_contr[:20] if len(all_contr) > 20 else all_contr)
    if default_contrasts is None:
        default_contrasts = all_contr[:15] if len(all_contr) > 15 else all_contr

    # Filter default_contrasts to only include valid contrasts
    default_contrasts = [c for c in default_contrasts if c in all_contr]

    contr_sel = st.sidebar.multiselect(
        "Select contrasts:",
        all_contr,
        default=default_contrasts
    )

    # Visualization options
    st.sidebar.subheader("Plot Options")
    plot_type = st.sidebar.radio(
        "Expression plot type:",
        ["violin", "box", "both"]
    )
    
    # Add pagination controls if needed
    if gene_sel:
        genes_per_page = 30
        total_pages = (len(gene_sel) + genes_per_page - 1) // genes_per_page
        
        if total_pages > 1:
            st.sidebar.markdown("---")
            st.sidebar.subheader("Gene Pagination")
            
            st.sidebar.info(f"Your selection contains {len(gene_sel)} genes, showing {genes_per_page} per page.")
            
            # Create three columns for pagination controls
            col1, col2, col3 = st.sidebar.columns([1, 2, 1])
            
            with col1:
                prev_disabled = 'page_num' in st.session_state and st.session_state.page_num <= 1
                if st.button("◀", disabled=prev_disabled, key="prev_page"):
                    st.session_state.page_num = max(1, st.session_state.get('page_num', 1) - 1)
                    try:
                        st.rerun()
                    except AttributeError:
                        # Fallback for older Streamlit versions
                        try:
                            st.experimental_rerun()
                        except AttributeError:
                            import streamlit.runtime.scriptrunner.magic as _m
                            _m._set_stop_thread(False)  # hidden API for very old versions
            
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
                next_disabled = 'page_num' in st.session_state and st.session_state.page_num >= total_pages
                if st.button("▶", disabled=next_disabled, key="next_page"):
                    st.session_state.page_num = min(total_pages, st.session_state.get('page_num', 1) + 1)
                    try:
                        st.rerun()
                    except AttributeError:
                        # Fallback for older Streamlit versions
                        try:
                            st.experimental_rerun()
                        except AttributeError:
                            import streamlit.runtime.scriptrunner.magic as _m
                            _m._set_stop_thread(False)  # hidden API for very old versions
                    
            # Show page indicator
            st.sidebar.caption(f"Page {st.session_state.page_num} of {total_pages}")
        else:
            if 'page_num' not in st.session_state:
                st.session_state.page_num = 1
    else:
        if 'page_num' not in st.session_state:
            st.session_state.page_num = 1

    # These options have been moved up in the UI

    # ---------- 2. main tabs ---------------------------------------
    tab1, tab2, tab3, tab4, tab5 = st.tabs(["Heat-map", "Expression", "Analysis Plots", "Datasets", "Contrasts"])

    with tab1:
        st.header("Differential Expression Heatmap")

        if not gene_sel:
            st.info("Please select at least one gene from the sidebar.")
        elif not contr_sel:
            st.info("Please select at least one contrast from the sidebar.")

        # Add option to show gene data directly
        if show_advanced and debug_mode:
            if st.checkbox("Show raw gene data", value=False):
                # Get the DEG data for selected genes and contrasts
                gene_data = {}
                for aid, cid in sel_pairs:
                    if aid in ri.deg_data and cid in ri.deg_data[aid]:
                        df = ri.deg_data[aid][cid]
                        if 'Gene' in df.columns and df[df['Gene'].isin(gene_sel)].shape[0] > 0:
                            gene_data[f"{aid}:{cid}"] = df[df['Gene'].isin(gene_sel)]

                # Display the data
                if gene_data:
                    for contrast, df in gene_data.items():
                        st.write(f"### Data for {contrast}")
                        st.dataframe(df)
                else:
                    st.warning("No matching gene data found in selected contrasts")
        else:
            # Use the contrast map to get the correct analysis_id and contrast_id pairs
            sel_pairs = [contrast_map[contr] for contr in contr_sel if contr in contrast_map]

            # Keep original contrasts for the heatmap
            # We'll handle the display in the heatmap function itself

            # Create and display the heatmap
            with st.spinner("Generating heatmap..."):
                try:
                    if show_advanced and debug_mode:
                        st.info("Debug mode: Showing detailed heatmap generation info")
                        import logging
                        import io

                        # Set up a string IO to capture log messages
                        log_stream = io.StringIO()
                        handler = logging.StreamHandler(log_stream)
                        handler.setLevel(logging.DEBUG)
                        formatter = logging.Formatter('%(levelname)s - %(message)s')
                        handler.setFormatter(formatter)

                        # Add the handler to the logger
                        logger = logging.getLogger("ResultsIntegration")
                        logger.setLevel(logging.DEBUG)
                        logger.addHandler(handler)

                    # Create heatmap with possibly modified parameters
                    fig = ri.create_lfc_heatmap(
                        genes=gene_sel,
                        contrasts=sel_pairs,
                        output_file=None
                    )

                    # Rename contrast labels to remove dataset prefix if there are multiple datasets
                    # The ResultsIntegrator already simplified the contrast labels
                    # No need to update x-axis labels here

                    if show_advanced and debug_mode:
                        # Remove the handler to avoid duplicates
                        logger.removeHandler(handler)

                        # Display the log
                        st.expander("Debug Log", expanded=True).code(log_stream.getvalue())

                    if fig:
                        # Add a note about hover functionality
                        st.info("💡 Hover over cells in the heatmap to see contrast descriptions and gene information.")
                        # Display the plot
                        st.plotly_chart(fig, use_container_width=True)
                    else:
                        st.error("Could not generate heatmap. Please check your selections.")
                except Exception as e:
                    st.error(f"Error generating heatmap: {str(e)}")

    with tab2:
        st.header("Gene Expression Plots")

        if not gene_sel:
            st.info("Please select at least one gene from the sidebar.")
        elif not ds_sel:
            st.info("Please select at least one dataset from the sidebar.")
        else:
            # Create and display the expression plots
            with st.spinner("Generating expression plots..."):
                try:
                    # Calculate gene slice for the current page
                    genes_per_page = 30
                    current_page = st.session_state.page_num
                    start_idx = (current_page - 1) * genes_per_page
                    end_idx = min(start_idx + genes_per_page, len(gene_sel))
                    current_genes = gene_sel[start_idx:end_idx]
                    
                    fig2 = ri.create_expression_plots(
                        genes=current_genes,
                        analyses=ds_sel,
                        plot_type=plot_type,
                        output_file=None,
                        hide_x_labels=hide_x_labels,
                        page_number=st.session_state.page_num
                    )
                    if fig2:
                        # Add page navigation info if we have multiple pages
                        if total_pages > 1:
                            # Create pagination controls at the top of the plot for convenience
                            cols = st.columns([2, 1, 1, 1, 2])
                            with cols[1]:
                                prev_disabled = st.session_state.page_num <= 1
                                if st.button("◀ Previous", disabled=prev_disabled, key="prev_main"):
                                    st.session_state.page_num = max(1, st.session_state.page_num - 1)
                                    try:
                                        st.rerun()
                                    except AttributeError:
                                        # Fallback for older Streamlit versions
                                        try:
                                            st.experimental_rerun()
                                        except AttributeError:
                                            import streamlit.runtime.scriptrunner.magic as _m
                                            _m._set_stop_thread(False)  # hidden API for very old versions
                            with cols[2]:
                                st.markdown(f"**Page {st.session_state.page_num}/{total_pages}**")
                            with cols[3]:
                                next_disabled = st.session_state.page_num >= total_pages
                                if st.button("Next ▶", disabled=next_disabled, key="next_main"):
                                    st.session_state.page_num = min(total_pages, st.session_state.page_num + 1)
                                    try:
                                        st.rerun()
                                    except AttributeError:
                                        # Fallback for older Streamlit versions
                                        try:
                                            st.experimental_rerun()
                                        except AttributeError:
                                            import streamlit.runtime.scriptrunner.magic as _m
                                            _m._set_stop_thread(False)  # hidden API for very old versions
                        
                        st.plotly_chart(fig2, use_container_width=True)
                    else:
                        st.error("Could not generate expression plots. Please check your selections.")
                except Exception as e:
                    st.error(f"Error generating expression plots: {str(e)}")

    with tab3:
        st.header("RNA-seq Analysis Plots")

        # Select a dataset first
        dataset_options = list(ri.cpm_data.keys()) if ri else []
        selected_dataset = st.selectbox("Select a dataset to view analysis plots:", dataset_options)

        if not selected_dataset:
            st.info("Please select a dataset to view RNA-seq analysis plots.")
        else:
            # Find QC plots for this dataset
            base_path = os.path.join(results_dir, selected_dataset, "RNAseqAnalysis")

            # Check if the path exists
            if not os.path.exists(base_path):
                st.error(f"Could not find RNAseqAnalysis folder for {selected_dataset}")
            else:
                # General QC plots
                st.subheader("Quality Control and Normalization")

                qc_plot_files = {
                    "MDS Plot": {"file": os.path.join(base_path, "MDS.png"),
                                "description": "**Multi-dimensional scaling plot** showing sample relationships based on their gene expression profiles. This is analogous to a PCA plot.\n\n- **What to look for**: Samples from the same experimental group should cluster together. Outliers may indicate technical issues.\n\n- **Technical details**: Distances on the plot represent leading log2-fold changes between samples (the average of the largest absolute log2-fold changes)."},
                    "Filtering Density": {"file": os.path.join(base_path, "filtering_density.png"),
                                         "description": "**Density plot** showing distribution of gene expression levels (log-CPM) before and after filtering low-expression genes.\n\n- **What to look for**: After filtering (red curve), the left peak of very low-expressed genes should be reduced or eliminated.\n\n- **Technical details**: Low-expression genes are filtered out because they provide little statistical power for differential expression analysis and can increase the multiple testing burden."},
                    "Normalization Boxplots": {"file": os.path.join(base_path, "normalization_boxplots.png"),
                                    "description": "**Boxplots** showing expression distributions before (left) and after (right) normalization.\n\n- **What to look for**: After normalization, all samples should have similar median expression and comparable distributions.\n\n- **Technical details**: Normalization adjusts for technical differences in sequencing depth and composition between libraries to make samples comparable."},
                    "Mean-Variance Relationship (voom)": {"file": os.path.join(base_path, "voom_mean_variance.png"),
                                                 "description": "**Voom transformation plot** showing how the variance of gene expression depends on the mean expression level.\n\n- **What to look for**: A smooth decreasing trend where variance is higher for low-expressed genes and stabilizes for highly-expressed genes.\n\n- **Technical details**: The voom method transforms count data to log-CPM values and estimates the mean-variance relationship to assign appropriate weights for linear modeling."},
                    "SA Plot": {"file": os.path.join(base_path, "sa_plot.png"),
                              "description": "**Sigma vs Average plot** showing the standard deviation of normalized expression values against their average expression.\n\n- **What to look for**: A smooth trend without unusual patterns or outliers.\n\n- **Technical details**: This diagnostic plot from limma helps visualize how gene-wise variances change with expression level after fitting the linear model."}
                }

                # Display QC plots in a structured layout with tabs
                plot_tabs = st.tabs(["Quality Overview", "Expression Filtering", "Normalization", "Variance Modeling"])

                # Tab 1: Quality Overview (MDS plot)
                with plot_tabs[0]:
                    if os.path.exists(qc_plot_files["MDS Plot"]["file"]):
                        st.image(qc_plot_files["MDS Plot"]["file"])
                        with st.expander("What does this plot mean?", expanded=True):
                            st.markdown(qc_plot_files["MDS Plot"]["description"])
                    else:
                        st.info("MDS plot not available for this dataset.")

                # Tab 2: Expression Filtering
                with plot_tabs[1]:
                    if os.path.exists(qc_plot_files["Filtering Density"]["file"]):
                        st.image(qc_plot_files["Filtering Density"]["file"])
                        with st.expander("What does this plot mean?", expanded=True):
                            st.markdown(qc_plot_files["Filtering Density"]["description"])
                    else:
                        st.info("Filtering density plot not available for this dataset.")

                # Tab 3: Normalization
                with plot_tabs[2]:
                    if os.path.exists(qc_plot_files["Normalization Boxplots"]["file"]):
                        st.image(qc_plot_files["Normalization Boxplots"]["file"])
                        with st.expander("What does this plot mean?", expanded=True):
                            st.markdown(qc_plot_files["Normalization Boxplots"]["description"])
                    else:
                        st.info("Normalization boxplots not available for this dataset.")

                # Tab 4: Variance Modeling
                with plot_tabs[3]:
                    col1, col2 = st.columns(2)

                    with col1:
                        if os.path.exists(qc_plot_files["Mean-Variance Relationship (voom)"]["file"]):
                            st.subheader("Voom Transformation")
                            st.image(qc_plot_files["Mean-Variance Relationship (voom)"]["file"])
                            with st.expander("What does this plot mean?"):
                                st.markdown(qc_plot_files["Mean-Variance Relationship (voom)"]["description"])

                    with col2:
                        if os.path.exists(qc_plot_files["SA Plot"]["file"]):
                            st.subheader("Variance Trend")
                            st.image(qc_plot_files["SA Plot"]["file"])
                            with st.expander("What does this plot mean?"):
                                st.markdown(qc_plot_files["SA Plot"]["description"])

                # Contrast-specific plots
                st.markdown("---")
                st.subheader("Differential Expression Results")

                # Find all contrast directories
                contrast_dirs = [d for d in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, d)) and d != "logs"]

                if contrast_dirs:
                    # Find contrast descriptions if available
                    contrast_descriptions = {}
                    for c_id in contrast_dirs:
                        if hasattr(ri, "contrast_info") and c_id in ri.contrast_info:
                            contrast_descriptions[c_id] = ri.contrast_info[c_id].get('description', '')

                    # Create a list of contrast options with descriptions
                    contrast_options = []
                    for c_id in contrast_dirs:
                        label = f"{c_id}"
                        if c_id in contrast_descriptions and contrast_descriptions[c_id]:
                            # Truncate description if too long
                            desc = contrast_descriptions[c_id]
                            if len(desc) > 80:
                                desc = desc[:77] + "..."
                            label += f" - {desc}"
                        contrast_options.append(label)

                    # Create a mapping from display string to actual contrast ID
                    contrast_display_to_id = {}
                    for i, opt in enumerate(contrast_options):
                        contrast_display_to_id[opt] = contrast_dirs[i]
                    
                    # Debug: Log available contrasts
                    if debug_mode:
                        st.write(f"Available contrasts: {', '.join(contrast_dirs)}")

                    selected_contrast_display = st.selectbox("Select a contrast:", contrast_options)
                    selected_contrast = contrast_display_to_id[selected_contrast_display]

                    if selected_contrast:
                        contrast_path = os.path.join(base_path, selected_contrast)

                        # Show the full contrast description if available
                        if selected_contrast in contrast_descriptions and contrast_descriptions[selected_contrast]:
                            st.markdown(f"**Description:** {contrast_descriptions[selected_contrast]}")

                        # Define plot information with improved descriptions
                        contrast_plot_files = {
                            "Heatmap of Top DEGs": {"file": os.path.join(contrast_path, "heatmap_top50.png"),
                                      "description": "**Hierarchical clustering heatmap** showing the top 50 differentially expressed genes for this contrast.\n\n- **What to look for**: Distinct expression patterns that separate the experimental groups. Groups of co-regulated genes (clusters) may represent functional modules.\n\n- **Technical details**: Genes are clustered based on their expression similarity. Color intensity represents expression level (red = high, blue = low)."},
                            "MA Plot": {"file": os.path.join(contrast_path, "ma_plot.png"),
                                      "description": "**MA plot** showing the relationship between average expression level (x-axis) and log2 fold change (y-axis).\n\n- **What to look for**: Symmetrical distribution around y=0 with significant genes (red points) distributed across expression levels. Asymmetry might indicate normalization issues.\n\n- **Technical details**: Each point represents a gene. M (y-axis) is the log-ratio of expression (log2 fold change), and A (x-axis) is the average expression. Red points are statistically significant DEGs."},
                            "Volcano Plot": {"file": os.path.join(contrast_path, "volcano_plot.png"),
                                          "description": "**Volcano plot** showing the relationship between statistical significance (-log10 p-value on y-axis) and biological significance (log2 fold change on x-axis).\n\n- **What to look for**: Genes in the upper left and right corners are both statistically significant and have large fold changes, making them the most interesting candidates.\n\n- **Technical details**: Each point represents a gene. Red points are statistically significant DEGs after multiple testing correction."}
                        }

                        # Create tabs for the contrast-specific plots
                        de_plot_tabs = st.tabs(["Volcano Plot", "MA Plot", "Heatmap"])

                        # Tab 1: Volcano Plot
                        with de_plot_tabs[0]:
                            if os.path.exists(contrast_plot_files["Volcano Plot"]["file"]):
                                st.image(contrast_plot_files["Volcano Plot"]["file"])
                                with st.expander("What does this plot mean?", expanded=True):
                                    st.markdown(contrast_plot_files["Volcano Plot"]["description"])
                            else:
                                st.info("Volcano plot not available for this contrast.")

                        # Tab 2: MA Plot
                        with de_plot_tabs[1]:
                            if os.path.exists(contrast_plot_files["MA Plot"]["file"]):
                                st.image(contrast_plot_files["MA Plot"]["file"])
                                with st.expander("What does this plot mean?", expanded=True):
                                    st.markdown(contrast_plot_files["MA Plot"]["description"])
                            else:
                                st.info("MA plot not available for this contrast.")

                        # Tab 3: Heatmap
                        with de_plot_tabs[2]:
                            if os.path.exists(contrast_plot_files["Heatmap of Top DEGs"]["file"]):
                                st.image(contrast_plot_files["Heatmap of Top DEGs"]["file"])
                                with st.expander("What does this plot mean?", expanded=True):
                                    st.markdown(contrast_plot_files["Heatmap of Top DEGs"]["description"])
                            else:
                                st.info("Heatmap not available for this contrast.")

                        # Add link to view DEG table
                        deg_file = os.path.join(contrast_path, "DEG.csv")
                        if os.path.exists(deg_file):
                            st.markdown("---")
                            with st.expander("View Differentially Expressed Genes Table", expanded=False):
                                try:
                                    deg_df = pd.read_csv(deg_file)

                                    # Create a copy for display, preserving original numeric types
                                    display_df = deg_df.copy()
                                    
                                    # Format p-values for display while keeping them numeric for sorting
                                    column_config = {}
                                    if 'adj.P.Val' in display_df.columns:
                                        # Ensure p-values are numeric for proper sorting
                                        try:
                                            display_df['adj.P.Val'] = pd.to_numeric(display_df['adj.P.Val'])
                                            column_config['adj.P.Val'] = st.column_config.NumberColumn(
                                                "adj.P.Val",
                                                format="%.2e",
                                                help="Adjusted p-value (corrected for multiple testing)"
                                            )
                                        except:
                                            st.warning("Could not convert adj.P.Val to numeric format for sorting.")
                                    if 'P.Value' in display_df.columns:
                                        # Ensure p-values are numeric for proper sorting
                                        try:
                                            display_df['P.Value'] = pd.to_numeric(display_df['P.Value'])
                                            column_config['P.Value'] = st.column_config.NumberColumn(
                                                "P.Value",
                                                format="%.2e",
                                                help="Unadjusted p-value"
                                            )
                                        except:
                                            st.warning("Could not convert P.Value to numeric format for sorting.")
                                    if 'logFC' in display_df.columns:
                                        column_config['logFC'] = st.column_config.NumberColumn(
                                            "logFC",
                                            format="%.2f",
                                            help="Log2 fold change"
                                        )
                                    if 'Gene' in display_df.columns:
                                        column_config['Gene'] = st.column_config.TextColumn(
                                            "Gene",
                                            width="auto",
                                            help="Gene symbol or identifier"
                                        )
                                    if 'Description' in display_df.columns:
                                        column_config['Description'] = st.column_config.TextColumn(
                                            "Description",
                                            width="medium",
                                            help="Contrast description"
                                        )

                                    # Show the dataframe with proper formatting
                                    st.dataframe(
                                        display_df, 
                                        use_container_width=True,
                                        column_config=column_config
                                    )

                                    # Add download button for DEG file
                                    csv = pd.read_csv(deg_file).to_csv(index=False)  # Use original values for download
                                    st.download_button(
                                        label="Download DEG Table as CSV",
                                        data=csv,
                                        file_name=f"{selected_dataset}_{selected_contrast}_DEG.csv",
                                        mime="text/csv"
                                    )
                                except Exception as e:
                                    st.error(f"Error loading DEG file: {str(e)}")
                else:
                    st.info("No contrast-specific plots found for this dataset.")

    with tab4:
        st.header("Dataset Information")

        # Create a DataFrame with dataset information
        dataset_info = []
        for analysis_id, info in ri.analysis_info.items():
            # Build a dataset info dictionary
            dataset_dict = {
                "Accession": info.get("accession", "Unknown"),
                "Organism": info.get("organism", "Unknown"),
                "Number of Samples": info.get("number_of_samples", 0),
                "Number of Contrasts": info.get("number_of_contrasts", 0),
                "Dataset ID": analysis_id  # Keep dataset ID but place it last
            }

            # Add dataset_info.txt content split into three columns if available
            if hasattr(ri, "dataset_info") and analysis_id in getattr(ri, "dataset_info", {}):
                # Remove any "Title:" prefix from the title field
                title = ri.dataset_info[analysis_id].get("title", "")
                if isinstance(title, str) and title.startswith("Title:"):
                    title = title[6:].strip()
                dataset_dict["Title"] = title
            
                dataset_dict["Summary"] = ri.dataset_info[analysis_id].get("summary", "")
                dataset_dict["Design"] = ri.dataset_info[analysis_id].get("design", "")

            dataset_info.append(dataset_dict)

        # Display the dataset information
        if dataset_info:
            df = pd.DataFrame(dataset_info)

            # Add filtering options
            st.subheader("Filter Datasets")
            col1, col2 = st.columns(2)

            with col1:
                organism_filter = st.multiselect(
                    "Filter by Organism",
                    options=sorted(df["Organism"].unique()),
                    default=[]
                )

            with col2:
                search_filter = st.text_input("Search Datasets", "")

            # Apply filters
            filtered_df = df
            if organism_filter:
                filtered_df = filtered_df[filtered_df["Organism"].isin(organism_filter)]

            if search_filter:
                search_mask = filtered_df.apply(
                    lambda row: any(search_filter.lower() in str(val).lower() for val in row),
                    axis=1
                )
                filtered_df = filtered_df[search_mask]

            # Display the filtered dataset information
            if not filtered_df.empty:
                # Display dataset information with dataframe for interactivity
                st.dataframe(filtered_df, use_container_width=True, height=400)

                # Show dataset details if available (always show by default)
                if any(col in filtered_df.columns for col in ["Title", "Summary", "Design"]):
                    for _, row in filtered_df.iterrows():
                        dataset_id = row.get("Dataset ID", row.get("Accession", "Unknown"))
                        accession = row.get("Accession", "")
                        with st.expander(f"Details for {dataset_id} {f'({accession})' if accession else ''}", expanded=True):
                            if "Title" in filtered_df.columns and pd.notna(row.get("Title")):
                                st.subheader("Title")
                                st.markdown(str(row["Title"]))
                            
                            if "Summary" in filtered_df.columns and pd.notna(row.get("Summary")):
                                st.subheader("Summary")
                                st.markdown(str(row["Summary"]))
                                
                            if "Design" in filtered_df.columns and pd.notna(row.get("Design")):
                                st.subheader("Overall Design")
                                st.markdown(str(row["Design"]))

                # Add dataset selection button
                if st.button("Select these datasets for analysis"):
                    datasets_to_select = filtered_df["Dataset ID"].tolist() if "Dataset ID" in filtered_df.columns else filtered_df["Accession"].tolist()
                    st.session_state['datasets_selected'] = datasets_to_select
                    # Reset page number when changing datasets
                    st.session_state.page_num = 1
                    st.success(f"Selected {len(datasets_to_select)} datasets for analysis!")
                    st.info("Switch to the Heat-map or Expression tab to view the analysis.")
            else:
                st.info("No datasets match the current filters.")
        else:
            st.info("No dataset information available.")

    with tab5:
        st.header("Contrast Information")

        # Create a DataFrame with contrast information
        # Display contrast information
        contrast_info = []
        for aid, contrasts in ri.deg_data.items():
            for contrast_id in contrasts.keys():
                # Get original contrast name and description
                original_name = contrast_id
                description = "No description available"
                
                # Try to get information from contrast_info
                if hasattr(ri, "contrast_info") and contrast_id in ri.contrast_info:
                    description = ri.contrast_info[contrast_id].get('description', "No description available")
                    # Use original name if available
                    if 'name' in ri.contrast_info[contrast_id]:
                        original_name = ri.contrast_info[contrast_id]['name']

                # Count DEGs for this contrast
                deg_count = 0
                if aid in ri.deg_data and contrast_id in ri.deg_data[aid]:
                    df = ri.deg_data[aid][contrast_id]
                    # Use exact column names from DEG.csv file - prefer adjusted p-value
                    p_value_col = None
                    if 'adj.P.Val' in df.columns:
                        p_value_col = 'adj.P.Val'  # Adjusted p-value (preferred)
                    elif 'P.Value' in df.columns:
                        p_value_col = 'P.Value'  # Unadjusted p-value (fallback)
                    
                    lfc_col = 'logFC' if 'logFC' in df.columns else None

                    if p_value_col and lfc_col:
                        deg_count = ((df[p_value_col] < pvalue_thresh) &
                                     (abs(df[lfc_col]) > lfc_thresh)).sum()

                contrast_info.append({
                    "Dataset": aid,
                    "Accession": ri.analysis_info.get(aid, {}).get("accession", "Unknown"),
                    "Contrast": original_name,
                    "Original ID": contrast_id,
                    "Description": description,
                    "DEGs": deg_count
                })

        # Display the contrast information
        if contrast_info:
            df = pd.DataFrame(contrast_info)
            # Format DEG count column as integer
            if 'DEGs' in df.columns:
                df['DEGs'] = df['DEGs'].astype(int)

            # Add filtering options
            st.subheader("Filter Contrasts")
            col1, col2, col3 = st.columns(3)

            with col1:
                dataset_column = "Dataset" if "Dataset" in df.columns else "Accession"
                dataset_filter = st.multiselect(
                    "Filter by Dataset",
                    options=sorted(df[dataset_column].unique()),
                    default=[]
                )

            with col2:
                min_degs = st.number_input("Minimum DEGs", min_value=0, value=0)

            with col3:
                search_filter = st.text_input("Search Contrasts", "")

            # Apply filters
            filtered_df = df
            if dataset_filter:
                # Get the column to filter by
                dataset_column = "Dataset" if "Dataset" in filtered_df.columns else "Accession"
                filtered_df = filtered_df[filtered_df[dataset_column].isin(dataset_filter)]

            if min_degs > 0:
                filtered_df = filtered_df[filtered_df["DEGs"] >= min_degs]

            if search_filter:
                search_mask = filtered_df.apply(
                    lambda row: any(search_filter.lower() in str(val).lower() for val in row),
                    axis=1
                )
                filtered_df = filtered_df[search_mask]

            # Display the filtered contrast information
            if not filtered_df.empty:
                # Sort by DEG count by default
                filtered_df = filtered_df.sort_values("DEGs", ascending=False)

                # Display contrast information with dataframe for interactivity
                st.dataframe(filtered_df, use_container_width=True, height=400)

                # Add contrast selection button
                if st.button("Select these contrasts for analysis"):
                    # Create list to store selected contrasts
                    contrasts_to_select = []
                    
                    # Get the original contrast ID if available
                    for _, row in filtered_df.iterrows():
                        aid = row['Dataset']
                    
                        # Use Original ID if available, otherwise use Contrast
                        if "Original ID" in row:
                            cid = row["Original ID"]
                        else:
                            cid = row["Contrast"]
                    
                        # Try different formats for matching
                        matched = False
                    
                        # Format 1: Try direct match
                        for display_string, (a_id, c_id) in contrast_map.items():
                            if a_id == aid and c_id == cid:
                                contrasts_to_select.append(display_string)
                                matched = True
                                break
                            
                        # Format 2: Try concatenated form "Dataset:Contrast"
                        if not matched:
                            combined_key = f"{aid}:{cid}"
                            if combined_key in contrast_map:
                                contrasts_to_select.append(combined_key)
                                matched = True
                    
                        # Format 3: Try searching for partial matches
                        if not matched:
                            for display_string in contrast_map.keys():
                                if cid in display_string and (aid in display_string or row.get('Accession', '') in display_string):
                                    contrasts_to_select.append(display_string)
                                    matched = True
                                    break
                                    
                    if contrasts_to_select:
                        st.session_state['contrasts_selected'] = contrasts_to_select
                        # Reset page number when changing contrasts
                        st.session_state.page_num = 1
                        st.success(f"Selected {len(contrasts_to_select)} contrasts for analysis!")
                        st.info("Switch to the Heat-map tab to view the analysis.")
                    else:
                        st.warning("Could not match selected contrasts to display strings. Please try selecting them from the sidebar.")
            else:
                st.info("No contrasts match the current filters.")
        else:
            st.info("No contrast information available.")

    # ---------- 3. additional features -------------------------------
    st.sidebar.divider()
    st.sidebar.subheader("Export Options")
    export_format = st.sidebar.selectbox("Export format:", ["HTML", "CSV"])

    if st.sidebar.button("Export Current View"):
            if export_format == "HTML":
                try:
                    # Create a temporary directory for the export
                    import tempfile
                    import shutil

                    with tempfile.TemporaryDirectory() as tmpdir:
                        # Generate the interactive report
                        output_dir = ri.create_integrated_report(
                            top_frequent=top_frequent_genes,
                            top_unique=top_unique_genes,
                            plot_type=plot_type,
                            gene_list=gene_sel,
                            max_genes=100,
                            min_unique_per_contrast=min_unique,
                            p_value_threshold=pvalue_thresh,
                            lfc_threshold=lfc_thresh,
                            max_contrasts_for_unique=max_contrasts_unique,
                            hide_x_labels=hide_x_labels,
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
                    st.sidebar.error(f"Error exporting HTML: {str(e)}")
            else:
                try:
                    # Create CSV exports
                    import io

                    # Create a buffer for the zip file
                    buffer = io.BytesIO()
                    import zipfile

                    with zipfile.ZipFile(buffer, 'w') as zf:
                        # Export gene expression data
                        for dataset_id in ds_sel:
                            if dataset_id in ri.cpm_data:
                                df = ri.cpm_data[dataset_id]
                                if 'Gene' in df.columns:
                                    # Filter to selected genes if any
                                    if gene_sel:
                                        df = df[df['Gene'].isin(gene_sel)]
                                    csv_data = df.to_csv(index=False)
                                    zf.writestr(f"{dataset_id}_expression.csv", csv_data)

                        # Export DEG data
                        for aid, cid in sel_pairs:
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
                    st.sidebar.error(f"Error exporting CSV: {str(e)}")

    # ---------- 4. housekeeping ------------------------------------
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
else:
    # Show help if no data is loaded
    st.info(
        """
        ## Welcome to UORCA Explorer

        This app allows you to interactively explore RNA-seq results from UORCA analyses.

        ### Getting Started
        1. Enter the path to your UORCA results directory in the sidebar
        2. The app will load your data and display interactive visualizations
        3. Use the sidebar controls to filter genes, datasets, and contrasts

        ### Troubleshooting
        - Make sure the path contains valid UORCA analysis results
        - Each analysis should have a directory structure with:
          - RNAseqAnalysis/ directory
          - CPM.csv file
          - contrasts.csv file
          - edger_analysis_samples.csv file
        """
    )

if __name__ == "__main__":
    # This is just here for documentation
    pass

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
    debug_mode = False
    
    if show_advanced:
        pvalue_thresh = st.sidebar.slider("P-value threshold", 0.0, 0.1, 0.05)
        lfc_thresh = st.sidebar.slider("Log2FC threshold", 0.0, 5.0, 1.0)
        min_unique = st.sidebar.slider("Min unique genes per contrast", 0, 10, 1)
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
            pvalue_thresh = st.sidebar.slider("P-value threshold", 0.0, 0.1, 0.05)
            lfc_thresh = st.sidebar.slider("abs(Log2FC) threshold", 0.0, 5.0, 1.0)
        
        # Get all significant genes
        top_genes = ri.identify_important_genes(
            top_frequent=30,
            top_unique=20,
            max_contrasts_for_unique=2,
            min_unique_per_contrast=min_unique if show_advanced else 1
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
            # Get contrast description if available
            desc = ""
            if hasattr(ri, "contrast_info") and cid in ri.contrast_info:
                desc = f" - {ri.contrast_info[cid]['description'][:30]}..."
            display_string = f"{aid}:{cid}{desc}"
            all_contr.append(display_string)
            contrast_map[display_string] = (aid, cid)
    
    # Use session state if available
    default_contrasts = st.session_state.get('contrasts_selected', 
                                           all_contr[:20] if len(all_contr) > 20 else all_contr)
    if default_contrasts is None:
        default_contrasts = all_contr[:5] if len(all_contr) > 5 else all_contr
    
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
    
    # These options have been moved up in the UI
    
    # ---------- 2. main tabs ---------------------------------------
    tab1, tab2, tab3, tab4 = st.tabs(["Heat-map", "Expression", "Datasets", "Contrasts"])
    
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
                    
                    if show_advanced and debug_mode:
                        # Remove the handler to avoid duplicates
                        logger.removeHandler(handler)
                        
                        # Display the log
                        st.expander("Debug Log", expanded=True).code(log_stream.getvalue())
                    if fig:
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
                    fig2 = ri.create_expression_plots(
                        genes=gene_sel,
                        analyses=ds_sel,
                        plot_type=plot_type,
                        output_file=None
                    )
                    if fig2:
                        st.plotly_chart(fig2, use_container_width=True)
                    else:
                        st.error("Could not generate expression plots. Please check your selections.")
                except Exception as e:
                    st.error(f"Error generating expression plots: {str(e)}")
    
    with tab3:
        st.header("Dataset Information")
        
        # Create a DataFrame with dataset information
        dataset_info = []
        for analysis_id, info in ri.analysis_info.items():
            dataset_info.append({
                "Dataset": analysis_id,
                "Accession": info.get("accession", "Unknown"),
                "Organism": info.get("organism", "Unknown"),
                "Number of Samples": info.get("number_of_samples", 0),
                "Number of Contrasts": info.get("number_of_contrasts", 0),
                "Analysis Column": info.get("analysis_column", "Unknown")
            })
        
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
                st.dataframe(filtered_df, use_container_width=True, height=400)
                
                # Add dataset selection button
                if st.button("Select these datasets for analysis"):
                    datasets_to_select = filtered_df["Dataset"].tolist()
                    st.session_state['datasets_selected'] = datasets_to_select
                    st.success(f"Selected {len(datasets_to_select)} datasets for analysis!")
                    st.info("Switch to the Heat-map or Expression tab to view the analysis.")
            else:
                st.info("No datasets match the current filters.")
        else:
            st.info("No dataset information available.")
    
    with tab4:
        st.header("Contrast Information")
        
        # Create a DataFrame with contrast information
        contrast_info = []
        for aid, contrasts in ri.deg_data.items():
            for cid in contrasts.keys():
                description = "No description available"
                if hasattr(ri, "contrast_info") and cid in ri.contrast_info:
                    description = ri.contrast_info[cid].get('description', "No description available")
                
                # Count DEGs for this contrast
                deg_count = 0
                if aid in ri.deg_data and cid in ri.deg_data[aid]:
                    df = ri.deg_data[aid][cid]
                    # Find p-value and logFC columns
                    p_value_col = next((col for col in ['adj.P.Val', 'padj', 'FDR', 'q_value', 'PValue', 'pvalue', 'P.Value'] 
                                      if col in df.columns), None)
                    lfc_col = next((col for col in ['logFC', 'log2FoldChange', 'log2FC', 'LogFC'] 
                                   if col in df.columns), None)
                    
                    if p_value_col and lfc_col:
                        deg_count = ((df[p_value_col] < pvalue_thresh) & 
                                     (abs(df[lfc_col]) > lfc_thresh)).sum()
                
                contrast_info.append({
                    "Dataset": aid,
                    "Contrast": cid,
                    "Description": description,
                    "DEGs": deg_count
                })
        
        # Display the contrast information
        if contrast_info:
            df = pd.DataFrame(contrast_info)
            
            # Add filtering options
            st.subheader("Filter Contrasts")
            col1, col2, col3 = st.columns(3)
            
            with col1:
                dataset_filter = st.multiselect(
                    "Filter by Dataset",
                    options=sorted(df["Dataset"].unique()),
                    default=[]
                )
            
            with col2:
                min_degs = st.number_input("Minimum DEGs", min_value=0, value=0)
            
            with col3:
                search_filter = st.text_input("Search Contrasts", "")
            
            # Apply filters
            filtered_df = df
            if dataset_filter:
                filtered_df = filtered_df[filtered_df["Dataset"].isin(dataset_filter)]
            
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
                
                st.dataframe(filtered_df, use_container_width=True, height=400)
                
                # Add contrast selection button
                if st.button("Select these contrasts for analysis"):
                    contrasts_to_select = [f"{row['Dataset']}:{row['Contrast']}" for _, row in filtered_df.iterrows()]
                    st.session_state['contrasts_selected'] = contrasts_to_select
                    st.success(f"Selected {len(contrasts_to_select)} contrasts for analysis!")
                    st.info("Switch to the Heat-map tab to view the analysis.")
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
            st.sidebar.info("HTML export is coming soon!")
        else:
            st.sidebar.info("CSV export is coming soon!")
    
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
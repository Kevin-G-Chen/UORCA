#!/usr/bin/env python3
"""
UORCA Landing Page Streamlit App

Interactive interface for generating AI-assisted landing pages with minimal user interaction.
"""

import os
import streamlit as st
import pandas as pd
import tempfile
from pathlib import Path
import logging

# Import the landing page generator
from landing_page_generator import LandingPageGenerator, LandingPageData

# Set page configuration
st.set_page_config(
    page_title="UORCA Landing Page Generator",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Add CSS for better styling
st.markdown("""
<style>
    .metric-card {
        background: linear-gradient(90deg, #667eea 0%, #764ba2 100%);
        padding: 1rem;
        border-radius: 0.5rem;
        color: white;
        text-align: center;
        margin: 0.5rem 0;
    }
    .narrative-box {
        background: #e8f4fd;
        padding: 1.5rem;
        border-radius: 0.5rem;
        border-left: 4px solid #3498db;
        margin: 1rem 0;
    }
    .threshold-box {
        background: #f0f8e8;
        padding: 1rem;
        border-radius: 0.5rem;
        border-left: 4px solid #27ae60;
        margin: 1rem 0;
    }
    .stDataFrame {
        font-size: 0.9rem;
    }
</style>
""", unsafe_allow_html=True)

# Title and description
st.title("🧬 UORCA AI-Assisted Landing Page Generator")
st.markdown("""
**Generate intelligent summaries of your RNA-seq results with minimal interaction.**

This tool automatically selects the most relevant contrasts, determines optimal statistical thresholds, 
and creates interpretive narratives for your differential expression analysis.
""")

# Sidebar controls
st.sidebar.title("⚙️ Configuration")

# Get default results directory
default_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))), "UORCA_results")
if not os.path.exists(default_dir):
    default_dir = ""

# Input controls
results_dir = st.sidebar.text_input(
    "Results Directory",
    value=default_dir,
    help="Path to directory containing UORCA analysis results"
)

biological_prompt = st.sidebar.text_area(
    "Research Context",
    value="General differential expression analysis",
    height=100,
    help="Describe your biological research question or context. This guides the AI in selecting relevant contrasts."
)

# Advanced options
with st.sidebar.expander("🔧 Advanced Options", expanded=False):
    max_contrasts = st.slider(
        "Maximum Contrasts",
        min_value=3,
        max_value=20,
        value=8,
        help="Maximum number of experimental contrasts to include in the analysis"
    )
    
    max_genes = st.slider(
        "Maximum Genes",
        min_value=20,
        max_value=200,
        value=50,
        help="Maximum number of top genes to display in visualizations"
    )

# Generate button
generate_button = st.sidebar.button(
    "🚀 Generate Landing Page",
    type="primary",
    use_container_width=True
)

# Main content area
if not results_dir:
    st.info("👆 Please specify a results directory in the sidebar to get started.")
    st.markdown("### 📋 How to use this tool:")
    st.markdown("""
    1. **Set Results Directory**: Point to your UORCA analysis results folder
    2. **Describe Research Context**: Provide biological context to guide AI selection
    3. **Adjust Options** (optional): Modify the number of contrasts and genes
    4. **Generate**: Click the button to create your landing page
    
    The AI will automatically:
    - Select the most biologically relevant contrasts
    - Determine optimal statistical thresholds  
    - Identify key differentially expressed genes
    - Generate interpretive narratives
    - Create interactive visualizations
    """)

elif not os.path.exists(results_dir):
    st.error(f"❌ Directory not found: {results_dir}")
    st.markdown("Please check the path and ensure the directory exists.")

elif generate_button:
    # Validate directory structure
    if not any(os.path.isdir(os.path.join(results_dir, item)) for item in os.listdir(results_dir)):
        st.error("❌ No analysis subdirectories found in the specified directory.")
        st.stop()
    
    # Generate landing page
    with st.spinner("🤖 AI is analyzing your data and generating the landing page..."):
        try:
            # Initialize generator
            generator = LandingPageGenerator(results_dir, biological_prompt)
            
            # Generate landing page data
            landing_data = generator.generate_landing_page(
                max_contrasts=max_contrasts,
                max_genes=max_genes,
                output_path=None  # Don't save to file yet
            )
            
            # Store in session state
            st.session_state.landing_data = landing_data
            st.session_state.generator = generator
            
        except Exception as e:
            st.error(f"❌ Error generating landing page: {str(e)}")
            st.stop()
    
    st.success("✅ Landing page generated successfully!")

# Display results if available
if hasattr(st.session_state, 'landing_data') and st.session_state.landing_data:
    landing_data = st.session_state.landing_data
    
    # Summary metrics
    st.markdown("### 📊 Analysis Summary")
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.metric("Selected Contrasts", len(landing_data.selected_contrasts))
    with col2:
        st.metric("Key Genes", len(landing_data.top_genes))
    with col3:
        st.metric("FDR Threshold", f"{landing_data.thresholds.fdr_cutoff}")
    with col4:
        st.metric("LogFC Threshold", f"{landing_data.thresholds.logfc_cutoff}")
    
    # Key findings narrative
    st.markdown("### 📝 Key Findings")
    st.markdown(f"""
    <div class="narrative-box">
    {landing_data.narrative}
    </div>
    """, unsafe_allow_html=True)
    
    # Analysis parameters
    st.markdown("### 🎯 Analysis Parameters")
    st.markdown(f"""
    <div class="threshold-box">
    <strong>Statistical Thresholds:</strong> FDR < {landing_data.thresholds.fdr_cutoff}, |log2FC| > {landing_data.thresholds.logfc_cutoff}, Min. frequency ≥ {landing_data.thresholds.min_frequency}<br>
    <strong>Justification:</strong> {landing_data.thresholds.justification}
    </div>
    """, unsafe_allow_html=True)
    
    # Tabs for different views
    tab1, tab2, tab3 = st.tabs(["🔍 Selected Contrasts", "🧬 Top Genes", "🌡️ Expression Heatmap"])
    
    with tab1:
        st.markdown("#### Selected Contrasts and AI Justifications")
        
        # Create contrast DataFrame
        contrast_df = pd.DataFrame([
            {
                "Dataset": c.analysis_id,
                "Contrast": c.contrast_id,
                "Relevance Score": f"{c.relevance_score:.1f}",
                "DEG Count": c.deg_count,
                "AI Justification": c.justification
            }
            for c in landing_data.selected_contrasts
        ])
        
        st.dataframe(
            contrast_df,
            use_container_width=True,
            column_config={
                "AI Justification": st.column_config.TextColumn(
                    "AI Justification",
                    width="large"
                )
            }
        )
    
    with tab2:
        st.markdown("#### Top Differentially Expressed Genes")
        
        if not landing_data.gene_table.empty:
            # Display gene table with better formatting
            display_df = landing_data.gene_table.copy()
            
            st.dataframe(
                display_df,
                use_container_width=True,
                column_config={
                    "Gene": st.column_config.TextColumn("Gene Symbol", width="medium"),
                    "Frequency": st.column_config.NumberColumn("Frequency", help="Number of contrasts where gene is significant"),
                    "Max_LogFC": st.column_config.NumberColumn("Max |Log2FC|", format="%.2f"),
                    "Min_AdjPVal": st.column_config.TextColumn("Min Adj. P-value"),
                    "Contrasts": st.column_config.TextColumn("Contrasts", width="large")
                }
            )
            
            # Gene list for easy copying
            with st.expander("📋 Gene List (for external tools)", expanded=False):
                gene_list_text = '\n'.join(landing_data.top_genes)
                st.text_area(
                    "Copy gene symbols:",
                    value=gene_list_text,
                    height=200,
                    help="Gene symbols, one per line. Copy and paste into other tools."
                )
        else:
            st.info("No genes met the significance criteria.")
    
    with tab3:
        st.markdown("#### Interactive Expression Heatmap")
        
        if landing_data.heatmap_fig:
            st.plotly_chart(landing_data.heatmap_fig, use_container_width=True)
            
            st.markdown("""
            💡 **How to interpret this heatmap:**
            - Each row represents a gene, each column represents a biological contrast
            - Colors indicate log2 fold change: red (upregulated), blue (downregulated), white (not significant)
            - Hover over cells to see detailed information
            - Genes and contrasts are clustered to reveal patterns
            """)
        else:
            st.warning("Could not generate heatmap. This may occur if no genes meet the significance criteria.")
    
    # Export options
    st.markdown("---")
    st.markdown("### 💾 Export Options")
    
    col1, col2 = st.columns(2)
    
    with col1:
        # Export to HTML
        if st.button("📄 Download HTML Report", type="secondary"):
            with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.html') as tmp_file:
                # Generate full HTML
                st.session_state.generator._generate_html(landing_data, tmp_file.name)
                
                # Read the HTML content
                with open(tmp_file.name, 'r', encoding='utf-8') as f:
                    html_content = f.read()
                
                # Provide download
                st.download_button(
                    label="📥 Download Landing Page HTML",
                    data=html_content,
                    file_name=f"uorca_landing_page_{pd.Timestamp.now().strftime('%Y%m%d_%H%M%S')}.html",
                    mime="text/html"
                )
                
                # Clean up
                os.unlink(tmp_file.name)
    
    with col2:
        # Export gene table as CSV
        if not landing_data.gene_table.empty:
            csv_data = landing_data.gene_table.to_csv(index=False)
            st.download_button(
                label="📊 Download Gene Table CSV",
                data=csv_data,
                file_name=f"top_genes_{pd.Timestamp.now().strftime('%Y%m%d_%H%M%S')}.csv",
                mime="text/csv"
            )

# Help and information
st.sidebar.markdown("---")
st.sidebar.markdown("### ❓ Help")
with st.sidebar.expander("🔍 How it works", expanded=False):
    st.markdown("""
    **AI-Assisted Analysis Pipeline:**
    
    1. **Contrast Selection**: AI scores all available contrasts based on biological relevance to your research question
    
    2. **Threshold Optimization**: Statistical thresholds are automatically selected based on dataset characteristics
    
    3. **Gene Ranking**: Genes are ranked by frequency across contrasts and effect size
    
    4. **Narrative Generation**: AI creates interpretive summaries in accessible language
    
    5. **Visualization**: Interactive heatmaps and tables are generated automatically
    """)

with st.sidebar.expander("⚙️ Technical Details", expanded=False):
    st.markdown("""
    **Requirements:**
    - UORCA analysis results directory
    - OpenAI API key (set as environment variable)
    - Results should contain:
      - RNAseqAnalysis folders with DEG.csv files
      - Metadata files with contrast descriptions
    
    **AI Models Used:**
    - GPT-4o-mini for contrast scoring and narrative generation
    - Biological knowledge integration for relevance assessment
    """)

# Footer
st.markdown("---")
st.markdown("""
<div style="text-align: center; color: #7f8c8d; font-size: 0.9em;">
🧬 UORCA Landing Page Generator • AI-Assisted Interpretation with Minimal Interactivity
</div>
""", unsafe_allow_html=True)
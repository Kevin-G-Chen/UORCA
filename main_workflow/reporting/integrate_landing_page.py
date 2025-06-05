#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Integration Guide for MCP-based Landing Page in UORCA Explorer

This file shows how to integrate the new MCP-based reporting agent
into the existing uorca_explorer.py Streamlit application.

The main changes are:
1. Replace the simple LLM-based functions with MCP agent calls
2. Add auto-start functionality when data is loaded
3. Maintain the same UI structure and data formats
"""

# ============================================================================
# STEP 1: Add new imports at the top of uorca_explorer.py
# ============================================================================

# Add these imports after the existing imports
from landing_page_integration import (
    generate_ai_landing_page,
    auto_analyze_on_load,
    LandingPageData,
    ContrastSelection,
    ThresholdSelection
)

# ============================================================================
# STEP 2: Replace the existing AI landing page functions
# ============================================================================

# REMOVE these existing functions from uorca_explorer.py:
# - generate_ai_landing_page()
# - select_contrasts_with_llm()
# - select_genes_simple()
# - create_gene_table_simple()
# - create_landing_heatmap()
# - generate_narrative_llm()

# The new landing_page_integration module provides all this functionality
# through the MCP-based agent system

# ============================================================================
# STEP 3: Add auto-start functionality when data is loaded
# ============================================================================

# In the sidebar data loading section, after successfully loading data:

# Add this after the line: status.update(label=f"✅ Loaded {len(ri.cpm_data)} datasets", state="complete")

# Check if this is the first time loading this directory
if 'previous_results_dir' not in st.session_state or st.session_state.previous_results_dir != results_dir:
    st.session_state.previous_results_dir = results_dir
    st.session_state.auto_analysis_complete = False
    st.session_state.auto_analysis_result = None

# Trigger auto-analysis if not done yet
if not st.session_state.get('auto_analysis_complete', False):
    with st.spinner("🤖 Automatically analyzing your dataset..."):
        auto_result = auto_analyze_on_load(ri)
        if auto_result and auto_result.get('success'):
            st.session_state.auto_analysis_result = auto_result
            st.session_state.auto_analysis_complete = True
            st.success("✅ Initial analysis complete! Check the AI Summary tab.")

# ============================================================================
# STEP 4: Update the AI landing page tab
# ============================================================================

# In the tab_landing section, update the generate button handler:

if generate_landing:
    # Check if required packages and API keys are available
    if not LANDING_PAGE_AVAILABLE:
        st.error("❌ AI landing page functionality requires OpenAI package.")
        st.info("💡 To enable: `pip install openai`")
    elif not api_key_available:
        st.error("❌ OpenAI API key required. Set OPENAI_API_KEY environment variable.")
    else:
        # Set flag to prevent unnecessary fragment execution during AI generation
        st.session_state.ai_generating = True
        
        # Generate AI-assisted landing page using the new MCP-based system
        with st.spinner("🤖 AI agent is analyzing your data..."):
            try:
                # Use the new generate_ai_landing_page function from landing_page_integration
                landing_data = generate_ai_landing_page(
                    integrator=ri,
                    biological_prompt=biological_prompt,
                    max_genes=50
                )
                
                if landing_data:
                    # Store in session state
                    st.session_state.landing_data = landing_data
                    st.success("✅ AI analysis complete!")
                else:
                    st.error("❌ Failed to generate landing page - no suitable data found.")
                    
            except Exception as e:
                st.error(f"❌ Error generating landing page: {str(e)}")
                if "api" in str(e).lower() or "key" in str(e).lower():
                    st.info("💡 Make sure your OpenAI API key is set correctly.")
            finally:
                # Clear the flag after generation completes
                if 'ai_generating' in st.session_state:
                    del st.session_state.ai_generating

# ============================================================================
# STEP 5: Display auto-analysis results if available
# ============================================================================

# Add this section at the beginning of the landing page tab:

# Display auto-analysis if available and no manual analysis has been run
if (not hasattr(st.session_state, 'landing_data') and 
    st.session_state.get('auto_analysis_result') and 
    st.session_state.auto_analysis_result.get('success')):
    
    st.info("🤖 Showing automatic dataset analysis. You can refine it with a specific research question below.")
    
    # Display the auto-analysis report
    auto_report = st.session_state.auto_analysis_result.get('report', '')
    
    st.markdown("### 📊 Dataset Overview")
    st.markdown(auto_report)
    
    # Suggest running a focused analysis
    st.markdown("---")
    st.markdown("### 🎯 Want a More Focused Analysis?")
    st.markdown("Enter a specific research question below to get targeted insights.")

# ============================================================================
# STEP 6: Update the Apply Selections button
# ============================================================================

# The "Apply These Selections to Other Tabs" button remains the same,
# as the data structures are compatible

# ============================================================================
# STEP 7: Environment Setup
# ============================================================================

# Ensure these environment variables are set:
# - OPENAI_API_KEY: Your OpenAI API key
# - UORCA_DEFAULT_RESULTS_DIR: Default results directory (optional)

# For local model support (Ollama), ensure Ollama is running:
# ollama serve

# ============================================================================
# STEP 8: Configuration
# ============================================================================

# The system will look for MCP servers in these locations:
# 1. main_workflow/reporting/mcp_servers/mcp_data_extractor.py
# 2. main_workflow/reporting/mcp_servers/mcp_analysis.py

# Make sure these files exist and are executable

# ============================================================================
# Example: Complete Integration for Generate Button
# ============================================================================

"""
# Replace the existing generate_landing button handler with:

if generate_landing:
    if not LANDING_PAGE_AVAILABLE:
        st.error("❌ AI functionality not available. Install required packages.")
    elif not api_key_available:
        st.error("❌ API key required. Set OPENAI_API_KEY environment variable.")
    else:
        with st.spinner("🤖 AI agent analyzing your RNA-seq data..."):
            try:
                # Progress callback for user feedback
                progress_bar = st.progress(0)
                def update_progress(progress, message):
                    progress_bar.progress(progress)
                    st.caption(message)
                
                # Generate landing page with MCP agent
                landing_data = generate_ai_landing_page(
                    integrator=ri,
                    biological_prompt=biological_prompt,
                    max_genes=50
                )
                
                if landing_data:
                    st.session_state.landing_data = landing_data
                    progress_bar.empty()
                    st.success("✅ Analysis complete!")
                    st.rerun()
                else:
                    st.error("Failed to generate analysis.")
                    
            except Exception as e:
                st.error(f"Error: {str(e)}")
                logger.error(f"Landing page generation failed: {e}")
"""

# ============================================================================
# Testing the Integration
# ============================================================================

"""
To test the integration:

1. Start the MCP servers:
   python main_workflow/reporting/test_mcp_setup.py

2. Run the Streamlit app:
   streamlit run main_workflow/reporting/uorca_explorer.py

3. Load a dataset and verify:
   - Auto-analysis runs on first load
   - Manual analysis works with research questions
   - Visualizations are generated correctly
   - All data structures are compatible
"""
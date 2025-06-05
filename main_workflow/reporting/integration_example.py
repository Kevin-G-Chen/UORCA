#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Minimal Example: Key Integration Changes for MCP-based Landing Page

This file shows the essential changes needed to integrate the new MCP-based
reporting agent into the existing uorca_explorer.py application.
"""

# ============================================================================
# KEY CHANGE 1: Import the new landing page integration module
# ============================================================================

# In uorca_explorer.py, replace the existing AI functions with:
from landing_page_integration import generate_ai_landing_page, auto_analyze_on_load

# ============================================================================
# KEY CHANGE 2: Auto-analysis on data load
# ============================================================================

def handle_data_load(ri, results_dir):
    """Example of auto-analysis integration when data is loaded."""
    
    # Check if this is a new directory
    if st.session_state.get('previous_results_dir') != results_dir:
        st.session_state.previous_results_dir = results_dir
        st.session_state.auto_analysis_done = False
        
    # Run auto-analysis if not done
    if not st.session_state.get('auto_analysis_done', False):
        with st.spinner("🤖 Analyzing dataset automatically..."):
            result = auto_analyze_on_load(ri)
            if result and result.get('success'):
                st.session_state.auto_analysis_result = result
                st.session_state.auto_analysis_done = True
                st.success("✅ Initial analysis ready!")

# ============================================================================
# KEY CHANGE 3: Replace the generate button handler
# ============================================================================

def handle_generate_button(ri, biological_prompt):
    """Example of handling the generate button with new system."""
    
    with st.spinner("🤖 AI agent analyzing RNA-seq data..."):
        try:
            # Call the new MCP-based function
            landing_data = generate_ai_landing_page(
                integrator=ri,
                biological_prompt=biological_prompt,
                max_genes=50
            )
            
            if landing_data:
                # Store in session state (same as before)
                st.session_state.landing_data = landing_data
                st.success("✅ Analysis complete!")
                return True
            else:
                st.error("Failed to generate analysis")
                return False
                
        except Exception as e:
            st.error(f"Error: {str(e)}")
            return False

# ============================================================================
# KEY CHANGE 4: Display auto-analysis results
# ============================================================================

def display_landing_page():
    """Example of displaying results in the landing page tab."""
    
    # Show auto-analysis if available and no manual analysis exists
    if (not hasattr(st.session_state, 'landing_data') and 
        st.session_state.get('auto_analysis_result')):
        
        result = st.session_state.auto_analysis_result
        if result.get('success'):
            st.info("🤖 Showing automatic analysis. Refine with a specific question below.")
            st.markdown("### Dataset Overview")
            st.markdown(result.get('report', 'No report available'))
            st.markdown("---")
    
    # Rest of the landing page UI remains the same
    # The generate_ai_landing_page function returns the same data structure
    # so all existing display code works without changes

# ============================================================================
# MINIMAL INTEGRATION CHECKLIST
# ============================================================================

"""
1. Import the new module:
   from landing_page_integration import generate_ai_landing_page, auto_analyze_on_load

2. Remove old functions:
   - select_contrasts_with_llm()
   - select_genes_simple()
   - create_gene_table_simple()
   - generate_narrative_llm()

3. Add auto-analysis after data load:
   if data_loaded_successfully:
       auto_analyze_on_load(ri)

4. Update generate button:
   landing_data = generate_ai_landing_page(ri, biological_prompt, max_genes=50)

5. Everything else stays the same!
   - Same data structures (LandingPageData, ContrastSelection, etc.)
   - Same visualization code
   - Same UI layout
"""

# ============================================================================
# ENVIRONMENT REQUIREMENTS
# ============================================================================

"""
Required environment variables:
- OPENAI_API_KEY: Your OpenAI API key

Required files:
- main_workflow/reporting/mcp_servers/mcp_data_extractor.py
- main_workflow/reporting/mcp_servers/mcp_analysis.py
- main_workflow/reporting/reporting_agent.py
- main_workflow/reporting/auto_start_manager.py
- main_workflow/reporting/landing_page_integration.py

The system will automatically:
1. Start MCP servers when needed
2. Run initial analysis when data is loaded
3. Generate focused analyses on demand
4. Clean up resources when done
"""
"""
AI Assistant Tab for UORCA Explorer.

This tab provides AI-powered analysis and exploration capabilities.
"""

import os
import json
import asyncio
import logging
import streamlit as st
from typing import Dict, Any, List, Optional


from .helpers import (
    check_ai_generating,
    setup_fragment_decorator,
    log_streamlit_tab,
    log_streamlit_function,
    log_streamlit_agent,
    log_streamlit_event
)
from ResultsIntegration import ResultsIntegrator

# Set up fragment decorator
setup_fragment_decorator()

logger = logging.getLogger(__name__)

# Check for AI functionality availability

try:
    from ai_agent_factory import create_uorca_agent
    UORCA_AGENT_AVAILABLE = True
except ImportError as e:
    logger.warning(f"UORCA agent not available: {e}")
    UORCA_AGENT_AVAILABLE = False

try:
    from contrast_relevance import run_contrast_relevance
    CONTRAST_RELEVANCE_AVAILABLE = True
except ImportError as e:
    logger.warning(f"Contrast relevance not available: {e}")
    CONTRAST_RELEVANCE_AVAILABLE = False


@log_streamlit_tab("AI Assistant")
def render_ai_assistant_tab(ri: ResultsIntegrator, results_dir: str):
    """
    Render the AI assistant tab.

    Args:
        ri: ResultsIntegrator instance
        results_dir: Path to the results directory
    """
    st.header("🤖 AI Assistant")
    st.markdown("**🤖 Interactive AI assistant for exploring your UORCA analysis results.** Ask questions, get insights, and explore your data using natural language.")

    # Check for OpenAI API key
    if not os.getenv("OPENAI_API_KEY"):
        log_streamlit_event("OpenAI API key not found for AI assistant")
        st.error("OPENAI_API_KEY environment variable not set")
        st.info("Please set your OpenAI API key to use AI features.")
        return

    # Render streamlined AI analysis workflow
    _render_streamlined_ai_workflow(ri, results_dir)


@log_streamlit_function
def _render_streamlined_ai_workflow(ri: ResultsIntegrator, results_dir: str):
    """Render the streamlined AI analysis workflow."""
    st.subheader("🧬 AI-Powered Gene Analysis")
    st.markdown("**Complete AI analysis of your RNA-seq data.** Enter your research question below and the AI will assess contrast relevance and identify key genes in one workflow.")

    # Research query input
    research_query = st.text_input(
        "Research Question",
        placeholder="e.g., What contrasts are most relevant to T cell activation and differentiation?",
        help="Describe your research question or area of interest. The AI will score each contrast based on how relevant it is to this query, then analyze key genes."
    )

    # Single workflow button
    col1, col2 = st.columns([1, 3])
    with col1:
        run_button = st.button("🚀 Run Complete AI Analysis", disabled=not research_query.strip(), type="primary")
    with col2:
        if not research_query.strip():
            st.info("💡 Enter a research question above to start AI analysis.")

    if run_button:
        log_streamlit_event(f"User started complete AI analysis: '{research_query.strip()}'")
        _run_complete_ai_analysis(ri, results_dir, research_query.strip())


@log_streamlit_agent
def _run_complete_ai_analysis(ri: ResultsIntegrator, results_dir: str, research_query: str):
    """Run the complete AI analysis workflow including contrast relevance and gene analysis."""
    if not CONTRAST_RELEVANCE_AVAILABLE:
        st.error("Contrast relevance assessment is not available. Please check your environment setup.")
        return

    if not UORCA_AGENT_AVAILABLE:
        st.error("UORCA agent is not available. Please check your environment setup.")
        return

    if not os.getenv("OPENAI_API_KEY"):
        st.error("OpenAI API key not found. Please set the OPENAI_API_KEY environment variable.")
        return

    # Step 1: Contrast Relevance Assessment
    with st.spinner("Step 1/2: Assessing contrast relevance... This may take a few minutes."):
        try:
            # Run contrast relevance assessment
            results_df = run_contrast_relevance(
                ri,
                query=research_query,
                repeats=3,
                batch_size=50,
                parallel_jobs=4
            )

            if not results_df.empty:
                _display_relevance_results(ri, results_df, research_query)

                # Store results for AI gene analysis
                st.session_state['selected_contrasts_for_ai'] = \
                    results_df[['analysis_id','contrast_id']].to_dict('records')
                st.session_state['research_query'] = research_query
            else:
                st.warning("No contrasts found for assessment.")
                return

        except Exception as e:
            st.error(f"Error during contrast relevance assessment: {str(e)}")
            with st.expander("🔍 Error Details", expanded=False):
                import traceback
                st.code(traceback.format_exc())
            return

    # Step 2: AI Gene Analysis
    st.markdown("---")
    with st.spinner("Step 2/2: AI is analyzing differential expression patterns... This may take several minutes."):
        try:
            # Set the RESULTS_DIR environment variable for the MCP server
            os.environ['RESULTS_DIR'] = results_dir

            selected_contrasts = st.session_state['selected_contrasts_for_ai']

            agent = create_uorca_agent()
            prompt = f"""
Research question: "{research_query}"

Available contrasts:
{json.dumps(selected_contrasts, indent=2)}

Please perform the analysis using your four tools, choose all thresholds reasonably, and return:
1) A structured summary showing the key genes identified for each contrast or gene set
2) A brief 2-3 sentence biological interpretation explaining your rationale and what patterns you discovered.
"""

            # Run the agent analysis
            result_text = _execute_ai_analysis(agent, prompt)

            # Display results
            _display_ai_analysis_results(result_text, research_query, selected_contrasts)

        except Exception as e:
            logger.error(f"Error in AI gene analysis: {str(e)}", exc_info=True)
            st.error(f"❌ Analysis failed: {str(e)}")

            with st.expander("🔍 Error Details", expanded=False):
                import traceback
                st.code(traceback.format_exc())


@log_streamlit_function
def _display_relevance_results(ri: ResultsIntegrator, results_df, research_query: str):
    """Display the contrast relevance assessment results."""
    # Sort by relevance score
    results_df = results_df.sort_values('RelevanceScore', ascending=False)

    # Add contrast descriptions for display
    results_df['Description'] = results_df.apply(
        lambda row: ri._get_contrast_description(row['analysis_id'], row['contrast_id']),
        axis=1
    )

    # Add accession info
    results_df['Accession'] = results_df['analysis_id'].map(
        lambda aid: ri.analysis_info.get(aid, {}).get('accession', aid)
    )

    # Add organism info
    results_df['organism'] = results_df['analysis_id'].map(
        lambda aid: ri.analysis_info.get(aid, {}).get('organism', 'Unknown')
    )

    log_streamlit_event(f"Contrast relevance assessment completed: {len(results_df)} contrasts scored")
    st.success(f"✅ Successfully assessed {len(results_df)} contrasts!")

    # Store results for AI gene analysis
    st.session_state['selected_contrasts_for_ai'] = \
        results_df[['analysis_id','contrast_id']].to_dict('records')
    st.session_state['research_query'] = research_query

    # Display results table
    st.subheader("Contrast Relevance Scores")

    # Configure column display
    display_columns = ['Accession', 'contrast_id', 'RelevanceScore', 'Description']
    if 'Run1Justification' in results_df.columns:
        display_columns.append('Run1Justification')

    st.dataframe(
        results_df[display_columns],
        use_container_width=True,
        column_config={
            "RelevanceScore": st.column_config.NumberColumn(
                "Relevance Score",
                format="%.2f",
                help="AI-assessed relevance score (0-1 scale)"
            ),
            "contrast_id": st.column_config.TextColumn(
                "Contrast",
                help="Contrast identifier"
            ),
            "Description": st.column_config.TextColumn(
                "Description",
                help="Contrast description"
            ),
            "Run1Justification": st.column_config.TextColumn(
                "AI Justification",
                help="AI explanation for the relevance score"
            ),
            "Accession": st.column_config.TextColumn(
                "Dataset",
                help="Dataset accession"
            )
        }
    )

    # Show summary statistics
    _display_relevance_summary(results_df)

    # Provide download option
    _provide_relevance_download(results_df)


@log_streamlit_function
def _display_relevance_summary(results_df):
    """Display summary statistics for relevance scores."""
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Highly Relevant (≥0.7)", len(results_df[results_df['RelevanceScore'] >= 0.7]))
    with col2:
        st.metric("Moderately Relevant (0.4-0.7)", len(results_df[(results_df['RelevanceScore'] >= 0.4) & (results_df['RelevanceScore'] < 0.7)]))
    with col3:
        st.metric("Low Relevance (<0.4)", len(results_df[results_df['RelevanceScore'] < 0.4]))


@log_streamlit_function
def _provide_relevance_download(results_df):
    """Provide download option for relevance results."""
    csv = results_df.to_csv(index=False)
    st.download_button(
        label="📥 Download Relevance Scores as CSV",
        data=csv,
        file_name=f"contrast_relevance_scores.csv",
        mime="text/csv"
    )





@log_streamlit_function
@log_streamlit_agent
def _execute_ai_analysis(agent, prompt: str) -> str:
    """Execute the AI analysis asynchronously."""
    loop = asyncio.new_event_loop()
    asyncio.set_event_loop(loop)

    async def run_analysis():
        async with agent.run_mcp_servers():
            result = await agent.run(prompt)
            return result.output if hasattr(result, 'output') else str(result)

    try:
        result_text = loop.run_until_complete(run_analysis())
        return result_text
    finally:
        loop.close()


@log_streamlit_function
def _display_ai_analysis_results(result_text: str, research_question: str, selected_contrasts: List[Dict]):
    """Display the AI gene analysis results."""
    log_streamlit_event("AI gene analysis completed successfully")
    st.success("✅ Analysis completed successfully!")

    st.subheader("📑 AI Gene Analysis Results")
    st.markdown(result_text)

    # Add download option
    analysis_report = f"""
# AI Gene Analysis Report

**Research Question:** {research_question}

**Contrasts Analyzed:** {len(selected_contrasts)}

## Analysis Results

{result_text}

---
*Generated by UORCA Explorer AI Assistant*
"""

    st.download_button(
        label="📄 Download Analysis Report",
        data=analysis_report,
        file_name="uorca_ai_gene_analysis.md",
        mime="text/markdown"
    )

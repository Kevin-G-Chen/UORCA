#!/usr/bin/env python3
"""
AI Key Findings Analysis for UORCA Results
==========================================

This module implements a comprehensive AI agent workflow that:
1. Takes contrast relevance assessment results
2. Uses MCP tools to select diverse, relevant contrasts
3. Analyzes genes across selected contrasts
4. Generates prose summaries of key findings

The agent uses structured tools to programmatically analyze RNA-seq data
and provide biological insights.
"""

import asyncio
import json
import logging
import pandas as pd
import streamlit as st
from typing import Dict, List, Any, Optional, Tuple
from pathlib import Path
import os

from ai_agent_factory import create_example_agent

logger = logging.getLogger(__name__)

async def analyze_key_findings(
    results_dir: str,
    contrast_relevance_df: pd.DataFrame,
    research_query: str,
    max_contrasts: int = 5,
    min_relevance_score: float = 0.6,
    p_value_thresh: float = 0.05,
    lfc_thresh: float = 1.0
) -> Dict[str, Any]:
    """
    Main workflow for AI-driven key findings analysis.

    Args:
        results_dir: Path to UORCA results directory
        contrast_relevance_df: DataFrame with contrast relevance scores
        research_query: Original research question
        max_contrasts: Maximum number of contrasts to analyze
        min_relevance_score: Minimum relevance score to consider
        p_value_thresh: P-value threshold for DEG filtering
        lfc_thresh: Log fold change threshold for DEG filtering

    Returns:
        Dictionary containing analysis results and AI response
    """
    try:
        # Initialize the agent
        agent = create_example_agent()

        # Prepare contrast data for the agent
        contrast_data = prepare_contrast_data(contrast_relevance_df, min_relevance_score)

        if not contrast_data:
            return {
                "success": False,
                "error": "No contrasts meet the minimum relevance criteria",
                "response": None
            }

        # Create the analysis prompt
        analysis_prompt = create_analysis_prompt(
            research_query=research_query,
            contrast_data=contrast_data,
            results_dir=results_dir,
            max_contrasts=max_contrasts,
            p_value_thresh=p_value_thresh,
            lfc_thresh=lfc_thresh
        )

        # Run the agent analysis
        async with agent.run_mcp_servers():
            result = await agent.run(analysis_prompt)
            response_text = result.output if hasattr(result, 'output') else str(result)

        return {
            "success": True,
            "response": response_text,
            "contrasts_available": len(contrast_data),
            "parameters": {
                "max_contrasts": max_contrasts,
                "min_relevance_score": min_relevance_score,
                "p_value_thresh": p_value_thresh,
                "lfc_thresh": lfc_thresh
            }
        }

    except Exception as e:
        logger.error(f"Error in key findings analysis: {str(e)}", exc_info=True)
        return {
            "success": False,
            "error": str(e),
            "response": None
        }

def prepare_contrast_data(
    contrast_relevance_df: pd.DataFrame,
    min_relevance_score: float
) -> List[Dict[str, Any]]:
    """
    Prepare contrast data for the AI agent.

    Args:
        contrast_relevance_df: DataFrame with contrast relevance scores
        min_relevance_score: Minimum score threshold

    Returns:
        List of contrast dictionaries for the agent
    """
    # Filter by minimum relevance score
    filtered_df = contrast_relevance_df[
        contrast_relevance_df['RelevanceScore'] >= min_relevance_score
    ].copy()

    if filtered_df.empty:
        return []

    # Sort by relevance score (descending)
    filtered_df = filtered_df.sort_values('RelevanceScore', ascending=False)

    # Convert to list of dictionaries
    contrast_data = []
    for _, row in filtered_df.iterrows():
        contrast_info = {
            "dataset_id": row.get('analysis_id', ''),
            "contrast_id": row.get('contrast_id', ''),
            "relevance_score": float(row.get('RelevanceScore', 0)),
            "description": row.get('description', ''),
            "accession": row.get('accession', ''),
            "organism": row.get('organism', 'Unknown')
        }

        # Add justification if available
        if 'Run1Justification' in row:
            contrast_info["ai_justification"] = row['Run1Justification']

        contrast_data.append(contrast_info)

    return contrast_data

def create_analysis_prompt(
    research_query: str,
    contrast_data: List[Dict[str, Any]],
    results_dir: str,
    max_contrasts: int,
    p_value_thresh: float,
    lfc_thresh: float
) -> str:
    """
    Create the analysis prompt for the AI agent.

    Args:
        research_query: Original research question
        contrast_data: List of contrast information
        results_dir: Path to results directory
        max_contrasts: Maximum contrasts to select
        p_value_thresh: P-value threshold
        lfc_thresh: Log fold change threshold

    Returns:
        Formatted prompt string
    """

    # Convert contrast data to JSON for the prompt
    contrasts_json = json.dumps(contrast_data, indent=2)

    prompt = f"""
I need you to analyze RNA-seq results and identify key findings. Here's the context:

RESEARCH QUESTION: "{research_query}"

RESULTS DIRECTORY: {results_dir}

AVAILABLE CONTRASTS (with relevance scores):
{contrasts_json}

TASK:
Please perform a comprehensive analysis to identify key findings by following these steps:

1. SELECT CONTRASTS:
   - Choose up to {max_contrasts} contrasts from the list above
   - Prioritize both HIGH RELEVANCE SCORES and DIVERSITY
   - Aim for diversity across organisms, experimental conditions, and biological contexts
   - Explain your selection rationale

2. ANALYZE GENES:
   - For each selected contrast, use filter_degs() with thresholds: p < {p_value_thresh}, |logFC| > {lfc_thresh}
   - Use count_gene_occurrences() to identify genes appearing across multiple contrasts
   - Focus on recurring patterns and notable individual findings
   - Use get_gene_stats() for specific genes of high interest

3. SYNTHESIZE FINDINGS:
   - Identify the most important recurring genes and their expression patterns
   - Note organism-specific vs cross-species patterns
   - Highlight consistent vs context-dependent gene behaviors
   - Provide biological interpretation where possible

4. GENERATE SUMMARY:
   Please provide a comprehensive prose response addressing "What are some of the key findings across these datasets?" Include:
   - Summary of selected contrasts and selection rationale
   - Top recurring genes with their expression patterns
   - Biological interpretation and significance
   - Notable context-specific observations
   - Any limitations or caveats

Please proceed with this analysis using the available MCP tools.
"""

    return prompt

def run_streamlit_key_findings_analysis(
    results_dir: str,
    contrast_relevance_df: pd.DataFrame,
    research_query: str
) -> None:
    """
    Streamlit interface for running key findings analysis.

    Args:
        results_dir: Path to UORCA results directory
        contrast_relevance_df: DataFrame with contrast relevance scores
        research_query: Original research question
    """

    st.subheader("🔍 AI Key Findings Analysis")
    st.markdown("**Generate comprehensive insights from your contrast relevance assessment using AI analysis.**")

    # Analysis parameters
    with st.expander("🔧 Analysis Parameters", expanded=False):
        col1, col2 = st.columns(2)

        with col1:
            max_contrasts = st.slider(
                "Maximum contrasts to analyze",
                min_value=2,
                max_value=10,
                value=5,
                help="AI will select up to this many contrasts for analysis"
            )

            min_relevance = st.slider(
                "Minimum relevance score",
                min_value=0.0,
                max_value=1.0,
                value=0.6,
                step=0.1,
                help="Only contrasts with relevance ≥ this threshold will be considered"
            )

        with col2:
            p_thresh = st.number_input(
                "P-value threshold",
                min_value=0.001,
                max_value=0.1,
                value=0.05,
                step=0.001,
                format="%.3f",
                help="P-value threshold for DEG filtering"
            )

            lfc_thresh = st.number_input(
                "Log fold change threshold",
                min_value=0.5,
                max_value=3.0,
                value=1.0,
                step=0.1,
                help="Absolute log fold change threshold for DEG filtering"
            )

    # Show preview of contrasts that will be considered
    eligible_contrasts = contrast_relevance_df[
        contrast_relevance_df['RelevanceScore'] >= min_relevance
    ]

    st.info(f"📊 {len(eligible_contrasts)} contrasts meet the minimum relevance threshold of {min_relevance}")

    if len(eligible_contrasts) < 2:
        st.warning("⚠️ At least 2 contrasts are needed for meaningful analysis. Consider lowering the minimum relevance threshold.")
        return

    # Analysis button
    if st.button("🚀 Generate Key Findings", type="primary"):
        if not os.getenv("OPENAI_API_KEY"):
            st.error("❌ OpenAI API key not found. Please set the OPENAI_API_KEY environment variable.")
            return

        # Run the analysis
        with st.spinner("🤖 AI is analyzing your data to identify key findings... This may take several minutes."):
            try:
                # Run the async analysis
                loop = asyncio.new_event_loop()
                asyncio.set_event_loop(loop)

                result = loop.run_until_complete(
                    analyze_key_findings(
                        results_dir=results_dir,
                        contrast_relevance_df=contrast_relevance_df,
                        research_query=research_query,
                        max_contrasts=max_contrasts,
                        min_relevance_score=min_relevance,
                        p_value_thresh=p_thresh,
                        lfc_thresh=lfc_thresh
                    )
                )

                loop.close()

                # Display results
                if result["success"]:
                    st.success("✅ Analysis completed successfully!")

                    # Show parameters used
                    with st.expander("📋 Analysis Parameters Used", expanded=False):
                        params = result["parameters"]
                        st.json(params)

                    # Show the AI response
                    st.subheader("🎯 Key Findings")
                    st.markdown(result["response"])

                    # Add download option
                    analysis_report = f"""
# Key Findings Analysis Report

**Research Question:** {research_query}

**Analysis Parameters:**
- Maximum contrasts analyzed: {result['parameters']['max_contrasts']}
- Minimum relevance score: {result['parameters']['min_relevance_score']}
- P-value threshold: {result['parameters']['p_value_thresh']}
- Log fold change threshold: {result['parameters']['lfc_thresh']}
- Contrasts available: {result['contrasts_available']}

## Key Findings

{result['response']}

---
*Generated by UORCA Explorer AI Assistant*
"""

                    st.download_button(
                        label="📄 Download Analysis Report",
                        data=analysis_report,
                        file_name="uorca_key_findings_report.md",
                        mime="text/markdown"
                    )

                else:
                    st.error(f"❌ Analysis failed: {result.get('error', 'Unknown error')}")

            except Exception as e:
                logger.error(f"Error in Streamlit key findings analysis: {str(e)}", exc_info=True)
                st.error(f"❌ Analysis failed: {str(e)}")

                with st.expander("🔍 Error Details", expanded=False):
                    import traceback
                    st.code(traceback.format_exc())

def validate_contrast_relevance_data(df: pd.DataFrame) -> Tuple[bool, str]:
    """
    Validate that the contrast relevance DataFrame has required columns.

    Args:
        df: Contrast relevance DataFrame

    Returns:
        Tuple of (is_valid, error_message)
    """
    required_columns = ['analysis_id', 'contrast_id', 'RelevanceScore']

    missing_columns = [col for col in required_columns if col not in df.columns]

    if missing_columns:
        return False, f"Missing required columns: {missing_columns}"

    if len(df) == 0:
        return False, "DataFrame is empty"

    # Check for valid relevance scores
    if not df['RelevanceScore'].between(0, 1).all():
        return False, "RelevanceScore values must be between 0 and 1"

    return True, ""

# Example usage function for testing
async def example_key_findings_analysis():
    """
    Example function showing how to use the key findings analysis.
    This is for testing and demonstration purposes.
    """
    # Mock data for testing
    mock_contrast_data = pd.DataFrame({
        'analysis_id': ['GSE111143', 'GSE123456', 'GSE789012'],
        'contrast_id': ['Day8_vs_Naive', 'Treatment_vs_Control', 'KO_vs_WT'],
        'RelevanceScore': [0.85, 0.72, 0.68],
        'description': [
            'Day 8 post-infection vs naive T cells',
            'Drug treatment vs vehicle control',
            'Knockout vs wildtype comparison'
        ],
        'accession': ['GSE111143', 'GSE123456', 'GSE789012'],
        'organism': ['Mus musculus', 'Homo sapiens', 'Mus musculus']
    })

    results_dir = "/path/to/results"
    research_query = "T cell activation and immune response"

    result = await analyze_key_findings(
        results_dir=results_dir,
        contrast_relevance_df=mock_contrast_data,
        research_query=research_query
    )

    print("Analysis Result:")
    print(json.dumps(result, indent=2))

if __name__ == "__main__":
    # Run example if called directly
    asyncio.run(example_key_findings_analysis())

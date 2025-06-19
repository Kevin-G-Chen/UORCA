"""
AI Assistant Tab for UORCA Explorer.

This tab provides AI-powered analysis and exploration capabilities.
"""

import os
import json
import asyncio
import logging
import streamlit as st
import pandas as pd
import pydantic
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
from ai_gene_schema import GeneAnalysisOutput

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
    from contrast_relevance import run_contrast_relevance, run_contrast_relevance_with_selection, SelectedContrast
    CONTRAST_RELEVANCE_AVAILABLE = True
    CONTRAST_RELEVANCE_WITH_SELECTION_AVAILABLE = True
except ImportError as e:
    logger.warning(f"Contrast relevance not available: {e}")
    CONTRAST_RELEVANCE_AVAILABLE = False
    CONTRAST_RELEVANCE_WITH_SELECTION_AVAILABLE = False


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

    # Step 1: Contrast Relevance Assessment with Selection
    with st.spinner("Step 1/2: Assessing contrast relevance and selecting optimal subset... This may take a few minutes."):
        try:
            # Try to use intelligent selection if available
            if CONTRAST_RELEVANCE_WITH_SELECTION_AVAILABLE:
                # Run contrast relevance assessment with intelligent selection
                results_df, selected_contrasts = run_contrast_relevance_with_selection(
                    ri,
                    query=research_query,
                    repeats=3,  # Reduced since we're doing more complex processing
                    batch_size=100,  # Larger batch size for better selection
                    parallel_jobs=4
                )

                if not results_df.empty and selected_contrasts:
                    _display_relevance_and_selection_results(
                        ri, results_df, selected_contrasts, research_query
                    )

                    # Store SELECTED contrasts for AI gene analysis (not all contrasts)
                    selected_contrast_dicts = [
                        {'analysis_id': sc.analysis_id, 'contrast_id': sc.contrast_id}
                        for sc in selected_contrasts
                    ]
                    st.session_state['selected_contrasts_for_ai'] = selected_contrast_dicts
                    st.session_state['research_query'] = research_query

                    # Show selection summary
                    st.success(f"✅ Selected {len(selected_contrasts)} optimal contrasts from {len(results_df)} total contrasts for detailed analysis!")

                else:
                    st.warning("No contrasts found for assessment.")
                    return
            else:
                # Fall back to original approach
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

    # Step 2: AI Gene Analysis (now using SELECTED contrasts only)
    st.markdown("---")
    with st.spinner("Step 2/2: AI is analyzing differential expression patterns in selected contrasts... This may take several minutes."):
        try:
            # Set the RESULTS_DIR environment variable for the MCP server
            os.environ['RESULTS_DIR'] = results_dir

            selected_contrast_dicts = st.session_state['selected_contrasts_for_ai']

            agent = create_uorca_agent()

            # Enhanced prompt that leverages the selection
            if CONTRAST_RELEVANCE_WITH_SELECTION_AVAILABLE and len(selected_contrast_dicts) <= 20:
                prompt = f"""
Research question: "{research_query}"

I have intelligently selected the following {len(selected_contrast_dicts)} contrasts from a larger set based on relevance, diversity, and analytical value:

{json.dumps(selected_contrast_dicts, indent=2)}

These contrasts were chosen to provide:
1. High relevance to the research question
2. Diversity of biological contexts
3. Appropriate controls and comparisons
4. Maximum analytical power for comparative analysis

Please perform comprehensive analysis using your four tools:
1. Use get_most_common_genes() to find genes frequently differentially expressed across these selected contrasts
2. Use filter_genes_by_contrast_sets() to find genes specific to subsets of contrasts (e.g., treatment-specific vs control-specific)
3. Use get_gene_contrast_stats() to drill into interesting genes
4. Use summarize_contrast() to understand individual contrast patterns

Choose reasonable thresholds and return:
1. Key genes identified with their biological significance
2. Patterns that distinguish different contrast categories
3. A brief biological interpretation of the findings
"""
            else:
                prompt = f"""
Research question: "{research_query}"

Available contrasts:
{json.dumps(selected_contrast_dicts, indent=2)}

Please perform the analysis using your four tools, choose all thresholds reasonably, and return:
1) A structured summary showing the key genes identified for each contrast or gene set
2) A brief 2-3 sentence biological interpretation explaining your rationale and what patterns you discovered.
"""

            # Run the agent analysis
            result_text = _execute_ai_analysis(agent, prompt)

            # Display results
            _display_ai_analysis_results(result_text, research_query, selected_contrast_dicts)

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

    # Add Selected column based on session state
    sel = st.session_state.get("selected_contrasts_for_ai", [])
    sel_set = {(item['analysis_id'], item['contrast_id']) for item in sel}
    results_df = results_df.copy()
    results_df["Selected"] = results_df.apply(
        lambda r: (r['analysis_id'], r['contrast_id']) in sel_set, axis=1
    )

    # Configure column display
    display_columns = ['Selected', 'Accession', 'contrast_id', 'RelevanceScore', 'Description']
    if 'Run1Justification' in results_df.columns:
        display_columns.append('Run1Justification')

    st.dataframe(
        results_df[display_columns],
        use_container_width=True,
        column_config={
            "Selected": st.column_config.CheckboxColumn("Selected"),
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

    # Provide download option
    _provide_relevance_download(results_df)





@log_streamlit_function
def _display_relevance_and_selection_results(
    ri: ResultsIntegrator,
    results_df: pd.DataFrame,
    selected_contrasts: List[SelectedContrast],
    research_query: str
):
    """Display both the full relevance results and the intelligent selection."""

    # Sort by relevance score
    results_df = results_df.sort_values('RelevanceScore', ascending=False)

    # Add contrast descriptions for display
    results_df['Description'] = results_df.apply(
        lambda row: ri._get_contrast_description(row['analysis_id'], row['contrast_id']),
        axis=1
    )

    # Add accession and organism info
    results_df['Accession'] = results_df['analysis_id'].map(
        lambda aid: ri.analysis_info.get(aid, {}).get('accession', aid)
    )
    results_df['organism'] = results_df['analysis_id'].map(
        lambda aid: ri.analysis_info.get(aid, {}).get('organism', 'Unknown')
    )

    log_streamlit_event(f"Contrast relevance with selection completed: {len(results_df)} assessed, {len(selected_contrasts)} selected")

    # Create tabs for full results vs selected subset
    full_tab, selected_tab = st.tabs(["📊 All Contrasts (Relevance Scores)", "🎯 Selected Contrasts (AI Chosen)"])

    with full_tab:
        st.subheader("All Contrast Relevance Scores")
        st.markdown("*Complete relevance assessment for all available contrasts*")

        # Mark which contrasts were selected
        selected_pairs = {(sc.analysis_id, sc.contrast_id) for sc in selected_contrasts}
        results_df['AI_Selected'] = results_df.apply(
            lambda r: (r['analysis_id'], r['contrast_id']) in selected_pairs, axis=1
        )

        # Configure column display for full results
        display_columns = ['AI_Selected', 'Accession', 'contrast_id', 'RelevanceScore', 'Description']
        if 'Run1Justification' in results_df.columns:
            display_columns.append('Run1Justification')

        st.dataframe(
            results_df[display_columns],
            use_container_width=True,
            column_config={
                "AI_Selected": st.column_config.CheckboxColumn("AI Selected", help="Selected by AI for analysis"),
                "RelevanceScore": st.column_config.NumberColumn(
                    "Relevance Score",
                    format="%.2f",
                    help="AI-assessed relevance score (0-1 scale)"
                ),
                "contrast_id": st.column_config.TextColumn("Contrast", help="Contrast identifier"),
                "Description": st.column_config.TextColumn("Description", help="Contrast description"),
                "Run1Justification": st.column_config.TextColumn("AI Justification", help="AI explanation for the relevance score"),
                "Accession": st.column_config.TextColumn("Dataset", help="Dataset accession")
            }
        )

    with selected_tab:
        st.subheader("AI-Selected Contrasts for Analysis")
        st.markdown("*Intelligently chosen subset optimizing relevance, diversity, and analytical power*")

        # Create DataFrame for selected contrasts
        selected_data = []
        for sc in selected_contrasts:
            # Find the relevance info
            relevance_row = results_df[
                (results_df['analysis_id'] == sc.analysis_id) &
                (results_df['contrast_id'] == sc.contrast_id)
            ]

            if not relevance_row.empty:
                selected_data.append({
                    'Accession': relevance_row.iloc[0]['Accession'],
                    'Contrast': sc.contrast_id,
                    'Relevance Score': f"{sc.RelevanceScore:.2f}",
                    'Category': sc.category.category,
                    'Category Reason': sc.category.justification,
                    'Selection Reason': sc.selection_justification,
                    'Description': relevance_row.iloc[0]['Description']
                })

        if selected_data:
            selected_df = pd.DataFrame(selected_data)

            # Display with better formatting
            st.dataframe(
                selected_df,
                use_container_width=True,
                column_config={
                    "Relevance Score": st.column_config.TextColumn("Relevance", help="AI relevance score"),
                    "Category": st.column_config.TextColumn("Role", help="Analytical role in the study"),
                    "Category Reason": st.column_config.TextColumn("Role Justification", help="Why assigned this role"),
                    "Selection Reason": st.column_config.TextColumn("Selection Justification", help="Why selected for analysis"),
                    "Description": st.column_config.TextColumn("Description", help="Contrast description")
                }
            )

            # Show selection strategy
            st.info("🎯 **Selection Strategy**: The AI selected these contrasts to maximize analytical power while ensuring diversity and appropriate controls for comparative analysis.")

        # Category breakdown
        if selected_contrasts:
            st.markdown("### Selection Summary")
            category_counts = {}
            for sc in selected_contrasts:
                cat = sc.category.category
                category_counts[cat] = category_counts.get(cat, 0) + 1

            cols = st.columns(len(category_counts))
            for i, (category, count) in enumerate(category_counts.items()):
                with cols[i]:
                    st.metric(f"{category.title()} Contrasts", count)

    # Provide download option for both
    _provide_relevance_download(results_df)


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
    """Execute the AI analysis asynchronously with proper loop hygiene."""
    async def run_analysis():
        async with agent.run_mcp_servers():
            result = await agent.run(prompt)
            return result.output if hasattr(result, 'output') else str(result)

    try:
        # Check if there's already a running event loop
        try:
            loop = asyncio.get_running_loop()
            # If we're in an existing loop, we need to run in a thread
            import concurrent.futures
            with concurrent.futures.ThreadPoolExecutor() as executor:
                future = executor.submit(asyncio.run, run_analysis())
                result_text = future.result()
        except RuntimeError:
            # No running loop, safe to use asyncio.run
            result_text = asyncio.run(run_analysis())

        return result_text
    except Exception as e:
        logger.error(f"Error in AI analysis execution: {str(e)}", exc_info=True)
        raise


@log_streamlit_function
def _display_ai_analysis_results(result_text: str, research_question: str, selected_contrasts: List[Dict]):
    """Display the AI gene analysis results with structured output."""
    log_streamlit_event("AI gene analysis completed successfully")
    st.success("✅ Analysis completed successfully!")

    # Parse structured JSON output
    try:
        parsed = result_text
    except pydantic.ValidationError as e:
        st.error(f"Agent output failed validation: {e}")
        st.subheader("📑 Raw AI Response")
        st.markdown(result_text)
        return

    # Get ResultsIntegrator from session state or create one
    ri = st.session_state.get('results_integrator')
    if not ri:
        try:
            from .helpers import get_integrator
            results_dir = st.session_state.get('results_dir', '/UORCA_results')
            ri, _ = get_integrator(results_dir)
        except ImportError:
            st.warning("Could not load data integrator. Some visualizations may not be available.")
            ri = None

    genes = set(parsed.genes)  # normalize + dedupe
    contrasts = [(item['analysis_id'], item['contrast_id']) for item in selected_contrasts]

    # --- 1. Gene list
    st.subheader("🧬 Genes Selected by AI")
    st.write(f"**Filters used:** Log2FC ≥ {parsed.filters.lfc_thresh}, P-value < {parsed.filters.p_thresh}")

    # Display genes in a nice format
    gene_cols = st.columns(min(len(genes), 4))
    for i, gene in enumerate(sorted(genes)):
        with gene_cols[i % len(gene_cols)]:
            st.code(gene)

    # --- 2. LFC table
    st.subheader("📊 Log Fold Change Table")

    if ri and ri.deg_data:
        rows, missing = [], []
        for gene in genes:
            hit = False
            for aid, cid in contrasts:
                if aid in ri.deg_data and cid in ri.deg_data[aid]:
                    df = ri.deg_data[aid][cid]
                    if 'Gene' in df.columns and gene in df['Gene'].values:
                        hit = True
                        rec = df.loc[df['Gene'] == gene].iloc[0]
                        rows.append({
                            'Gene': gene,
                            'Dataset': aid,
                            'Contrast': cid,
                            'logFC': rec.get("logFC", 0),
                            'P-value': rec.get("adj.P.Val") or rec.get("P.Value", 1)
                        })
            if not hit:
                missing.append(gene)

        if rows:
            lfc_df = pd.DataFrame(rows)
            st.dataframe(
                lfc_df,
                use_container_width=True,
                column_config={
                    "logFC": st.column_config.NumberColumn("Log2FC", format="%.2f"),
                    "P-value": st.column_config.NumberColumn("P-value", format="%.2e")
                }
            )
        else:
            st.info("No expression data found for the selected genes in the chosen contrasts.")

        if missing:
            st.info(f"ℹ️ {len(missing)} gene(s) not present in the selected contrasts: {', '.join(missing)}")

        # --- 3. Heatmap
        if rows and len(genes - set(missing)) > 0:
            st.subheader("🌡️ Expression Heatmap")
            try:
                fig = ri.create_lfc_heatmap(
                    genes=list(genes - set(missing)),
                    contrasts=contrasts
                )
                if fig:
                    st.plotly_chart(fig, use_container_width=True)
                else:
                    st.info("Could not generate heatmap for the selected genes.")
            except Exception as e:
                st.warning(f"Error generating heatmap: {str(e)}")
    else:
        st.warning("No expression data available for visualization.")

    # --- 4. Interpretation
    st.subheader("🧠 AI-Driven Interpretation")
    st.markdown(parsed.interpretation)

    # Add download option
    analysis_report = f"""
# AI Gene Analysis Report

**Research Question:** {research_question}

**Contrasts Analyzed:** {len(selected_contrasts)}

## Selected Genes
{', '.join(sorted(genes))}

## Filter Parameters
- Log2 Fold Change Threshold: {parsed.filters.lfc_thresh}
- P-value Threshold: {parsed.filters.p_thresh}

## AI Interpretation

{parsed.interpretation}

---
*Generated by UORCA Explorer AI Assistant*
"""

    st.download_button(
        label="📄 Download Analysis Report",
        data=analysis_report,
        file_name="uorca_ai_gene_analysis.md",
        mime="text/markdown"
    )

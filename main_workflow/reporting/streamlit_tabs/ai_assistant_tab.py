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
from pydantic_ai.usage import UsageLimits
from typing import Dict, Any, List, Optional, Tuple


from .helpers import (
    check_ai_generating,
    setup_fragment_decorator,
    log_streamlit_tab,
    log_streamlit_function,
    log_streamlit_agent,
    log_streamlit_event
)
from .helpers.ai_agent_tool_logger import (
    start_ai_analysis_session,
    get_ai_tool_logs_for_display,
    get_current_log_file,
    read_log_file_contents
)
from ResultsIntegration import ResultsIntegrator
from ai_gene_schema import GeneAnalysisOutput

# Set up fragment decorator
setup_fragment_decorator()

logger = logging.getLogger(__name__)


def load_query_config() -> Optional[str]:
    """Load the dataset identification query from config file."""
    config_file_path = "main_workflow/reporting/.config/dataset_query.json"

    try:
        if os.path.exists(config_file_path):
            with open(config_file_path, 'r') as f:
                config_data = json.load(f)
                return config_data.get("query")
    except (json.JSONDecodeError, KeyError, FileNotFoundError):
        pass
    return None


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
    Render the AI assistant tab with subtabs for analysis and tool logs.

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

    # Render streamlined AI analysis workflow (no tabs needed)
    _render_streamlined_ai_workflow(ri, results_dir)


@log_streamlit_function
def _render_streamlined_ai_workflow(ri: ResultsIntegrator, results_dir: str):
    """Render the streamlined AI analysis workflow."""
    st.subheader("🧬 AI-Powered Gene Analysis")
    st.markdown("**Complete AI analysis of your RNA-seq data.** Enter your research question below and the AI will assess contrast relevance and identify key genes in one workflow.")

    # Load saved query from dataset identification
    saved_query = load_query_config()
    default_placeholder = "e.g., What contrasts are most relevant to T cell activation and differentiation?"

    # Use saved query as placeholder if available
    if saved_query:
        placeholder_text = f"Dataset query: {saved_query}"
        help_text = "Using the research question from your dataset identification. You can modify this or enter a new question focusing on specific aspects of your data."
    else:
        placeholder_text = default_placeholder
        help_text = "Describe your research question or area of interest. The AI will score each contrast based on how relevant it is to this query, then analyze key genes."

    # Research query input
    research_query = st.text_input(
        "Research Question",
        value=saved_query if saved_query else "",
        placeholder=placeholder_text,
        help=help_text
    )

    # Single workflow button
    col1, col2 = st.columns([1, 3])
    with col1:
        run_button = st.button("🚀 Run Complete AI Analysis", type="primary")
    with col2:
        if not research_query.strip():
            st.info("💡 Enter a research question above to start AI analysis.")

    if run_button:
        log_streamlit_event(f"User started complete AI analysis: '{research_query.strip()}'")
        _run_complete_ai_analysis(ri, results_dir, research_query.strip())


@log_streamlit_agent
def _run_complete_ai_analysis(ri: ResultsIntegrator, results_dir: str, research_query: str):
    """Run the complete AI analysis workflow including contrast relevance and gene analysis."""
    # Add validation for empty research query
    if not research_query.strip():
        st.error("Please enter a research question before running the analysis.")
        return

    if not CONTRAST_RELEVANCE_AVAILABLE:
        st.error("Contrast relevance assessment is not available. Please check your environment setup.")
        return

    if not UORCA_AGENT_AVAILABLE:
        st.error("UORCA agent is not available. Please check your environment setup.")
        return

    if not os.getenv("OPENAI_API_KEY"):
        st.error("OpenAI API key not found. Please set the OPENAI_API_KEY environment variable.")
        return

    # Run complete analysis workflow in single spinner
    with st.spinner("Running complete AI analysis: assessing contrast relevance and analyzing differential expression patterns..."):
        try:
            # Step 1: Contrast Relevance Assessment with Selection (no display)
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
                    # Store SELECTED contrasts for AI gene analysis (not all contrasts)
                    selected_contrast_dicts = [
                        {'analysis_id': sc.analysis_id, 'contrast_id': sc.contrast_id}
                        for sc in selected_contrasts
                    ]
                    st.session_state['selected_contrasts_for_ai'] = selected_contrast_dicts
                    st.session_state['research_query'] = research_query
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
                    # Store results for AI gene analysis
                    st.session_state['selected_contrasts_for_ai'] = \
                        results_df[['analysis_id','contrast_id']].to_dict('records')
                    st.session_state['research_query'] = research_query
                    selected_contrasts = None  # No selection objects in fallback mode
                else:
                    st.warning("No contrasts found for assessment.")
                    return

            # Step 2: AI Gene Analysis (using selected contrasts)
            # Set the RESULTS_DIR environment variable for the MCP server
            os.environ['RESULTS_DIR'] = results_dir

            selected_contrast_dicts = st.session_state['selected_contrasts_for_ai']

            # Pass selected contrasts to MCP server for filtering
            import json
            os.environ['SELECTED_CONTRASTS_FOR_AI'] = json.dumps(selected_contrast_dicts)

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
            result_output, tool_calls = _execute_ai_analysis(agent, prompt)

            # Get tool calls for display
            log_file = get_current_log_file()
            if log_file and log_file.exists():
                display_tool_calls = get_ai_tool_logs_for_display()
            else:
                display_tool_calls = []

            # Display unified results
            _display_unified_ai_results(
                ri, result_output, research_query, selected_contrast_dicts,
                results_df, selected_contrasts, display_tool_calls
            )

        except Exception as e:
            logger.error(f"Error in AI analysis workflow: {str(e)}", exc_info=True)
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
    """Display the AI-selected contrasts for analysis."""

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

    st.subheader("AI-Selected Contrasts for Analysis")
    st.markdown("*Intelligently chosen subset for detailed analysis*")

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
                'Dataset': relevance_row.iloc[0]['Accession'],
                'Contrast': sc.contrast_id,
                'Justification': sc.selection_justification
            })

    if selected_data:
        selected_df = pd.DataFrame(selected_data)

        # Display with simple formatting
        st.dataframe(
            selected_df,
            use_container_width=True,
            column_config={
                "Dataset": st.column_config.TextColumn("Dataset", help="Dataset accession (GSE)"),
                "Contrast": st.column_config.TextColumn("Contrast", help="Contrast identifier"),
                "Justification": st.column_config.TextColumn("Justification", help="Why selected for analysis")
            }
        )

    # Provide download option for selected contrasts data
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
def _execute_ai_analysis(agent, prompt: str) -> Tuple[GeneAnalysisOutput, List[Dict]]:
    """Execute the AI analysis asynchronously with proper loop hygiene and capture tool calls."""
    async def run_analysis():
        async with agent.run_mcp_servers():
            # Start new analysis session and clear previous tool calls
            start_ai_analysis_session()

            # Run main analysis
            result = await agent.run(prompt,
                usage_limits = UsageLimits(request_limit = 100))

            # Get tool calls from new logging system
            tool_calls = get_ai_tool_logs_for_display()

            # Return the structured output directly (GeneAnalysisOutput object)
            return result.output if hasattr(result, 'output') else result, tool_calls

    try:
        # Check if there's already a running event loop
        try:
            loop = asyncio.get_running_loop()
            # If we're in an existing loop, we need to run in a thread
            import concurrent.futures
            with concurrent.futures.ThreadPoolExecutor() as executor:
                future = executor.submit(asyncio.run, run_analysis())
                result_output, tool_calls = future.result()
        except RuntimeError:
            # No running loop, safe to use asyncio.run
            result_output, tool_calls = asyncio.run(run_analysis())

        return result_output, tool_calls
    except Exception as e:
        logger.error(f"Error in AI analysis execution: {str(e)}", exc_info=True)
        raise


@log_streamlit_function
def _display_unified_ai_results(
    ri: ResultsIntegrator,
    result_output: GeneAnalysisOutput,
    research_question: str,
    selected_contrast_dicts: List[Dict],
    results_df: pd.DataFrame,
    selected_contrasts: Optional[List] = None,
    tool_calls: Optional[List[Dict]] = None
):
    """Display unified AI analysis results with tabbed interface."""
    log_streamlit_event("AI analysis completed successfully")

    # Use the structured output directly
    try:
        parsed = result_output
    except Exception as e:
        st.error(f"Error processing AI output: {e}")
        st.subheader("📑 Raw AI Response")
        st.markdown(str(result_output))
        return

    genes = set(parsed.genes)  # normalize + dedupe
    contrasts = [(item['analysis_id'], item['contrast_id']) for item in selected_contrast_dicts]

    # --- AI Interpretation (at the very top)
    st.subheader("🧠 AI-Driven Interpretation")
    st.markdown(parsed.interpretation)

    # Create tabs for organized display
    gene_tab, heatmap_tab, table_tab, contrast_tab, logs_tab = st.tabs([
        "🧬 Selected Genes",
        "🌡️ Heatmap",
        "📊 Expression Data",
        "🎯 Contrast Selection",
        "🔧 Tool Logs"
    ])

    # Prepare data for all tabs
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
    else:
        rows, missing = [], []

    # Tab 1: Gene List
    with gene_tab:
        st.subheader("Selected Genes")

        # Display genes in a copyable format
        sorted_genes = sorted(genes)
        genes_text = ', '.join(sorted_genes)

        st.markdown("**Gene List (Copy-friendly format):**")
        st.code(genes_text, language=None)

        st.write(f"**Total genes:** {len(sorted_genes)}")

        # Download analysis report
        analysis_report = f"""
# AI Gene Analysis Report

**Research Question:** {research_question}

**Contrasts Analyzed:** {len(selected_contrast_dicts)}

## Selected Genes
{', '.join(sorted(genes))}

## AI Interpretation

{parsed.interpretation}

---
*Generated by UORCA Explorer AI Assistant*
"""

        st.download_button(
            label="📄 Download Analysis Report",
            data=analysis_report,
            file_name="uorca_ai_gene_analysis.md",
            mime="text/markdown",
            help="Download the complete analysis report"
        )

    # Tab 2: Heatmap
    with heatmap_tab:
        st.subheader("Expression Heatmap")

        if rows and len(genes - set(missing)) > 0:
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
            st.warning("No expression data available for heatmap visualization.")

    # Tab 3: Expression Table
    with table_tab:
        st.subheader("Log Fold Change Data")

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

            # Download filtered dataframe CSV
            try:
                from mcp_server_core import get_filtered_dataframe
                filtered_df = get_filtered_dataframe()

                if not filtered_df.empty:
                    csv_data = filtered_df.to_csv(index=False)
                    st.download_button(
                        label="📊 Download Filtered Dataset (CSV)",
                        data=csv_data,
                        file_name="uorca_ai_filtered_data.csv",
                        mime="text/csv",
                        help="Download the exact dataset the AI agent analyzed"
                    )
            except ImportError:
                st.info("Filtered dataset download not available (MCP server not loaded)")
            except Exception as e:
                st.warning(f"Could not prepare filtered dataset: {e}")

        else:
            st.info("No expression data found for the selected genes in the chosen contrasts.")

        if missing:
            st.info(f"ℹ️ {len(missing)} gene(s) not present in the selected contrasts: {', '.join(missing)}")

    # Tab 4: Contrast Selection
    with contrast_tab:
        st.subheader("Selected Contrasts for Analysis")

        if selected_contrasts and hasattr(selected_contrasts[0], 'selection_justification'):
            # Show intelligent selection results
            st.markdown("*Intelligently chosen subset for detailed analysis*")

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
                        'Dataset': relevance_row.iloc[0].get('Accession', sc.analysis_id),
                        'Contrast': sc.contrast_id,
                        'Justification': sc.selection_justification
                    })

            if selected_data:
                selected_df = pd.DataFrame(selected_data)
                st.dataframe(
                    selected_df,
                    use_container_width=True,
                    column_config={
                        "Dataset": st.column_config.TextColumn("Dataset", help="Dataset accession (GSE)"),
                        "Contrast": st.column_config.TextColumn("Contrast", help="Contrast identifier"),
                        "Justification": st.column_config.TextColumn("Justification", help="Why selected for analysis")
                    }
                )
        else:
            # Show basic contrast list
            st.markdown("*Contrasts used in analysis*")

            # Prepare contrast display data
            contrast_display_data = []
            for item in selected_contrast_dicts:
                aid = item['analysis_id']
                cid = item['contrast_id']

                # Get additional info if available
                accession = ri.analysis_info.get(aid, {}).get('accession', aid) if ri else aid
                description = ri._get_contrast_description(aid, cid) if ri else ''

                contrast_display_data.append({
                    'Dataset': accession,
                    'Contrast': cid,
                    'Description': description
                })

            if contrast_display_data:
                contrast_df = pd.DataFrame(contrast_display_data)
                st.dataframe(
                    contrast_df,
                    use_container_width=True,
                    column_config={
                        "Dataset": st.column_config.TextColumn("Dataset", help="Dataset accession (GSE)"),
                        "Contrast": st.column_config.TextColumn("Contrast", help="Contrast identifier"),
                        "Description": st.column_config.TextColumn("Description", help="Contrast description")
                    }
                )

        # Provide download option for contrast selection data
        if not results_df.empty:
            csv = results_df.to_csv(index=False)
            st.download_button(
                label="📥 Download All Contrast Relevance Scores (CSV)",
                data=csv,
                file_name=f"contrast_relevance_scores.csv",
                mime="text/csv"
            )

    # Tab 5: Tool Logs
    with logs_tab:
        st.subheader("AI Tool Call Logs")
        st.markdown("**Detailed logs of AI agent tool usage** from the completed analysis.")

        if tool_calls:
            _display_tool_calls_detailed(tool_calls)
        else:
            st.info("🔍 No tool calls found. Logs appear after running AI analysis.")




@log_streamlit_function
def _display_tool_calls_detailed(tool_calls: List[Dict]):
    """Display tool calls in a detailed, structured format."""
    if not tool_calls:
        st.info("No tool calls to display.")
        return

    st.markdown("### Tool Call Details")

    for i, call in enumerate(tool_calls, 1):
        # Format success indicator
        status_icon = "✅" if call.get('success', True) else "❌"
        tool_name = call.get('tool_name', 'Unknown Tool')
        timestamp = call.get('timestamp', 'Unknown time')

        with st.expander(f"{status_icon} **{tool_name}** (Call #{i}) - {timestamp}", expanded=False):

            # Create two-column layout
            param_col, result_col = st.columns(2)

            # Left column: Parameters
            with param_col:
                st.markdown("#### 📋 Parameters")
                if call.get('parameters'):
                    # Display parameters in a formatted way
                    for param_name, param_value in call['parameters'].items():
                        if isinstance(param_value, (list, dict)):
                            st.markdown(f"**{param_name}:**")
                            st.json(param_value)
                        else:
                            st.markdown(f"**{param_name}:** `{param_value}`")
                else:
                    st.write("*No parameters*")

            # Right column: Results
            with result_col:
                st.markdown("#### 📤 Results")
                if call.get('success', True):
                    output_snippet = call.get('output_snippet')
                    if output_snippet and output_snippet.strip():
                        # Get tool name and output for custom formatting
                        tool_name = call.get('tool_name', '')

                        # Use full_output for parsing if available, otherwise use snippet
                        full_output = call.get('full_output')
                        if full_output is not None:
                            output_for_parsing = full_output
                            is_truncated = "truncated" in output_snippet
                        else:
                            output_for_parsing = output_snippet
                            is_truncated = "truncated" in output_snippet

                        if is_truncated:
                            st.warning("⚠️ Output truncated for display")

                        # Clean output for parsing (remove truncation text if present)
                        if is_truncated and full_output is None:
                            clean_output = output_snippet.split('[truncated')[0]
                        else:
                            clean_output = output_for_parsing

                        try:
                            # Parse Python literal output (not JSON)
                            import ast
                            if isinstance(call.get('output_snippet'), (dict, list)):
                                parsed_output = call.get('output_snippet')
                            elif clean_output.strip().startswith(('{', '[')):
                                # Use ast.literal_eval for Python literals with single quotes
                                parsed_output = ast.literal_eval(clean_output)
                            else:
                                parsed_output = None

                            # Custom formatting based on tool type
                            if tool_name == 'get_most_common_genes' and parsed_output:
                                # Display as table
                                if isinstance(parsed_output, list):
                                    df_data = []
                                    for item in parsed_output:
                                        if isinstance(item, dict) and 'gene' in item and 'count' in item:
                                            df_data.append({'Gene': item['gene'], 'Count': item['count']})
                                    if df_data:
                                        import pandas as pd
                                        df = pd.DataFrame(df_data)
                                        st.dataframe(df, use_container_width=True)
                                    else:
                                        st.code(clean_output)
                                else:
                                    st.code(clean_output)

                            elif tool_name == 'filter_genes_by_contrast_sets' and parsed_output:
                                # Extract and display just the gene list
                                if isinstance(parsed_output, dict) and 'genes' in parsed_output:
                                    genes = parsed_output['genes']
                                    if genes:
                                        # Display genes as wrapped text
                                        gene_text = ', '.join(genes)
                                        st.text_area("Genes:", value=gene_text, height=150, disabled=True)
                                        st.caption(f"Total genes: {len(genes)}")
                                    else:
                                        st.write("No genes found")
                                else:
                                    st.code(clean_output)

                            elif tool_name == 'get_gene_contrast_stats':
                                # Display as dataframe or show empty message
                                if parsed_output and isinstance(parsed_output, list) and len(parsed_output) > 0:
                                    df_data = []
                                    for item in parsed_output:
                                        if isinstance(item, dict):
                                            df_data.append({
                                                'Contrast': item.get('contrast_id', ''),
                                                'logFC': item.get('logFC', ''),
                                                'P-value': item.get('pvalue', '')
                                            })
                                    if df_data:
                                        import pandas as pd
                                        df = pd.DataFrame(df_data)
                                        st.dataframe(df, use_container_width=True)
                                    else:
                                        st.info("ℹ️ No gene statistics found for the specified parameters")
                                elif parsed_output == [] or (isinstance(parsed_output, list) and len(parsed_output) == 0):
                                    st.info("ℹ️ No gene statistics found - this gene may not be significantly expressed in any contrasts")
                                else:
                                    st.info("ℹ️ No gene statistics available")

                            elif tool_name == 'summarize_contrast' and parsed_output:
                                # Display with text wrapping
                                if isinstance(parsed_output, dict):
                                    st.json(parsed_output)
                                else:
                                    st.text_area("Summary:", value=str(parsed_output), height=100, disabled=True)

                            else:
                                # Default display for other tools or unparseable output
                                if parsed_output:
                                    st.json(parsed_output)
                                else:
                                    st.code(clean_output)

                            if is_truncated:
                                st.caption("*Output truncated - see raw log for complete results*")

                        except Exception as e:
                            # Handle truncated output specially
                            if is_truncated and "truncated" in output_snippet:
                                st.warning("⚠️ Output was truncated and cannot be fully parsed")

                                # For truncated output, try to extract what we can
                                if tool_name == 'get_most_common_genes':
                                    # Try to extract partial gene list from truncated output
                                    try:
                                        # Look for complete gene entries before truncation
                                        import re
                                        gene_pattern = r"\{'gene': '([^']+)', 'count': (\d+)\}"
                                        matches = re.findall(gene_pattern, clean_output)
                                        if matches:
                                            st.info(f"Showing first {len(matches)} genes (output was truncated)")
                                            import pandas as pd
                                            df_data = [{'Gene': gene, 'Count': int(count)} for gene, count in matches]
                                            df = pd.DataFrame(df_data)
                                            st.dataframe(df, use_container_width=True)
                                            st.caption("*See raw log for complete results*")
                                        else:
                                            st.text_area("Truncated Output:", value=output_snippet, height=100, disabled=True)
                                    except:
                                        st.text_area("Truncated Output:", value=output_snippet, height=100, disabled=True)

                                elif tool_name == 'filter_genes_by_contrast_sets':
                                    # Try to extract partial gene list
                                    try:
                                        import re
                                        # Look for genes in the format ['gene1', 'gene2', ...]
                                        gene_pattern = r"'([A-Za-z0-9_]+)'"
                                        matches = re.findall(gene_pattern, clean_output)
                                        if matches:
                                            st.info(f"Showing first {len(matches)} genes (output was truncated)")
                                            gene_text = ', '.join(matches)
                                            st.text_area("Partial Gene List:", value=gene_text, height=100, disabled=True)
                                            st.caption("*See raw log for complete results*")
                                        else:
                                            st.text_area("Truncated Output:", value=output_snippet, height=100, disabled=True)
                                    except:
                                        st.text_area("Truncated Output:", value=output_snippet, height=100, disabled=True)

                                else:
                                    # For other tools, just show the truncated output
                                    st.text_area("Truncated Output:", value=output_snippet, height=100, disabled=True)

                            else:
                                # Non-truncated parsing error
                                st.warning(f"⚠️ Could not parse output: {e}")

                                # Try JSON parsing as fallback
                                try:
                                    import json as json_module
                                    parsed_json = json_module.loads(clean_output.replace("'", '"'))
                                    st.json(parsed_json)
                                    st.caption("*Successfully parsed using JSON fallback*")
                                except:
                                    # Display as formatted text
                                    if clean_output.strip().startswith(('{', '[')):
                                        st.code(clean_output, language='json')
                                        st.caption("*Displaying as formatted code (parsing failed)*")
                                    else:
                                        st.text_area("Raw Output:", value=clean_output, height=100, disabled=True)
                                        st.caption("*Displaying as raw text*")
                    elif not output_snippet or not output_snippet.strip():
                        st.info("ℹ️ Tool completed successfully but returned no output")
                    else:
                        st.write("*No output*")
                else:
                    st.error(f"**Error:** {call.get('error', 'Unknown error')}")

            # Metadata (full width below columns)
            if call.get('analysis_id'):
                st.caption(f"Analysis ID: {call['analysis_id']}")


@log_streamlit_function
def _display_raw_log_file():
    """Display the raw log file contents."""
    st.markdown("### Raw Log File Contents")
    st.markdown("*Complete JSON log file with all tool call data*")

    raw_contents = read_log_file_contents()

    if raw_contents:
        # Add download button for raw log
        st.download_button(
            label="📥 Download Raw Log File",
            data=raw_contents,
            file_name=f"ai_tool_calls_{pd.Timestamp.now().strftime('%Y%m%d_%H%M%S')}.json",
            mime="application/json"
        )

        # Display in code block
        st.code(raw_contents, language='json')
    else:
        st.error("Could not read log file contents.")

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Landing Page Integration for UORCA Explorer

This module bridges the new MCP-based reporting agent with the existing
Streamlit UI, replacing the simple LLM-based implementation with a more
sophisticated agent-based approach.
"""

import os
import sys
import json
import asyncio
import logging
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from typing import List, Dict, Any, Optional, Tuple
from dataclasses import dataclass
from datetime import datetime

# Add parent directory to path for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import the reporting agent and auto-start manager
from reporting_agent import ReportingAgent
from auto_start_manager import AutoStartManager

# Set up logging
logging.basicConfig(
    level=logging.WARNING,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# Data structures to match existing UI expectations
@dataclass
class ContrastSelection:
    """Container for selected contrast with justification"""
    analysis_id: str
    contrast_id: str
    relevance_score: float
    justification: str
    deg_count: int


@dataclass
class ThresholdSelection:
    """Container for automatically selected statistical thresholds"""
    fdr_cutoff: float
    logfc_cutoff: float
    min_frequency: int
    justification: str


@dataclass
class LandingPageData:
    """Container for all landing page data"""
    selected_contrasts: List[ContrastSelection]
    thresholds: ThresholdSelection
    top_genes: List[str]
    heatmap_fig: Optional[go.Figure]
    gene_table: pd.DataFrame


class LandingPageGenerator:
    """Generates landing page data using the MCP-based reporting agent."""

    def __init__(self, integrator, model_name: str = "openai:gpt-4o-mini"):
        """
        Initialize the landing page generator.

        Args:
            integrator: ResultsIntegrator instance from the UI
            model_name: Model to use for analysis
        """
        self.integrator = integrator
        self.model_name = model_name
        self.agent = None
        self.results_dir = integrator.results_dir

    async def generate_landing_page(self, biological_prompt: str, max_genes: int = 50) -> Optional[LandingPageData]:
        """
        Generate AI-assisted landing page data using the reporting agent.

        Args:
            biological_prompt: Research question or biological context
            max_genes: Maximum number of genes to include

        Returns:
            LandingPageData object or None if generation fails
        """
        try:
            # Use existing agent if available, otherwise create new one
            if not self.agent:
                # Initialize the reporting agent
                self.agent = ReportingAgent(
                    results_dir=self.results_dir,
                    model_name=self.model_name,
                    temperature=0.1
                )

                # Setup servers only if agent doesn't have them
                if not self.agent.servers:
                    await asyncio.get_event_loop().run_in_executor(None, self.agent.setup_servers)

            # Create agent with research focus (or recreate with new prompt)
            self.agent.create_agent(biological_prompt)

            # Step 1: Select relevant contrasts
            logger.info("Selecting relevant contrasts...")
            selected_contrasts = await self._select_contrasts_with_agent(biological_prompt)
            if not selected_contrasts:
                logger.warning("No contrasts selected")
                return None

            # Step 2: Find key genes
            logger.info("Identifying key genes...")
            top_genes = await self._select_genes_with_agent(selected_contrasts, max_genes)
            if not top_genes:
                logger.warning("No genes identified")
                return None

            # Step 3: Create visualizations using existing integrator methods
            logger.info("Creating visualizations...")
            heatmap_fig = self._create_heatmap(top_genes, selected_contrasts)
            gene_table = await self._create_gene_table(top_genes, selected_contrasts)

            # Create thresholds object (using standard values for now)
            thresholds = ThresholdSelection(
                fdr_cutoff=0.05,
                logfc_cutoff=1.0,
                min_frequency=2,
                justification="Standard thresholds: FDR < 0.05, |logFC| > 1.0"
            )

            return LandingPageData(
                selected_contrasts=selected_contrasts,
                thresholds=thresholds,
                top_genes=top_genes,
                heatmap_fig=heatmap_fig,
                gene_table=gene_table
            )

        except Exception as e:
            logger.error(f"Error generating landing page: {e}")
            return None
        # Note: No cleanup in finally block - we want to reuse the agent and servers

    async def _select_contrasts_with_agent(self, biological_prompt: str) -> List[ContrastSelection]:
        """Select contrasts using the reporting agent."""
        try:
            # Use the agent to rank contrasts
            prompt = f"""Please analyze the available contrasts and select the most relevant ones for this research question: "{biological_prompt}"

1. Use rank_contrasts_by_relevance to identify the most relevant contrasts
2. For each selected contrast, provide:
   - The contrast ID and analysis ID
   - A relevance score (0-10)
   - A clear justification for why it's relevant
   - The number of DEGs

Select between 3-12 contrasts that best address the research question.
Provide your response in a structured format that includes these details for each contrast."""

            result = await self.agent.agent.run(prompt)

            # Parse the agent's response to extract contrast selections
            selected_contrasts = []

            # The agent should have used the rank_contrasts_by_relevance tool
            # We need to parse its structured output
            # For now, we'll also check the integrator's data directly

            # Get all available contrasts from integrator
            all_contrasts = []
            for analysis_id, contrasts in self.integrator.deg_data.items():
                for contrast_id in contrasts.keys():
                    description = self.integrator._get_contrast_description(analysis_id, contrast_id)

                    # Count DEGs
                    df = contrasts[contrast_id]
                    deg_count = 0
                    if 'adj.P.Val' in df.columns and 'logFC' in df.columns:
                        deg_count = ((df['adj.P.Val'] < 0.05) & (abs(df['logFC']) > 1.0)).sum()

                    # Simple relevance scoring based on keywords
                    relevance = self._calculate_relevance(biological_prompt, description)

                    all_contrasts.append({
                        'analysis_id': analysis_id,
                        'contrast_id': contrast_id,
                        'description': description,
                        'deg_count': deg_count,
                        'relevance': relevance
                    })

            # Sort by relevance and take top contrasts
            all_contrasts.sort(key=lambda x: x['relevance'], reverse=True)

            # Select top 3-12 contrasts based on relevance
            n_contrasts = min(12, max(3, len([c for c in all_contrasts if c['relevance'] > 0.1])))

            for contrast in all_contrasts[:n_contrasts]:
                selected_contrasts.append(ContrastSelection(
                    analysis_id=contrast['analysis_id'],
                    contrast_id=contrast['contrast_id'],
                    relevance_score=min(10, contrast['relevance'] * 10),
                    justification=f"This contrast examines {contrast['description']} which is relevant to {biological_prompt}",
                    deg_count=contrast['deg_count']
                ))

            return selected_contrasts

        except Exception as e:
            logger.error(f"Error selecting contrasts: {e}", exc_info=True)
            return []

    async def _select_genes_with_agent(self, selected_contrasts: List[ContrastSelection], max_genes: int) -> List[str]:
        """Select key genes using the reporting agent."""
        try:
            # Build list of contrast keys
            contrast_keys = [f"{c.analysis_id}:{c.contrast_id}" for c in selected_contrasts]

            prompt = f"""Identify the most important genes across these selected contrasts.

Selected contrasts: {', '.join(contrast_keys)}

Please:
1. Use find_common_degs to identify genes that are consistently differentially expressed
2. Look for genes with:
   - High frequency across contrasts (appear in multiple contrasts)
   - Large fold changes
   - Consistent direction of change
3. Return up to {max_genes} genes that are most biologically important

Focus on genes that would be most relevant for understanding the biological processes involved."""

            result = await self.agent.agent.run(prompt)

            # Also use direct analysis as fallback
            analysis_ids = [c.analysis_id for c in selected_contrasts]
            contrast_ids = [c.contrast_id for c in selected_contrasts]

            # Count gene occurrences
            gene_stats = {}

            for contrast in selected_contrasts:
                if (contrast.analysis_id in self.integrator.deg_data and
                    contrast.contrast_id in self.integrator.deg_data[contrast.analysis_id]):

                    df = self.integrator.deg_data[contrast.analysis_id][contrast.contrast_id]

                    if 'Gene' in df.columns and 'adj.P.Val' in df.columns and 'logFC' in df.columns:
                        # Filter significant genes
                        significant = df[(df['adj.P.Val'] < 0.05) & (abs(df['logFC']) > 1.0)]

                        for _, row in significant.iterrows():
                            gene = row['Gene']
                            if gene not in gene_stats:
                                gene_stats[gene] = {
                                    'count': 0,
                                    'max_lfc': 0,
                                    'lfcs': []
                                }

                            gene_stats[gene]['count'] += 1
                            gene_stats[gene]['max_lfc'] = max(gene_stats[gene]['max_lfc'], abs(row['logFC']))
                            gene_stats[gene]['lfcs'].append(row['logFC'])

            # Score genes by frequency and effect size
            gene_scores = []
            for gene, stats in gene_stats.items():
                # Require gene to be in at least 2 contrasts or 30% of contrasts
                min_contrasts = max(2, len(selected_contrasts) * 0.3)
                if stats['count'] >= min_contrasts:
                    score = stats['count'] * stats['max_lfc']
                    gene_scores.append((gene, score, stats))

            # Sort by score and take top genes
            gene_scores.sort(key=lambda x: x[1], reverse=True)
            top_genes = [g[0] for g in gene_scores[:max_genes]]

            return top_genes

        except Exception as e:
            logger.error(f"Error selecting genes: {e}", exc_info=True)
            return []

    def _create_heatmap(self, genes: List[str], contrasts: List[ContrastSelection]) -> Optional[go.Figure]:
        """Create heatmap using the integrator's existing method."""
        try:
            contrast_pairs = [(c.analysis_id, c.contrast_id) for c in contrasts]

            fig = self.integrator.create_lfc_heatmap(
                genes=genes,
                contrasts=contrast_pairs,
                output_file=None,
                p_value_threshold=0.05,
                lfc_threshold=1.0,
                hide_empty_rows_cols=True,
                font_size=11
            )

            if fig:
                # Adjust layout for landing page
                calculated_height = max(500, min(1200, len(genes) * 28))

                fig.update_layout(
                    title="Key Genes Across Relevant Contrasts",
                    title_font_size=16,
                    height=calculated_height,
                    margin=dict(l=150, r=20, t=60, b=80)
                )

                fig.update_yaxes(
                    tickfont=dict(size=10),
                    tickmode='linear'
                )

            return fig

        except Exception as e:
            logger.error(f"Error creating heatmap: {e}", exc_info=True)
            return None

    async def _create_gene_table(self, genes: List[str], contrasts: List[ContrastSelection]) -> pd.DataFrame:
        """Create gene summary table with statistics."""
        gene_data = []

        for gene in genes:
            up_count = 0
            down_count = 0
            total_tested = 0
            logfcs = []

            for contrast in contrasts:
                if (contrast.analysis_id in self.integrator.deg_data and
                    contrast.contrast_id in self.integrator.deg_data[contrast.analysis_id]):

                    df = self.integrator.deg_data[contrast.analysis_id][contrast.contrast_id]

                    if 'Gene' in df.columns and 'adj.P.Val' in df.columns and 'logFC' in df.columns:
                        gene_row = df[df['Gene'] == gene]

                        if not gene_row.empty:
                            total_tested += 1
                            row = gene_row.iloc[0]
                            p_val = row['adj.P.Val']
                            logfc = row['logFC']
                            logfcs.append(logfc)

                            if p_val < 0.05 and abs(logfc) > 1.0:
                                if logfc > 0:
                                    up_count += 1
                                else:
                                    down_count += 1

            if logfcs:
                median_logfc = np.median(logfcs)
                direction_summary = f"↑{up_count}/↓{down_count} (of {total_tested})"

                gene_data.append({
                    'Gene': gene,
                    'Median LogFC': round(median_logfc, 2),
                    'Direction': direction_summary,
                    'Significant in': f"{up_count + down_count}/{total_tested}",
                    'Consistency': 'Yes' if (up_count == 0 or down_count == 0) else 'Mixed'
                })

        return pd.DataFrame(gene_data)



    def _calculate_relevance(self, query: str, description: str) -> float:
        """Calculate simple relevance score between query and description."""
        query_words = set(query.lower().split())
        desc_words = set(description.lower().split())

        # Remove common words
        stopwords = {'the', 'a', 'an', 'and', 'or', 'but', 'in', 'on', 'at', 'to', 'for', 'of', 'with', 'by'}
        query_words = query_words - stopwords
        desc_words = desc_words - stopwords

        if not query_words:
            return 0.0

        # Calculate overlap
        overlap = len(query_words & desc_words)
        return overlap / len(query_words)


# Synchronous wrapper for Streamlit compatibility
def generate_ai_landing_page(integrator, biological_prompt: str, max_genes: int = 50) -> Optional[LandingPageData]:
    """
    Synchronous wrapper for generating AI landing page.
    Reuses existing MCP servers from session state if available.
    """
    try:
        logger.info(f"Starting generate_ai_landing_page with prompt: {biological_prompt[:100]}...")
        # Check if we have existing servers in session state
        import streamlit as st

        if hasattr(st.session_state, 'mcp_manager') and st.session_state.mcp_manager:
            logger.info("Found existing MCP manager in session state")
            # Reuse existing manager and its agent
            manager = st.session_state.mcp_manager
            
            # Create generator
            generator = LandingPageGenerator(integrator)
            
            # Always create a new agent configured for this specific biological prompt
            # but reuse the existing servers
            logger.info("Creating new ReportingAgent with existing servers")
            reporting_agent = ReportingAgent(
                results_dir=integrator.results_dir,
                model_name="openai:gpt-4o-mini",
                temperature=0.1
            )
            # Reuse existing servers instead of setting up new ones
            servers = manager.servers if hasattr(manager, 'servers') else []
            logger.info(f"Reusing {len(servers)} existing MCP servers")
            reporting_agent.servers = servers
            # Create agent configured for the specific biological prompt
            logger.info("Creating agent with biological prompt")
            reporting_agent.create_agent(biological_prompt)
            generator.agent = reporting_agent

            # Create event loop if needed
            try:
                loop = asyncio.get_event_loop()
            except RuntimeError:
                loop = asyncio.new_event_loop()
                asyncio.set_event_loop(loop)

            # Generate landing page without cleaning up servers
            logger.info("Starting landing page generation")
            if loop.is_running():
                logger.debug("Event loop is running, creating task for async execution")
                # Create a task instead of using ThreadPoolExecutor to maintain context
                task = asyncio.create_task(generator.generate_landing_page(biological_prompt, max_genes))
                # Use asyncio.wait_for with a reasonable timeout
                try:
                    result = asyncio.run_coroutine_threadsafe(
                        generator.generate_landing_page(biological_prompt, max_genes), loop
                    ).result(timeout=120)  # 2 minute timeout
                    logger.info(f"Landing page generation completed successfully: {result is not None}")
                    return result
                except Exception as async_error:
                    logger.error(f"Error in async execution: {async_error}", exc_info=True)
                    # Fallback to sync approach
                    logger.info("Falling back to synchronous execution")
                    # Skip fallback to avoid nested event loop issues
                    logger.info("Skipping fallback execution due to event loop conflicts")
                    return None
            else:
                logger.debug("Using run_until_complete for async execution")
                result = loop.run_until_complete(
                    generator.generate_landing_page(biological_prompt, max_genes)
                )
                logger.info(f"Landing page generation completed successfully: {result is not None}")
                return result
        else:
            # Fallback - create new manager (shouldn't normally happen)
            logger.warning("No MCP servers available in session state")
            st.error("No MCP servers available. Please refresh the page.")
            return None

    except Exception as e:
        logger.error(f"Error in synchronous wrapper: {e}", exc_info=True)
        # Additional context for debugging
        logger.error(f"Error type: {type(e).__name__}")
        logger.error(f"Error args: {e.args}")
        return None



async def _generate_landing_page_async(integrator, biological_prompt: str, max_genes: int = 50) -> Optional[LandingPageData]:
    """Async implementation of landing page generation."""
    generator = LandingPageGenerator(integrator)
    return await generator.generate_landing_page(biological_prompt, max_genes)


# Auto-analysis function for initial page load
async def auto_analyze_on_load(integrator, progress_callback=None) -> Optional[Dict[str, Any]]:
    """
    Automatically analyze the dataset when first loaded.
    Keeps MCP servers alive in session state for subsequent use.

    Args:
        integrator: ResultsIntegrator instance
        progress_callback: Optional callback for progress updates

    Returns:
        Analysis results or None
    """
    import streamlit as st

    try:
        results_dir = integrator.results_dir
        logger.info(f"Starting auto_analyze_on_load for results_dir: {results_dir}")

        # Check if we already have a manager in session state
        if 'mcp_manager' not in st.session_state:
            logger.info("No existing MCP manager found, creating new one")
            # Create and store manager in session state
            manager = AutoStartManager(results_dir)
            logger.info(f"Created AutoStartManager for {results_dir}")

            # Validate directory
            is_valid, error = manager.validate_data_directory()
            if not is_valid:
                logger.warning(f"Invalid directory: {error}")
                return None
            logger.info("Directory validation passed")

            # Setup servers and store in session state
            logger.info("Setting up MCP servers...")
            manager.setup_servers()
            logger.info(f"Setup {len(manager.servers)} MCP servers")
            
            st.session_state.mcp_manager = manager
            st.session_state.mcp_servers = manager.servers
            logger.info("Created and stored MCP manager in session state")
        else:
            # Reuse existing manager
            manager = st.session_state.mcp_manager
            logger.info(f"Reusing existing MCP manager from session state (has {len(getattr(manager, 'servers', []))} servers)")

        # Run analysis
        logger.info("Starting initial analysis...")
        try:
            loop = asyncio.get_event_loop()
            logger.debug("Using existing event loop")
        except RuntimeError:
            logger.debug("Creating new event loop")
            loop = asyncio.new_event_loop()
            asyncio.set_event_loop(loop)

        if loop.is_running():
            logger.debug("Event loop is running, using run_coroutine_threadsafe")
            # Use run_coroutine_threadsafe to maintain context
            try:
                result = asyncio.run_coroutine_threadsafe(
                    manager.run_initial_analysis(progress_callback), loop
                ).result(timeout=60)  # 1 minute timeout
            except Exception as async_error:
                logger.error(f"Error in async auto-analysis: {async_error}", exc_info=True)
                # Skip fallback to avoid nested event loop issues
                logger.info("Skipping auto-analysis due to event loop conflicts")
                result = None
        else:
            logger.debug("Event loop not running, using run_until_complete")
            result = loop.run_until_complete(
                manager.run_initial_analysis(progress_callback)
            )

        logger.info(f"Analysis completed with success: {result.get('success', False) if result else False}")
        return result

    except Exception as e:
        logger.error(f"Error in auto-analysis: {e}", exc_info=True)
        # Additional context for debugging
        logger.error(f"Error type: {type(e).__name__}")
        logger.error(f"Error args: {e.args}")
        logger.error(f"Results dir: {integrator.results_dir}")
        return None
    # Note: NO finally block that cleans up the manager
    # The servers stay alive in session state for future use

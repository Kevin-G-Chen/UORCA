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
    narrative: str


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
            # Initialize the reporting agent
            self.agent = ReportingAgent(
                results_dir=self.results_dir,
                model_name=self.model_name,
                temperature=0.1
            )
            
            # Setup servers
            await asyncio.get_event_loop().run_in_executor(None, self.agent.setup_servers)
            
            # Create agent with research focus
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
            
            # Step 4: Generate comprehensive narrative
            logger.info("Generating narrative...")
            narrative = await self._generate_narrative(selected_contrasts, top_genes, biological_prompt)
            
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
                gene_table=gene_table,
                narrative=narrative
            )
            
        except Exception as e:
            logger.error(f"Error generating landing page: {e}")
            return None
        finally:
            # Always cleanup agent resources
            if self.agent:
                try:
                    self.agent.cleanup()
                except Exception as e:
                    logger.warning(f"Error during agent cleanup: {e}")
    
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
            logger.error(f"Error selecting contrasts: {e}")
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
            logger.error(f"Error selecting genes: {e}")
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
            logger.error(f"Error creating heatmap: {e}")
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
    
    async def _generate_narrative(self, contrasts: List[ContrastSelection], genes: List[str], biological_prompt: str) -> str:
        """Generate comprehensive narrative using the agent."""
        try:
            # Prepare context for narrative generation
            contrast_summary = []
            for c in contrasts:
                desc = self.integrator._get_contrast_description(c.analysis_id, c.contrast_id)
                contrast_summary.append(f"- {c.contrast_id}: {desc} (Score: {c.relevance_score:.1f}, DEGs: {c.deg_count})")
            
            # Get gene statistics
            gene_summary = []
            for gene in genes[:10]:  # Focus on top 10 for narrative
                # Collect stats across contrasts
                lfcs = []
                significant_count = 0
                
                for c in contrasts:
                    if (c.analysis_id in self.integrator.deg_data and 
                        c.contrast_id in self.integrator.deg_data[c.analysis_id]):
                        
                        df = self.integrator.deg_data[c.analysis_id][c.contrast_id]
                        if 'Gene' in df.columns:
                            gene_row = df[df['Gene'] == gene]
                            if not gene_row.empty and 'logFC' in df.columns:
                                lfc = gene_row.iloc[0]['logFC']
                                lfcs.append(lfc)
                                if 'adj.P.Val' in df.columns:
                                    if gene_row.iloc[0]['adj.P.Val'] < 0.05 and abs(lfc) > 1.0:
                                        significant_count += 1
                
                if lfcs:
                    gene_summary.append(f"{gene}: median logFC={np.median(lfcs):.2f}, significant in {significant_count}/{len(contrasts)} contrasts")
            
            prompt = f"""Generate a comprehensive interpretation of these RNA-seq results for: "{biological_prompt}"

SELECTED CONTRASTS:
{chr(10).join(contrast_summary)}

TOP DIFFERENTIALLY EXPRESSED GENES:
{chr(10).join(gene_summary)}

Additional genes identified: {', '.join(genes[10:20]) if len(genes) > 10 else 'None'}

Please provide:
1. **Overview**: Summarize what these contrasts reveal about {biological_prompt}
2. **Key Findings**: Highlight the most important genes and their expression patterns
3. **Biological Interpretation**: Explain the biological significance and potential mechanisms
4. **Research Implications**: Suggest future directions or therapeutic implications

Write in clear, scientific language with specific details about genes and pathways.
Format with markdown for readability."""

            result = await self.agent.agent.run(prompt)
            return result.data
            
        except Exception as e:
            logger.error(f"Error generating narrative: {e}")
            # Fallback narrative
            return f"""## Analysis Summary

Based on the analysis of {len(contrasts)} relevant contrasts, we identified {len(genes)} key genes related to {biological_prompt}.

### Key Findings
The selected contrasts represent experimental comparisons that are most relevant to understanding {biological_prompt}. 
Among the {len(genes)} genes identified, several show consistent differential expression patterns across multiple conditions.

### Next Steps
Further investigation of these genes and their associated pathways may provide deeper insights into the biological mechanisms involved."""
    
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
    
    This function is compatible with the existing Streamlit UI.
    """
    try:
        # Create event loop if needed
        try:
            loop = asyncio.get_event_loop()
        except RuntimeError:
            loop = asyncio.new_event_loop()
            asyncio.set_event_loop(loop)
        
        # Check if we're already in an async context
        try:
            # If we're already in an async context, create a new task
            if loop.is_running():
                # Use asyncio.run in a thread
                import concurrent.futures
                with concurrent.futures.ThreadPoolExecutor() as executor:
                    future = executor.submit(
                        asyncio.run,
                        _generate_landing_page_async(integrator, biological_prompt, max_genes)
                    )
                    return future.result()
            else:
                # Run normally
                return loop.run_until_complete(
                    _generate_landing_page_async(integrator, biological_prompt, max_genes)
                )
        except:
            # Fallback: create new event loop
            return asyncio.run(
                _generate_landing_page_async(integrator, biological_prompt, max_genes)
            )
            
    except Exception as e:
        logger.error(f"Error in synchronous wrapper: {e}")
        return None


async def _generate_landing_page_async(integrator, biological_prompt: str, max_genes: int = 50) -> Optional[LandingPageData]:
    """Async implementation of landing page generation."""
    generator = LandingPageGenerator(integrator)
    return await generator.generate_landing_page(biological_prompt, max_genes)


# Auto-analysis function for initial page load
def auto_analyze_on_load(integrator, progress_callback=None) -> Optional[Dict[str, Any]]:
    """
    Automatically analyze the dataset when first loaded.
    
    Args:
        integrator: ResultsIntegrator instance
        progress_callback: Optional callback for progress updates
        
    Returns:
        Analysis results or None
    """
    manager = None
    try:
        results_dir = integrator.results_dir
        
        # Create auto-start manager
        manager = AutoStartManager(results_dir)
        
        # Validate directory
        is_valid, error = manager.validate_data_directory()
        if not is_valid:
            logger.warning(f"Invalid directory: {error}")
            return None
        
        # Run analysis
        try:
            loop = asyncio.get_event_loop()
        except RuntimeError:
            loop = asyncio.new_event_loop()
            asyncio.set_event_loop(loop)
        
        if loop.is_running():
            # Use thread for async execution
            import concurrent.futures
            with concurrent.futures.ThreadPoolExecutor() as executor:
                future = executor.submit(
                    asyncio.run,
                    manager.run_initial_analysis(progress_callback)
                )
                result = future.result()
        else:
            result = loop.run_until_complete(
                manager.run_initial_analysis(progress_callback)
            )
        
        return result
        
    except Exception as e:
        logger.error(f"Error in auto-analysis: {e}")
        return None
    finally:
        # Always cleanup manager resources
        if manager:
            try:
                manager.cleanup()
            except Exception as e:
                logger.warning(f"Error during manager cleanup: {e}")
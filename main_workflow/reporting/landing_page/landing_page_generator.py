#!/usr/bin/env python3
"""
UORCA Landing Page Generator - AI-Assisted Interpretation with Minimal Interactivity

Creates a self-contained landing page presenting the most relevant, high-confidence 
differential-expression results automatically selected for a given biological prompt.
"""

import os
import json
import logging
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Any, Optional, Tuple, Set
from dataclasses import dataclass
from openai import OpenAI
from dotenv import load_dotenv

# Import existing components
from ..ResultsIntegration import ResultsIntegrator

# Set up logging
logger = logging.getLogger(__name__)

# Load environment
load_dotenv()

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
    heatmap_fig: go.Figure
    gene_table: pd.DataFrame
    narrative: str

class LandingPageGenerator:
    """
    Generates AI-assisted landing pages for UORCA results with minimal user interaction.
    
    Automatically selects relevant contrasts, appropriate thresholds, and generates
    interpretive narratives for RNA-seq differential expression results.
    """
    
    def __init__(self, results_dir: str, biological_prompt: str = None):
        """
        Initialize the landing page generator.
        
        Parameters:
        -----------
        results_dir : str
            Directory containing UORCA analysis results
        biological_prompt : str, optional
            Biological question or context to guide analysis selection
        """
        self.results_dir = results_dir
        self.biological_prompt = biological_prompt or "General differential expression analysis"
        
        # Initialize OpenAI client
        self.client = OpenAI()
        
        # Initialize the results integrator
        self.integrator = ResultsIntegrator(results_dir=results_dir)
        self.integrator.load_data()
        
        # Cache for expensive operations
        self._contrast_embeddings = None
        
        logger.info(f"Initialized LandingPageGenerator for {results_dir}")
    
    def generate_landing_page(self, 
                            max_contrasts: int = 8,
                            max_genes: int = 50,
                            output_path: str = None) -> LandingPageData:
        """
        Generate a complete landing page with automatic selections.
        
        Parameters:
        -----------
        max_contrasts : int
            Maximum number of contrasts to include
        max_genes : int
            Maximum number of genes to include in visualizations
        output_path : str, optional
            Path to save the landing page HTML
            
        Returns:
        --------
        LandingPageData containing all components
        """
        logger.info("Starting automatic landing page generation")
        
        # Step 1: Automatic contrast selection
        selected_contrasts = self._select_relevant_contrasts(max_contrasts)
        logger.info(f"Selected {len(selected_contrasts)} relevant contrasts")
        
        # Step 2: Automatic threshold selection
        thresholds = self._select_optimal_thresholds(selected_contrasts)
        logger.info(f"Selected thresholds: FDR={thresholds.fdr_cutoff}, logFC={thresholds.logfc_cutoff}")
        
        # Step 3: Aggregate and rank DEGs
        top_genes = self._aggregate_and_rank_genes(selected_contrasts, thresholds, max_genes)
        logger.info(f"Identified {len(top_genes)} top differentially expressed genes")
        
        # Step 4: Create visualizations
        heatmap_fig = self._create_landing_heatmap(top_genes, selected_contrasts, thresholds)
        gene_table = self._create_gene_table(top_genes, selected_contrasts, thresholds)
        
        # Step 5: Generate narrative interpretation
        narrative = self._generate_narrative(selected_contrasts, thresholds, top_genes)
        
        # Create landing page data
        landing_data = LandingPageData(
            selected_contrasts=selected_contrasts,
            thresholds=thresholds,
            top_genes=top_genes,
            heatmap_fig=heatmap_fig,
            gene_table=gene_table,
            narrative=narrative
        )
        
        # Generate HTML if output path specified
        if output_path:
            self._generate_html(landing_data, output_path)
            logger.info(f"Landing page saved to {output_path}")
        
        return landing_data
    
    def _select_relevant_contrasts(self, max_contrasts: int) -> List[ContrastSelection]:
        """
        Automatically select the most relevant contrasts using AI scoring.
        """
        logger.info("Selecting relevant contrasts using AI scoring")
        
        # Get all available contrasts
        all_contrasts = []
        for analysis_id, contrasts in self.integrator.deg_data.items():
            for contrast_id in contrasts.keys():
                # Get basic info about each contrast
                description = self.integrator._get_contrast_description(analysis_id, contrast_id)
                
                # Count potential DEGs (using lenient thresholds for initial filtering)
                df = contrasts[contrast_id]
                deg_count = 0
                if 'adj.P.Val' in df.columns and 'logFC' in df.columns:
                    deg_count = ((df['adj.P.Val'] < 0.1) & (abs(df['logFC']) > 0.5)).sum()
                
                all_contrasts.append({
                    'analysis_id': analysis_id,
                    'contrast_id': contrast_id,
                    'description': description,
                    'deg_count': deg_count
                })
        
        if not all_contrasts:
            logger.warning("No contrasts found")
            return []
        
        # Score contrasts using LLM
        scored_contrasts = self._score_contrasts_with_llm(all_contrasts)
        
        # Select top contrasts
        top_contrasts = sorted(scored_contrasts, key=lambda x: x.relevance_score, reverse=True)[:max_contrasts]
        
        return top_contrasts
    
    def _score_contrasts_with_llm(self, contrasts: List[Dict]) -> List[ContrastSelection]:
        """Score contrasts using LLM based on biological relevance."""
        
        # Prepare prompt for contrast scoring
        contrast_descriptions = []
        for i, contrast in enumerate(contrasts):
            contrast_descriptions.append(
                f"{i+1}. {contrast['analysis_id']}_{contrast['contrast_id']}: "
                f"{contrast['description']} (Potential DEGs: {contrast['deg_count']})"
            )
        
        prompt = f"""
You are an expert RNA-seq analyst. Score the biological relevance of these differential expression contrasts for the research context: "{self.biological_prompt}"

Available contrasts:
{chr(10).join(contrast_descriptions)}

For each contrast, provide:
1. Relevance score (0-10, where 10 is most relevant)
2. Brief justification (1-2 sentences)

Consider:
- Biological significance for the research question
- Number of potential DEGs (more is generally better)
- Clarity of experimental design
- Scientific interest and interpretability

Respond in JSON format:
{{
  "scored_contrasts": [
    {{
      "contrast_number": 1,
      "relevance_score": 8.5,
      "justification": "This contrast directly addresses the research question by comparing..."
    }},
    ...
  ]
}}
"""
        
        try:
            response = self.client.chat.completions.create(
                model="gpt-4o-mini",
                messages=[{"role": "user", "content": prompt}],
                response_format={"type": "json_object"}
            )
            
            result = json.loads(response.choices[0].message.content)
            scored_contrasts = []
            
            for score_data in result.get("scored_contrasts", []):
                idx = score_data["contrast_number"] - 1
                if 0 <= idx < len(contrasts):
                    contrast = contrasts[idx]
                    scored_contrasts.append(ContrastSelection(
                        analysis_id=contrast['analysis_id'],
                        contrast_id=contrast['contrast_id'],
                        relevance_score=score_data["relevance_score"],
                        justification=score_data["justification"],
                        deg_count=contrast['deg_count']
                    ))
            
            return scored_contrasts
            
        except Exception as e:
            logger.warning(f"LLM scoring failed: {e}, using fallback scoring")
            return self._fallback_contrast_scoring(contrasts)
    
    def _fallback_contrast_scoring(self, contrasts: List[Dict]) -> List[ContrastSelection]:
        """Fallback scoring based on DEG counts and description length."""
        scored_contrasts = []
        
        for contrast in contrasts:
            # Simple heuristic: score based on DEG count and description informativeness
            base_score = min(10, contrast['deg_count'] / 100 * 8)  # Up to 8 points for DEG count
            desc_score = min(2, len(contrast['description'].split()) / 20 * 2)  # Up to 2 points for description
            
            total_score = base_score + desc_score
            justification = f"Selected based on {contrast['deg_count']} potential DEGs and experimental design clarity."
            
            scored_contrasts.append(ContrastSelection(
                analysis_id=contrast['analysis_id'],
                contrast_id=contrast['contrast_id'],
                relevance_score=total_score,
                justification=justification,
                deg_count=contrast['deg_count']
            ))
        
        return scored_contrasts
    
    def _select_optimal_thresholds(self, selected_contrasts: List[ContrastSelection]) -> ThresholdSelection:
        """Automatically select optimal statistical thresholds using LLM."""
        
        # Calculate dataset statistics
        total_contrasts = len(selected_contrasts)
        total_deg_counts = [c.deg_count for c in selected_contrasts]
        median_degs = np.median(total_deg_counts) if total_deg_counts else 0
        
        # Get sample size information
        sample_info = []
        for contrast in selected_contrasts:
            if contrast.analysis_id in self.integrator.analysis_info:
                n_samples = self.integrator.analysis_info[contrast.analysis_id].get('number_of_samples', 0)
                sample_info.append(n_samples)
        
        median_samples = np.median(sample_info) if sample_info else 0
        
        prompt = f"""
You are an expert statistician selecting optimal thresholds for RNA-seq differential expression analysis.

Dataset characteristics:
- Number of contrasts: {total_contrasts}
- Median potential DEGs per contrast: {median_degs:.0f}
- Median sample size: {median_samples:.0f}
- Research context: {self.biological_prompt}

Select appropriate thresholds considering:
1. FDR cutoff: Stricter if few contrasts or small samples (0.01-0.1 range)
2. Log fold change cutoff: Higher for well-powered studies (0.5-2.0 range)
3. Minimum frequency: How many contrasts a gene must appear in (1 to {max(1, total_contrasts//3)})

Provide scientific justification for each choice.

Respond in JSON format:
{{
  "fdr_cutoff": 0.05,
  "logfc_cutoff": 1.0,
  "min_frequency": 2,
  "justification": "Selected FDR=0.05 for balanced sensitivity/specificity given the dataset size..."
}}
"""
        
        try:
            response = self.client.chat.completions.create(
                model="gpt-4o-mini",
                messages=[{"role": "user", "content": prompt}],
                response_format={"type": "json_object"}
            )
            
            result = json.loads(response.choices[0].message.content)
            return ThresholdSelection(
                fdr_cutoff=result["fdr_cutoff"],
                logfc_cutoff=result["logfc_cutoff"],
                min_frequency=result["min_frequency"],
                justification=result["justification"]
            )
            
        except Exception as e:
            logger.warning(f"LLM threshold selection failed: {e}, using fallback")
            return self._fallback_threshold_selection(total_contrasts, median_samples)
    
    def _fallback_threshold_selection(self, n_contrasts: int, median_samples: float) -> ThresholdSelection:
        """Fallback threshold selection using heuristics."""
        
        # Heuristic rules
        if n_contrasts < 3 or median_samples < 6:
            fdr_cutoff = 0.01  # Stricter for small studies
            logfc_cutoff = 1.5
        elif median_samples > 20:
            fdr_cutoff = 0.05  # Standard for well-powered
            logfc_cutoff = 1.0
        else:
            fdr_cutoff = 0.05
            logfc_cutoff = 1.2
        
        min_frequency = max(1, min(3, n_contrasts // 3))
        
        justification = f"Heuristic selection: FDR={fdr_cutoff} and logFC={logfc_cutoff} based on {n_contrasts} contrasts and median {median_samples:.0f} samples."
        
        return ThresholdSelection(
            fdr_cutoff=fdr_cutoff,
            logfc_cutoff=logfc_cutoff,
            min_frequency=min_frequency,
            justification=justification
        )
    
    def _aggregate_and_rank_genes(self, 
                                selected_contrasts: List[ContrastSelection],
                                thresholds: ThresholdSelection,
                                max_genes: int) -> List[str]:
        """Aggregate DEGs across contrasts and rank by frequency and effect size."""
        
        gene_stats = {}  # gene -> {frequency, max_logfc, contrasts}
        
        for contrast in selected_contrasts:
            if contrast.analysis_id in self.integrator.deg_data:
                if contrast.contrast_id in self.integrator.deg_data[contrast.analysis_id]:
                    df = self.integrator.deg_data[contrast.analysis_id][contrast.contrast_id]
                    
                    if 'Gene' in df.columns and 'adj.P.Val' in df.columns and 'logFC' in df.columns:
                        # Filter significant genes
                        sig_genes = df[
                            (df['adj.P.Val'] < thresholds.fdr_cutoff) & 
                            (abs(df['logFC']) > thresholds.logfc_cutoff)
                        ]
                        
                        for _, row in sig_genes.iterrows():
                            gene = row['Gene']
                            logfc = abs(row['logFC'])
                            
                            if gene not in gene_stats:
                                gene_stats[gene] = {
                                    'frequency': 0,
                                    'max_logfc': 0,
                                    'contrasts': []
                                }
                            
                            gene_stats[gene]['frequency'] += 1
                            gene_stats[gene]['max_logfc'] = max(gene_stats[gene]['max_logfc'], logfc)
                            gene_stats[gene]['contrasts'].append(f"{contrast.analysis_id}_{contrast.contrast_id}")
        
        # Filter by minimum frequency
        filtered_genes = {
            gene: stats for gene, stats in gene_stats.items() 
            if stats['frequency'] >= thresholds.min_frequency
        }
        
        # Rank by frequency then by max log fold change
        ranked_genes = sorted(
            filtered_genes.items(),
            key=lambda x: (x[1]['frequency'], x[1]['max_logfc']),
            reverse=True
        )
        
        top_genes = [gene for gene, _ in ranked_genes[:max_genes]]
        
        logger.info(f"Ranked {len(filtered_genes)} genes meeting frequency threshold, selected top {len(top_genes)}")
        
        return top_genes
    
    def _create_landing_heatmap(self, 
                              top_genes: List[str],
                              selected_contrasts: List[ContrastSelection],
                              thresholds: ThresholdSelection) -> go.Figure:
        """Create an interactive heatmap for the landing page."""
        
        # Use the integrator's heatmap creation with selected parameters
        contrast_pairs = [(c.analysis_id, c.contrast_id) for c in selected_contrasts]
        
        fig = self.integrator.create_lfc_heatmap(
            genes=top_genes,
            contrasts=contrast_pairs,
            output_file=None,
            p_value_threshold=thresholds.fdr_cutoff,
            lfc_threshold=thresholds.logfc_cutoff,
            hide_empty_rows_cols=True,
            font_size=11
        )
        
        if fig:
            # Update title for landing page
            fig.update_layout(
                title="Key Differentially Expressed Genes Across Selected Contrasts",
                title_font_size=16,
                height=max(400, min(800, len(top_genes) * 20))
            )
        
        return fig
    
    def _create_gene_table(self, 
                         top_genes: List[str],
                         selected_contrasts: List[ContrastSelection],
                         thresholds: ThresholdSelection) -> pd.DataFrame:
        """Create a summary table of top genes."""
        
        gene_data = []
        
        for gene in top_genes:
            # Collect statistics across all selected contrasts
            appearances = 0
            max_logfc = 0
            min_pval = 1.0
            contrast_list = []
            
            for contrast in selected_contrasts:
                if contrast.analysis_id in self.integrator.deg_data:
                    if contrast.contrast_id in self.integrator.deg_data[contrast.analysis_id]:
                        df = self.integrator.deg_data[contrast.analysis_id][contrast.contrast_id]
                        
                        if 'Gene' in df.columns:
                            gene_row = df[df['Gene'] == gene]
                            if not gene_row.empty:
                                row = gene_row.iloc[0]
                                
                                # Check if meets significance criteria
                                if ('adj.P.Val' in df.columns and 'logFC' in df.columns and
                                    row['adj.P.Val'] < thresholds.fdr_cutoff and
                                    abs(row['logFC']) > thresholds.logfc_cutoff):
                                    
                                    appearances += 1
                                    max_logfc = max(max_logfc, abs(row['logFC']))
                                    min_pval = min(min_pval, row['adj.P.Val'])
                                    contrast_list.append(f"{contrast.analysis_id}_{contrast.contrast_id}")
            
            if appearances > 0:
                gene_data.append({
                    'Gene': gene,
                    'Frequency': appearances,
                    'Max_LogFC': round(max_logfc, 2),
                    'Min_AdjPVal': f"{min_pval:.2e}" if min_pval < 0.01 else f"{min_pval:.3f}",
                    'Contrasts': '; '.join(contrast_list[:3]) + ('...' if len(contrast_list) > 3 else '')
                })
        
        df = pd.DataFrame(gene_data)
        df = df.sort_values(['Frequency', 'Max_LogFC'], ascending=[False, False])
        
        return df
    
    def _generate_narrative(self, 
                          selected_contrasts: List[ContrastSelection],
                          thresholds: ThresholdSelection,
                          top_genes: List[str]) -> str:
        """Generate interpretive narrative using LLM."""
        
        # Prepare context for narrative generation
        contrast_summaries = []
        for contrast in selected_contrasts:
            desc = self.integrator._get_contrast_description(contrast.analysis_id, contrast.contrast_id)
            contrast_summaries.append(f"- {contrast.contrast_id}: {desc} (Score: {contrast.relevance_score:.1f})")
        
        gene_list_str = ', '.join(top_genes[:10])
        if len(top_genes) > 10:
            gene_list_str += f' and {len(top_genes)-10} others'
        
        prompt = f"""
You are an expert computational biologist writing a clear, accessible summary of RNA-seq differential expression results.

Research Context: {self.biological_prompt}

Analysis Summary:
- {len(selected_contrasts)} biologically relevant contrasts were automatically selected
- Statistical thresholds: FDR < {thresholds.fdr_cutoff}, |log2FC| > {thresholds.logfc_cutoff}
- {len(top_genes)} key differentially expressed genes identified

Selected Contrasts:
{chr(10).join(contrast_summaries)}

Top Genes: {gene_list_str}

Threshold Justification: {thresholds.justification}

Write a 2-3 paragraph narrative that:
1. Summarizes the key findings in accessible language
2. Highlights the most important genes and patterns
3. Provides biological context and potential significance
4. Maintains scientific accuracy while being readable

Focus on biological insights rather than technical details.
"""
        
        try:
            response = self.client.chat.completions.create(
                model="gpt-4o-mini",
                messages=[{"role": "user", "content": prompt}],
                max_tokens=1000
            )
            
            return response.choices[0].message.content.strip()
            
        except Exception as e:
            logger.warning(f"Narrative generation failed: {e}, using fallback")
            return self._fallback_narrative(selected_contrasts, thresholds, top_genes)
    
    def _fallback_narrative(self, 
                           selected_contrasts: List[ContrastSelection],
                           thresholds: ThresholdSelection,
                           top_genes: List[str]) -> str:
        """Fallback narrative generation."""
        
        narrative = f"""
This analysis automatically identified {len(selected_contrasts)} biologically relevant experimental contrasts related to {self.biological_prompt}. Using statistically rigorous thresholds (FDR < {thresholds.fdr_cutoff}, |log2FC| > {thresholds.logfc_cutoff}), we discovered {len(top_genes)} key genes that show consistent differential expression patterns.

The most frequently dysregulated genes include {', '.join(top_genes[:5])}{'and others' if len(top_genes) > 5 else ''}. These genes appear across multiple experimental conditions, suggesting they may play central roles in the biological processes under investigation.

The results provide a foundation for understanding the molecular mechanisms underlying the observed phenotypes and suggest potential targets for further experimental validation.
"""
        
        return narrative
    
    def _generate_html(self, landing_data: LandingPageData, output_path: str):
        """Generate the final HTML landing page."""
        
        # Convert plotly figure to HTML
        heatmap_html = ""
        if landing_data.heatmap_fig:
            heatmap_html = landing_data.heatmap_fig.to_html(include_plotlyjs='cdn', div_id='heatmap')
        
        # Create contrast justification table
        contrast_rows = []
        for contrast in landing_data.selected_contrasts:
            contrast_rows.append(
                f"<tr><td>{contrast.analysis_id}</td><td>{contrast.contrast_id}</td>"
                f"<td>{contrast.relevance_score:.1f}</td><td>{contrast.justification}</td></tr>"
            )
        
        # Create gene table HTML
        gene_table_html = ""
        if not landing_data.gene_table.empty:
            gene_table_html = landing_data.gene_table.to_html(index=False, classes='gene-table', table_id='gene-table')
        
        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>UORCA Landing Page - AI-Assisted Interpretation</title>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <style>
        body {{ font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; line-height: 1.6; margin: 0; padding: 20px; color: #333; background-color: #f8f9fa; }}
        .container {{ max-width: 1200px; margin: 0 auto; background: white; padding: 30px; border-radius: 8px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }}
        h1 {{ color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }}
        h2 {{ color: #34495e; margin-top: 30px; }}
        .narrative {{ background: #e8f4fd; padding: 20px; border-radius: 5px; border-left: 4px solid #3498db; margin: 20px 0; }}
        .threshold-info {{ background: #f0f8e8; padding: 15px; border-radius: 5px; border-left: 4px solid #27ae60; margin: 20px 0; }}
        table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}
        th, td {{ border: 1px solid #ddd; padding: 12px; text-align: left; }}
        th {{ background-color: #f2f2f2; font-weight: 600; }}
        .gene-table {{ font-size: 14px; }}
        .timestamp {{ color: #7f8c8d; font-size: 0.9em; margin-top: 20px; }}
        .stats {{ display: flex; gap: 20px; margin: 20px 0; }}
        .stat-box {{ background: #ecf0f1; padding: 15px; border-radius: 5px; text-align: center; flex: 1; }}
        .stat-number {{ font-size: 2em; font-weight: bold; color: #2c3e50; }}
        .stat-label {{ color: #7f8c8d; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>🧬 UORCA AI-Assisted Interpretation</h1>
        <p><strong>Research Context:</strong> {self.biological_prompt}</p>
        
        <div class="stats">
            <div class="stat-box">
                <div class="stat-number">{len(landing_data.selected_contrasts)}</div>
                <div class="stat-label">Selected Contrasts</div>
            </div>
            <div class="stat-box">
                <div class="stat-number">{len(landing_data.top_genes)}</div>
                <div class="stat-label">Key Genes</div>
            </div>
            <div class="stat-box">
                <div class="stat-number">{landing_data.thresholds.fdr_cutoff}</div>
                <div class="stat-label">FDR Threshold</div>
            </div>
        </div>
        
        <div class="narrative">
            <h2>📝 Key Findings</h2>
            {landing_data.narrative}
        </div>
        
        <div class="threshold-info">
            <h3>🎯 Analysis Parameters</h3>
            <p><strong>Statistical Thresholds:</strong> FDR < {landing_data.thresholds.fdr_cutoff}, |log2FC| > {landing_data.thresholds.logfc_cutoff}, Min. frequency ≥ {landing_data.thresholds.min_frequency}</p>
            <p><strong>Justification:</strong> {landing_data.thresholds.justification}</p>
        </div>
        
        <h2>🔍 Selected Contrasts and Justifications</h2>
        <table>
            <tr>
                <th>Dataset</th>
                <th>Contrast</th>
                <th>Relevance Score</th>
                <th>Justification</th>
            </tr>
            {''.join(contrast_rows)}
        </table>
        
        <h2>🧬 Top Differentially Expressed Genes</h2>
        {gene_table_html}
        
        <h2>🌡️ Interactive Expression Heatmap</h2>
        <div id="heatmap-container">
            {heatmap_html}
        </div>
        
        <div class="timestamp">
            <p>Generated on {datetime.now().strftime('%Y-%m-%d at %H:%M:%S')} using UORCA AI-assisted interpretation</p>
        </div>
    </div>
</body>
</html>
"""
        
        # Write to file
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(html_content)

def main():
    """Command line interface for the landing page generator."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Generate UORCA AI-assisted landing page")
    parser.add_argument("--results_dir", required=True, help="Directory containing UORCA results")
    parser.add_argument("--output", required=True, help="Output HTML file path")
    parser.add_argument("--prompt", default="General differential expression analysis", 
                       help="Biological research context/prompt")
    parser.add_argument("--max_contrasts", type=int, default=8, 
                       help="Maximum number of contrasts to include")
    parser.add_argument("--max_genes", type=int, default=50, 
                       help="Maximum number of genes to include")
    
    args = parser.parse_args()
    
    # Set up logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    
    # Generate landing page
    generator = LandingPageGenerator(args.results_dir, args.prompt)
    landing_data = generator.generate_landing_page(
        max_contrasts=args.max_contrasts,
        max_genes=args.max_genes,
        output_path=args.output
    )
    
    print(f"Landing page generated successfully: {args.output}")
    print(f"Selected {len(landing_data.selected_contrasts)} contrasts and {len(landing_data.top_genes)} genes")

if __name__ == "__main__":
    main()
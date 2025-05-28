#!/usr/bin/env python3
"""
Integration script for adding landing page generation to UORCA master workflow.

This script provides a clean interface for the master agent to generate 
AI-assisted landing pages after analysis completion.
"""

import os
import logging
import json
from pathlib import Path
from typing import Optional, Dict, Any
from datetime import datetime

from landing_page_generator import LandingPageGenerator, LandingPageData

logger = logging.getLogger(__name__)

class MasterWorkflowIntegration:
    """
    Integration class for adding landing page generation to UORCA master workflow.
    """
    
    def __init__(self, output_dir: str, accession: str):
        """
        Initialize integration for a specific analysis run.
        
        Parameters:
        -----------
        output_dir : str
            Base output directory containing analysis results
        accession : str
            Dataset accession being analyzed
        """
        self.output_dir = output_dir
        self.accession = accession
        self.analysis_dir = os.path.join(output_dir, accession)
        
    def should_generate_landing_page(self) -> bool:
        """
        Check if landing page generation is appropriate for this analysis.
        
        Returns:
        --------
        bool : True if landing page should be generated
        """
        # Check if analysis was successful
        analysis_info_paths = [
            os.path.join(self.analysis_dir, "metadata", "analysis_info.json"),
            os.path.join(self.analysis_dir, "analysis_info.json")
        ]
        
        for info_path in analysis_info_paths:
            if os.path.exists(info_path):
                try:
                    with open(info_path, 'r') as f:
                        analysis_info = json.load(f)
                    
                    analysis_success = analysis_info.get('analysis_success', False)
                    if analysis_success:
                        logger.info(f"Analysis successful for {self.accession}, landing page generation recommended")
                        return True
                    else:
                        logger.info(f"Analysis not successful for {self.accession}, skipping landing page")
                        return False
                        
                except (json.JSONDecodeError, FileNotFoundError) as e:
                    logger.warning(f"Could not read analysis info: {e}")
        
        # Fallback: check if basic files exist
        rnaseq_dir = os.path.join(self.analysis_dir, "RNAseqAnalysis")
        if os.path.exists(rnaseq_dir):
            # Look for at least one DEG file
            for item in os.listdir(rnaseq_dir):
                deg_file = os.path.join(rnaseq_dir, item, "DEG.csv")
                if os.path.exists(deg_file):
                    logger.info(f"Found DEG files for {self.accession}, proceeding with landing page")
                    return True
        
        logger.info(f"No suitable results found for {self.accession}, skipping landing page")
        return False
    
    def extract_biological_context(self) -> str:
        """
        Extract biological context from dataset information for AI prompt.
        
        Returns:
        --------
        str : Biological context prompt for AI
        """
        # Try to get dataset information from various sources
        context_parts = []
        
        # 1. Check analysis_info.json for dataset_information
        analysis_info_paths = [
            os.path.join(self.analysis_dir, "metadata", "analysis_info.json"),
            os.path.join(self.analysis_dir, "analysis_info.json")
        ]
        
        for info_path in analysis_info_paths:
            if os.path.exists(info_path):
                try:
                    with open(info_path, 'r') as f:
                        analysis_info = json.load(f)
                    
                    dataset_info = analysis_info.get('dataset_information', '')
                    organism = analysis_info.get('organism', '')
                    accession = analysis_info.get('accession', self.accession)
                    
                    if dataset_info:
                        context_parts.append(f"Dataset {accession} ({organism}): {dataset_info}")
                        break
                        
                except (json.JSONDecodeError, FileNotFoundError):
                    continue
        
        # 2. Fallback: create basic context from accession and organism
        if not context_parts:
            # Try to get organism from analysis_info
            organism = "unknown organism"
            for info_path in analysis_info_paths:
                if os.path.exists(info_path):
                    try:
                        with open(info_path, 'r') as f:
                            analysis_info = json.load(f)
                        organism = analysis_info.get('organism', organism)
                        break
                    except:
                        continue
            
            context_parts.append(f"Differential expression analysis of {organism} RNA-seq data from dataset {self.accession}")
        
        # 3. Add contrast information if available
        contrasts_file = os.path.join(self.analysis_dir, "metadata", "contrasts.csv")
        if os.path.exists(contrasts_file):
            try:
                import pandas as pd
                contrasts_df = pd.read_csv(contrasts_file)
                if 'description' in contrasts_df.columns:
                    descriptions = contrasts_df['description'].tolist()[:3]  # First 3 descriptions
                    if descriptions and any(desc for desc in descriptions if desc and len(str(desc)) > 10):
                        context_parts.append(f"Experimental contrasts include: {'; '.join(str(d)[:100] for d in descriptions if str(d) != 'nan')}")
            except Exception as e:
                logger.debug(f"Could not extract contrast descriptions: {e}")
        
        # Combine into final prompt
        biological_prompt = ". ".join(context_parts)
        
        # Clean and truncate if too long
        biological_prompt = biological_prompt.replace('\n', ' ').replace('\r', ' ')
        if len(biological_prompt) > 1000:
            biological_prompt = biological_prompt[:1000] + "..."
        
        logger.info(f"Generated biological context: {biological_prompt[:200]}...")
        return biological_prompt
    
    def generate_landing_page(self, 
                            biological_context: Optional[str] = None,
                            max_contrasts: int = 8,
                            max_genes: int = 50) -> Optional[LandingPageData]:
        """
        Generate landing page for the analysis results.
        
        Parameters:
        -----------
        biological_context : str, optional
            Biological research context. If None, will be extracted automatically.
        max_contrasts : int
            Maximum number of contrasts to include
        max_genes : int
            Maximum number of genes to include
            
        Returns:
        --------
        LandingPageData or None : Generated landing page data or None if failed
        """
        try:
            # Extract biological context if not provided
            if biological_context is None:
                biological_context = self.extract_biological_context()
            
            logger.info(f"Generating landing page for {self.accession}")
            
            # Initialize landing page generator
            generator = LandingPageGenerator(
                results_dir=self.output_dir,  # Use parent directory containing all datasets
                biological_prompt=biological_context
            )
            
            # Create output path
            landing_pages_dir = os.path.join(self.analysis_dir, "landing_pages")
            os.makedirs(landing_pages_dir, exist_ok=True)
            
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            output_path = os.path.join(landing_pages_dir, f"landing_page_{timestamp}.html")
            
            # Generate landing page
            landing_data = generator.generate_landing_page(
                max_contrasts=max_contrasts,
                max_genes=max_genes,
                output_path=output_path
            )
            
            # Also create a "latest" symlink or copy
            latest_path = os.path.join(landing_pages_dir, "landing_page_latest.html")
            try:
                if os.path.exists(latest_path):
                    os.remove(latest_path)
                
                # Copy instead of symlink for better compatibility
                import shutil
                shutil.copy2(output_path, latest_path)
                
            except Exception as e:
                logger.warning(f"Could not create latest landing page link: {e}")
            
            logger.info(f"Landing page generated successfully: {output_path}")
            return landing_data
            
        except Exception as e:
            logger.error(f"Failed to generate landing page for {self.accession}: {str(e)}")
            return None

def add_landing_page_tool_to_master():
    """
    Function to add landing page generation tool to master agent.
    This should be called during master agent initialization.
    """
    from pydantic_ai import Agent
    from shared import RNAseqCoreContext
    from shared.workflow_logging import log_tool
    
    def create_landing_page_tool(master_agent: Agent):
        """Create the landing page tool for the master agent."""
        
        @master_agent.tool
        @log_tool
        async def generate_landing_page(ctx, 
                                      biological_context: Optional[str] = None,
                                      max_contrasts: int = 8,
                                      max_genes: int = 50) -> str:
            """
            Generate an AI-assisted landing page summarizing the analysis results.
            
            This tool creates an intelligent summary of differential expression results with:
            - Automatic contrast selection based on biological relevance  
            - AI-optimized statistical thresholds
            - Interactive visualizations with hover details
            - Interpretive narratives in accessible language
            
            Parameters:
            -----------
            biological_context : str, optional
                Research context to guide AI selections. If not provided, will be extracted
                from dataset metadata automatically.
            max_contrasts : int
                Maximum number of experimental contrasts to include (default: 8)
            max_genes : int  
                Maximum number of top genes to display (default: 50)
                
            Returns:
            --------
            str : Summary of landing page generation results
            """
            try:
                logger.info("🎨 Generating AI-assisted landing page")
                
                # Initialize integration
                integration = MasterWorkflowIntegration(
                    output_dir=ctx.deps.output_dir,
                    accession=ctx.deps.accession
                )
                
                # Check if landing page should be generated
                if not integration.should_generate_landing_page():
                    return "Landing page generation skipped - analysis not successful or suitable results not found."
                
                # Generate landing page
                landing_data = integration.generate_landing_page(
                    biological_context=biological_context,
                    max_contrasts=max_contrasts,
                    max_genes=max_genes
                )
                
                if landing_data is None:
                    return "Landing page generation failed due to errors in processing."
                
                # Return summary
                summary = f"""
Landing page generated successfully!

Key Results:
- Selected {len(landing_data.selected_contrasts)} most relevant experimental contrasts
- Identified {len(landing_data.top_genes)} top differentially expressed genes  
- Applied AI-optimized thresholds: FDR < {landing_data.thresholds.fdr_cutoff}, |logFC| > {landing_data.thresholds.logfc_cutoff}
- Generated interpretive narrative and interactive visualizations

The landing page provides:
✓ Automated contrast selection with AI justifications
✓ Optimal statistical threshold selection  
✓ Interactive heatmap with hover details
✓ Ranked gene table with frequency statistics
✓ Biological interpretation in accessible language

Access your landing page at:
{os.path.join(ctx.deps.output_dir, ctx.deps.accession, "landing_pages", "landing_page_latest.html")}
"""
                
                return summary
                
            except Exception as e:
                error_msg = f"Error generating landing page: {str(e)}"
                logger.error(error_msg)
                return error_msg
        
        return generate_landing_page
    
    return create_landing_page_tool

def integrate_with_master_workflow():
    """
    Main integration function to add landing page generation to master workflow.
    
    This function should be called in master.py to add the landing page tool.
    """
    tool_creator = add_landing_page_tool_to_master()
    logger.info("Landing page tool integration ready for master agent")
    return tool_creator

# Standalone usage for testing
def main():
    """Test the integration with a sample dataset."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Test landing page integration")
    parser.add_argument("--output_dir", required=True, help="UORCA results directory")
    parser.add_argument("--accession", required=True, help="Dataset accession to process")
    parser.add_argument("--context", help="Biological context (optional)")
    
    args = parser.parse_args()
    
    # Set up logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    
    # Test integration
    integration = MasterWorkflowIntegration(args.output_dir, args.accession)
    
    if integration.should_generate_landing_page():
        print(f"✅ {args.accession} is suitable for landing page generation")
        
        biological_context = args.context or integration.extract_biological_context()
        print(f"📝 Biological context: {biological_context[:200]}...")
        
        landing_data = integration.generate_landing_page(biological_context)
        
        if landing_data:
            print(f"🎉 Landing page generated successfully!")
            print(f"   - Contrasts: {len(landing_data.selected_contrasts)}")
            print(f"   - Genes: {len(landing_data.top_genes)}")
            print(f"   - Thresholds: FDR={landing_data.thresholds.fdr_cutoff}, FC={landing_data.thresholds.logfc_cutoff}")
        else:
            print(f"❌ Landing page generation failed")
    else:
        print(f"⏭️ {args.accession} not suitable for landing page generation")

if __name__ == "__main__":
    main()
#!/usr/bin/env python3
"""
Example usage of the UORCA Landing Page Generator

This script demonstrates how to use the landing page generator programmatically
with different biological contexts and parameters.
"""

import os
import sys
import logging
from pathlib import Path

# Add the parent directory to the path to import the generator
sys.path.append(str(Path(__file__).parent))
from landing_page_generator import LandingPageGenerator

def setup_logging():
    """Set up logging for the example script."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def example_neuroblastoma_analysis():
    """Example: Generate landing page for neuroblastoma research."""
    print("=" * 60)
    print("EXAMPLE 1: Neuroblastoma Research")
    print("=" * 60)
    
    # Example results directory (adjust path as needed)
    results_dir = "../../../UORCA_results"
    
    if not os.path.exists(results_dir):
        print(f"❌ Results directory not found: {results_dir}")
        return
    
    # Define biological context
    biological_prompt = """
    Neuroblastoma tumor development and progression. I am interested in understanding
    the molecular mechanisms driving neuroblastoma pathogenesis, including MYCN amplification
    effects, tumor vs normal tissue differences, and potential therapeutic targets.
    Key genes of interest include MYCN, ALK, PHOX2B, and p53 pathway components.
    """
    
    # Initialize generator
    generator = LandingPageGenerator(results_dir, biological_prompt)
    
    # Generate landing page
    output_path = "neuroblastoma_landing_page.html"
    landing_data = generator.generate_landing_page(
        max_contrasts=10,
        max_genes=60,
        output_path=output_path
    )
    
    print(f"✅ Generated landing page: {output_path}")
    print(f"   - Selected {len(landing_data.selected_contrasts)} contrasts")
    print(f"   - Identified {len(landing_data.top_genes)} key genes")
    print(f"   - Thresholds: FDR < {landing_data.thresholds.fdr_cutoff}, |logFC| > {landing_data.thresholds.logfc_cutoff}")
    
    # Print top genes
    print(f"\n🧬 Top 10 genes: {', '.join(landing_data.top_genes[:10])}")
    
    return landing_data

def example_cancer_immunology():
    """Example: Generate landing page for cancer immunology research."""
    print("\n" + "=" * 60)
    print("EXAMPLE 2: Cancer Immunology")
    print("=" * 60)
    
    results_dir = "../../../UORCA_results"
    
    if not os.path.exists(results_dir):
        print(f"❌ Results directory not found: {results_dir}")
        return
    
    biological_prompt = """
    Cancer immunology and immune response modulation. Focus on T cell activation,
    immune checkpoint pathways, tumor microenvironment, and immunotherapy responses.
    Interested in PD-1/PD-L1 axis, CTLA-4, interferon signaling, and immune cell infiltration markers.
    """
    
    generator = LandingPageGenerator(results_dir, biological_prompt)
    
    output_path = "immunology_landing_page.html"
    landing_data = generator.generate_landing_page(
        max_contrasts=6,
        max_genes=40,
        output_path=output_path
    )
    
    print(f"✅ Generated landing page: {output_path}")
    print(f"   - Selected {len(landing_data.selected_contrasts)} contrasts")
    print(f"   - Identified {len(landing_data.top_genes)} key genes")
    
    return landing_data

def example_minimal_analysis():
    """Example: Minimal analysis with default parameters."""
    print("\n" + "=" * 60)
    print("EXAMPLE 3: Minimal Analysis")
    print("=" * 60)
    
    results_dir = "../../../UORCA_results"
    
    if not os.path.exists(results_dir):
        print(f"❌ Results directory not found: {results_dir}")
        return
    
    # Simple prompt
    biological_prompt = "General differential expression analysis"
    
    generator = LandingPageGenerator(results_dir, biological_prompt)
    
    # Use default parameters
    output_path = "minimal_landing_page.html"
    landing_data = generator.generate_landing_page(output_path=output_path)
    
    print(f"✅ Generated minimal landing page: {output_path}")
    print(f"   - Selected {len(landing_data.selected_contrasts)} contrasts")
    print(f"   - Identified {len(landing_data.top_genes)} key genes")
    
    return landing_data

def example_custom_thresholds():
    """Example: How to access detailed results without generating HTML."""
    print("\n" + "=" * 60)
    print("EXAMPLE 4: Accessing Detailed Results")
    print("=" * 60)
    
    results_dir = "../../../UORCA_results"
    
    if not os.path.exists(results_dir):
        print(f"❌ Results directory not found: {results_dir}")
        return
    
    biological_prompt = "Developmental biology and cell differentiation"
    
    generator = LandingPageGenerator(results_dir, biological_prompt)
    
    # Generate without saving HTML
    landing_data = generator.generate_landing_page(
        max_contrasts=5,
        max_genes=30,
        output_path=None  # Don't save HTML
    )
    
    print("📊 Detailed Analysis Results:")
    print(f"   Research context: {generator.biological_prompt}")
    
    print("\n🔍 Selected Contrasts:")
    for i, contrast in enumerate(landing_data.selected_contrasts, 1):
        print(f"   {i}. {contrast.analysis_id}_{contrast.contrast_id}")
        print(f"      Score: {contrast.relevance_score:.1f}")
        print(f"      Justification: {contrast.justification[:100]}...")
    
    print(f"\n🎯 Automatic Threshold Selection:")
    print(f"   FDR cutoff: {landing_data.thresholds.fdr_cutoff}")
    print(f"   LogFC cutoff: {landing_data.thresholds.logfc_cutoff}")
    print(f"   Min frequency: {landing_data.thresholds.min_frequency}")
    print(f"   Justification: {landing_data.thresholds.justification}")
    
    print(f"\n🧬 Gene Table Shape: {landing_data.gene_table.shape}")
    if not landing_data.gene_table.empty:
        print("   Top 5 genes by frequency:")
        for _, row in landing_data.gene_table.head().iterrows():
            print(f"     - {row['Gene']}: freq={row['Frequency']}, maxFC={row['Max_LogFC']}")
    
    return landing_data

def compare_different_prompts():
    """Example: Compare results with different biological prompts."""
    print("\n" + "=" * 60)
    print("EXAMPLE 5: Comparing Different Research Contexts")
    print("=" * 60)
    
    results_dir = "../../../UORCA_results"
    
    if not os.path.exists(results_dir):
        print(f"❌ Results directory not found: {results_dir}")
        return
    
    prompts = [
        ("metabolism", "Metabolic pathways and energy production"),
        ("apoptosis", "Cell death and apoptosis signaling"),
        ("development", "Developmental biology and morphogenesis")
    ]
    
    results = {}
    
    for name, prompt in prompts:
        print(f"\n🔬 Analyzing: {name}")
        generator = LandingPageGenerator(results_dir, prompt)
        landing_data = generator.generate_landing_page(
            max_contrasts=4,
            max_genes=20,
            output_path=None
        )
        
        results[name] = {
            'contrasts': len(landing_data.selected_contrasts),
            'genes': len(landing_data.top_genes),
            'fdr_threshold': landing_data.thresholds.fdr_cutoff,
            'top_genes': landing_data.top_genes[:5]
        }
        
        print(f"   Contrasts: {results[name]['contrasts']}")
        print(f"   Genes: {results[name]['genes']}")
        print(f"   Top genes: {', '.join(results[name]['top_genes'])}")
    
    print(f"\n📊 Summary Comparison:")
    for name, data in results.items():
        print(f"   {name:12}: {data['contrasts']} contrasts, {data['genes']} genes, FDR={data['fdr_threshold']}")

def main():
    """Run all examples."""
    setup_logging()
    
    print("🧬 UORCA Landing Page Generator - Example Usage")
    print("This script demonstrates various ways to use the landing page generator.")
    
    # Check if OpenAI API key is set
    if not os.getenv("OPENAI_API_KEY"):
        print("⚠️ Warning: OPENAI_API_KEY environment variable not set.")
        print("   AI features will fall back to heuristic methods.")
    
    try:
        # Run examples
        example_neuroblastoma_analysis()
        example_cancer_immunology()
        example_minimal_analysis()
        example_custom_thresholds()
        compare_different_prompts()
        
        print("\n" + "=" * 60)
        print("✅ All examples completed successfully!")
        print("Check the generated HTML files for interactive landing pages.")
        print("=" * 60)
        
    except Exception as e:
        print(f"\n❌ Error running examples: {str(e)}")
        logging.exception("Error in example script")

if __name__ == "__main__":
    main()
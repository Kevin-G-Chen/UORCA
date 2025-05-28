# UORCA Landing Page Generator

**AI-Assisted Interpretation with Minimal Interactivity**

Generate intelligent, self-contained summaries of your RNA-seq differential expression results with minimal user input. The landing page automatically selects the most relevant contrasts, determines optimal statistical thresholds, and creates interpretive narratives using AI.

## 🎯 Overview

The UORCA Landing Page Generator implements the design goals outlined in the [landing page design document](../../../scratch/Notes/landing_page_design.md):

- **Automated Analysis**: AI selects contrasts and thresholds based on biological relevance
- **Minimal Interaction**: Users provide only a research context prompt
- **Intelligent Interpretation**: LLM-generated narratives explain findings in accessible language
- **Interactive Visualizations**: Hover-enabled heatmaps and gene tables
- **Reproducible Results**: All selections are documented with AI justifications

## 🚀 Quick Start

### Prerequisites

- UORCA analysis results directory
- OpenAI API key (set as environment variable `OPENAI_API_KEY`)
- Python packages: `pandas`, `plotly`, `openai`, `streamlit` (optional)

### Command Line Usage

```bash
# Basic usage
python landing_page_generator.py \
    --results_dir /path/to/UORCA_results \
    --output landing_page.html \
    --prompt "neuroblastoma tumor vs normal tissue"

# Advanced usage
python landing_page_generator.py \
    --results_dir /path/to/UORCA_results \
    --output detailed_analysis.html \
    --prompt "immune checkpoint inhibitor response in melanoma" \
    --max_contrasts 12 \
    --max_genes 75
```

### Streamlit Interface

```bash
# Launch interactive web interface
streamlit run landing_page_app.py

# Access at http://localhost:8501
```

### Programmatic Usage

```python
from landing_page_generator import LandingPageGenerator

# Initialize with research context
generator = LandingPageGenerator(
    results_dir="/path/to/UORCA_results",
    biological_prompt="MYCN amplification in neuroblastoma"
)

# Generate landing page
landing_data = generator.generate_landing_page(
    max_contrasts=8,
    max_genes=50,
    output_path="mycn_analysis.html"
)

# Access results programmatically
print(f"Selected {len(landing_data.selected_contrasts)} contrasts")
print(f"Top genes: {', '.join(landing_data.top_genes[:10])}")
```

## 🧠 How It Works

### 1. Automatic Contrast Selection

- **AI Scoring**: GPT-4o-mini evaluates all available contrasts for biological relevance
- **Hybrid Approach**: Combines embedding similarity with LLM re-ranking
- **Justification**: Each selected contrast includes AI-generated rationale
- **Fallback**: Heuristic scoring based on DEG counts if AI unavailable

### 2. Intelligent Threshold Selection

- **Context-Aware**: Considers dataset size, sample numbers, and research context
- **Dynamic**: FDR and LogFC thresholds adapted to study characteristics
- **Justified**: AI provides scientific reasoning for threshold choices
- **Robust**: Fallback heuristics ensure reliable operation

### 3. Gene Aggregation & Ranking

- **Frequency Scoring**: Genes ranked by appearance across contrasts
- **Effect Size**: Weighted by maximum absolute log fold change
- **Statistical Rigor**: Only genes meeting threshold criteria included
- **Comprehensive**: Tracks gene behavior across all selected contrasts

### 4. Narrative Generation

- **Accessible Language**: Technical findings translated to readable summaries
- **Biological Context**: Interpretations grounded in research question
- **Key Insights**: Highlights most important genes and patterns
- **Scientific Accuracy**: Maintains precision while improving accessibility

## 📊 Output Components

### Landing Page Sections

1. **Summary Metrics**: Key statistics and thresholds
2. **Key Findings**: AI-generated narrative interpretation
3. **Analysis Parameters**: Threshold justifications and methodology
4. **Selected Contrasts**: Table with relevance scores and AI justifications
5. **Top Genes**: Ranked table with frequency and effect size data
6. **Interactive Heatmap**: Clustered visualization with hover details

### Data Objects

```python
@dataclass
class LandingPageData:
    selected_contrasts: List[ContrastSelection]  # AI-selected contrasts
    thresholds: ThresholdSelection              # Optimal thresholds
    top_genes: List[str]                        # Ranked gene list
    heatmap_fig: go.Figure                      # Interactive plotly figure
    gene_table: pd.DataFrame                    # Summary statistics
    narrative: str                              # AI-generated interpretation
```

## 🔧 Configuration Options

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `max_contrasts` | 8 | Maximum number of contrasts to include |
| `max_genes` | 50 | Maximum number of genes in visualizations |
| `biological_prompt` | Required | Research context for AI guidance |
| `results_dir` | Required | Path to UORCA analysis results |

### Environment Variables

- `OPENAI_API_KEY`: Required for AI features
- `OPENAI_MODEL`: Optional, defaults to "gpt-4o-mini"

## 📁 Required Directory Structure

Your UORCA results directory should contain:

```
UORCA_results/
├── dataset1/
│   ├── RNAseqAnalysis/
│   │   ├── contrast1/
│   │   │   └── DEG.csv
│   │   ├── contrast2/
│   │   │   └── DEG.csv
│   │   └── CPM.csv
│   └── metadata/
│       ├── contrasts.csv
│       └── analysis_info.json
├── dataset2/
│   └── ... (similar structure)
```

### Required Files

- **DEG.csv**: Differential expression results with columns `Gene`, `logFC`, `adj.P.Val`
- **contrasts.csv**: Contrast descriptions with columns `name`, `description`
- **analysis_info.json**: Analysis metadata including organism and sample counts

## 🎮 Examples

### Research Context Examples

```python
# Oncology research
prompt = """
MYCN amplification effects in neuroblastoma. Interested in understanding
how MYCN drives tumor progression and identifying potential therapeutic targets.
"""

# Immunology research  
prompt = """
T cell activation and immune checkpoint pathways. Focus on PD-1/PD-L1 axis
and response to immunotherapy in cancer patients.
"""

# Developmental biology
prompt = """
Neural development and differentiation. Studying transcriptional programs
controlling neuronal maturation and synaptic formation.
"""

# Metabolism research
prompt = """
Metabolic reprogramming in cancer cells. Understanding how tumor cells
alter glucose and lipid metabolism for growth and survival.
"""
```

### Batch Processing

```python
# Process multiple research contexts
contexts = [
    ("oncology", "tumor progression and metastasis"),
    ("immunology", "immune response modulation"),
    ("metabolism", "cellular energy production")
]

for name, prompt in contexts:
    generator = LandingPageGenerator(results_dir, prompt)
    landing_data = generator.generate_landing_page(
        output_path=f"{name}_landing_page.html"
    )
    print(f"Generated {name}: {len(landing_data.top_genes)} genes")
```

## 🛠️ Advanced Usage

### Custom AI Models

```python
# Use different OpenAI models
import os
os.environ["OPENAI_MODEL"] = "gpt-4o"  # More powerful but slower

generator = LandingPageGenerator(results_dir, prompt)
```

### Accessing Detailed Results

```python
# Generate without saving HTML
landing_data = generator.generate_landing_page(output_path=None)

# Access AI justifications
for contrast in landing_data.selected_contrasts:
    print(f"{contrast.contrast_id}: {contrast.justification}")

# Examine threshold selection
print(f"AI selected FDR={landing_data.thresholds.fdr_cutoff}")
print(f"Reasoning: {landing_data.thresholds.justification}")

# Export gene rankings
gene_df = landing_data.gene_table
gene_df.to_csv("ranked_genes.csv", index=False)
```

### Integration with Existing Workflows

```python
# Combine with UORCA Explorer
from ..uorca_explorer import ResultsIntegrator

integrator = ResultsIntegrator(results_dir)
generator = LandingPageGenerator(results_dir, prompt)

# Use landing page selections in explorer
landing_data = generator.generate_landing_page(output_path=None)
selected_genes = landing_data.top_genes

# Create detailed visualizations
heatmap = integrator.create_lfc_heatmap(genes=selected_genes)
expression_plots = integrator.create_expression_plots(genes=selected_genes)
```

## 🐛 Troubleshooting

### Common Issues

**No contrasts selected**
- Check that DEG.csv files exist in contrast subdirectories
- Verify contrast descriptions are available in contrasts.csv
- Try a more general biological prompt

**AI features not working**
- Ensure `OPENAI_API_KEY` environment variable is set
- Check OpenAI API quota and billing status
- Fallback heuristics will be used automatically

**Empty gene table**
- Lower statistical thresholds may be too strict
- Check that DEG.csv files have required columns
- Verify data quality in source files

**Heatmap not displaying**
- Ensure genes meet significance criteria
- Check that multiple contrasts are selected
- Verify plotly installation and compatibility

### Debug Mode

```python
import logging
logging.basicConfig(level=logging.DEBUG)

generator = LandingPageGenerator(results_dir, prompt)
# Detailed logging will show AI API calls and fallback usage
```

## 📚 API Reference

### Classes

#### `LandingPageGenerator`

Main class for generating AI-assisted landing pages.

**Methods:**
- `__init__(results_dir, biological_prompt)`: Initialize generator
- `generate_landing_page(max_contrasts, max_genes, output_path)`: Generate complete landing page
- `_select_relevant_contrasts(max_contrasts)`: AI-powered contrast selection
- `_select_optimal_thresholds(selected_contrasts)`: Intelligent threshold selection
- `_aggregate_and_rank_genes(contrasts, thresholds, max_genes)`: Gene ranking
- `_generate_narrative(contrasts, thresholds, genes)`: AI narrative generation

#### Data Classes

- `ContrastSelection`: Selected contrast with AI justification
- `ThresholdSelection`: Optimal thresholds with reasoning
- `LandingPageData`: Complete landing page results

## 🔄 Integration

### With UORCA Pipeline

Add landing page generation to your analysis workflow:

```python
# After running UORCA analysis
from main_workflow.reporting.landing_page import LandingPageGenerator

# Generate landing page for each research context
contexts = ["tumor progression", "immune response", "drug resistance"]

for context in contexts:
    generator = LandingPageGenerator(output_dir, context)
    generator.generate_landing_page(
        output_path=f"{output_dir}/landing_pages/{context.replace(' ', '_')}.html"
    )
```

### With Streamlit Apps

```python
# In your Streamlit app
from landing_page_generator import LandingPageGenerator

if st.button("Generate AI Summary"):
    generator = LandingPageGenerator(results_dir, user_prompt)
    landing_data = generator.generate_landing_page(output_path=None)
    
    st.markdown(landing_data.narrative)
    st.plotly_chart(landing_data.heatmap_fig)
    st.dataframe(landing_data.gene_table)
```

## 🤝 Contributing

### Adding New Features

1. **Custom Scoring Methods**: Extend `_score_contrasts_with_llm()` for domain-specific scoring
2. **Visualization Options**: Add new plot types in `_create_landing_heatmap()`
3. **Export Formats**: Extend `_generate_html()` for PDF, PowerPoint, etc.
4. **AI Models**: Support for different LLMs in threshold and narrative generation

### Testing

```bash
# Run example usage script
python example_usage.py

# Test with different research contexts
python -c "
from landing_page_generator import LandingPageGenerator
generator = LandingPageGenerator('./test_data', 'test prompt')
result = generator.generate_landing_page(output_path=None)
print(f'Success: {len(result.top_genes)} genes identified')
"
```

## 📄 License

Part of the UORCA (Unified Omics Reference Corpus of Analyses) project. See main repository for license details.

## 🆘 Support

For issues, feature requests, or questions:

1. Check troubleshooting section above
2. Review example usage in `example_usage.py`
3. Test with the Streamlit interface for debugging
4. Consult the main UORCA documentation

---

*Generated with ❤️ by the UORCA AI-Assisted Landing Page Generator*
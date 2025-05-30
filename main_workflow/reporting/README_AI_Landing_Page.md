# UORCA Explorer - AI Landing Page Integration

**AI-Assisted RNA-seq Interpretation with Minimal User Interaction**

The AI Landing Page is now seamlessly integrated into UORCA Explorer as the **first tab**, providing researchers with intelligent, automated analysis of their RNA-seq results. Generate publication-ready insights with just a biological research prompt.

## 🎯 Overview

The AI Landing Page transforms how researchers interact with RNA-seq differential expression results by:

- **Automatically selecting** the most biologically relevant contrasts from your data
- **Optimizing statistical thresholds** based on dataset characteristics and research context
- **Ranking genes** by biological importance across multiple experimental conditions
- **Generating interpretive narratives** in accessible, publication-ready language
- **Creating interactive visualizations** with detailed hover information and clustering

**No manual parameter tuning required** - just describe your research and let AI handle the complexity.

## 🚀 Quick Start

### 1. Launch UORCA Explorer

```bash
# Navigate to the reporting directory
cd UORCA/main_workflow/reporting

# Launch with the integrated launcher (recommended)
python launch_uorca_explorer.py

# Or launch directly with uv
uv run streamlit run uorca_explorer.py --server.port 8501
```

### 2. Load Your Data

1. Enter your UORCA results directory path in the sidebar
2. Wait for data loading to complete (you'll see status updates)
3. Navigate to the **"🎨 Landing Page"** tab (first tab)

### 3. Generate AI Analysis

1. **Enter Research Context**: Describe your biological question
   ```
   Example: "MYCN amplification effects in neuroblastoma tumors"
   ```

2. **Adjust Options** (optional): Set maximum contrasts and genes

3. **Generate**: Click "🚀 Generate AI Landing Page"

4. **Review Results**: 
   - AI-generated narrative summary displayed side-by-side with interactive heatmap
   - Automatically selected contrasts with justifications
   - Top differentially expressed genes with easy export options
   - Interactive expression heatmap with improved gene label visibility

### 4. Explore and Export

- Apply AI selections to other tabs for detailed analysis
- Export gene lists, summaries, and visualizations
- Use interactive features for deeper exploration

## ✨ Key Features

### 🤖 AI-Powered Analysis

| Feature | Description | Benefit |
|---------|-------------|---------|
| **Smart Contrast Selection** | AI evaluates all contrasts for biological relevance | No manual filtering needed |
| **Optimal Thresholds** | Automatically determines best statistical cutoffs | Scientifically rigorous results |
| **Gene Prioritization** | Ranks genes by frequency and effect size | Focus on most important findings |
| **Biological Narratives** | Generates accessible research summaries | Publication-ready interpretations |
| **Side-by-Side Display** | View AI interpretation alongside interactive heatmap | Enhanced cross-referencing and analysis |
| **Improved Visualization** | Better gene label visibility with adjustable heatmap height | Clear, readable gene identification |

### 🔄 Seamless Integration

- **First Tab Experience**: AI analysis is the first thing you see
- **Unified Interface**: No switching between different applications  
- **Selection Transfer**: Apply AI choices to detailed exploration tabs
- **Consistent Data**: Same robust data loading as rest of UORCA Explorer
- **Flexible Layout**: Choose between side-by-side or stacked display modes
- **Interactive Controls**: Adjustable heatmap height for optimal gene label visibility

### 🎯 Intelligent Automation

- **Context-Aware**: Adapts analysis to your specific research focus
- **Statistically Rigorous**: Uses appropriate thresholds for your dataset
- **Biologically Relevant**: Prioritizes meaningful contrasts and genes
- **Reproducible**: All selections documented with AI justifications
- **User-Friendly Interface**: Side-by-side layout for seamless interpretation
- **Customizable Display**: Adjustable heatmap height and layout preferences

## 📋 Requirements

### Essential Dependencies

```bash
# Managed automatically with uv
streamlit >= 1.45.1
plotly >= 6.1.1  
pandas >= 2.0.0
numpy >= 1.24.0
```

### AI Features (Optional)

```bash
# For full AI capabilities
openai >= 1.70.0

# Environment variable
export OPENAI_API_KEY=your_api_key_here
```

### Data Requirements

Your UORCA results directory should contain:

```
UORCA_results/
├── dataset1/
│   ├── RNAseqAnalysis/
│   │   ├── contrast1/
│   │   │   └── DEG.csv                    # Required: Gene, logFC, adj.P.Val columns
│   │   ├── contrast2/
│   │   │   └── DEG.csv
│   │   └── CPM.csv                        # Required: Gene column + sample columns
│   └── metadata/
│       ├── contrasts.csv                  # Required: name, description columns
│       ├── analysis_info.json             # Required: organism, sample counts
│       └── edger_analysis_samples.csv     # Optional: for expression plots
```

## 🎨 Using the AI Landing Page

### Tab Structure

UORCA Explorer now has **7 tabs**:

1. **🤖 View AI Summary** *(NEW)* - AI-assisted analysis and interpretation with side-by-side layout
2. **☑️ Select Data & Contrasts** - Manual dataset and contrast selection  
3. **🌡️ Explore Heatmap** - Interactive differential expression heatmaps
4. **📈 Plot Gene Expression** - Gene expression plots across samples
5. **🧑‍🔬 Analyze Experiments** - Quality control and diagnostic plots
6. **📋 View Dataset Info** - Study metadata and descriptions
7. **🔍 View Contrast Info** - Detailed contrast information

### Display Modes

The AI Landing Page offers two layout options:

#### **Side-by-Side Mode** (Default)
- **Left Panel**: AI-generated narrative in a scrollable container
- **Right Panel**: Interactive heatmap with adjustable height
- **Collapsible Sections**: Detailed contrast and gene tables expand on demand
- **Real-time Cross-Reference**: Compare AI insights with visualization simultaneously

#### **Stacked Mode** (Traditional)
- **Top Section**: Full-width AI narrative
- **Tabbed Results**: Separate tabs for contrasts, genes, and heatmap
- **Comprehensive View**: Each section gets full screen space

### AI Operation Modes

#### 🎯 Full AI Mode (with API key)
- **Smart Contrast Selection**: GPT-4o-mini evaluates biological relevance
- **Optimal Thresholds**: AI determines best statistical parameters
- **Rich Narratives**: Detailed biological interpretations with side-by-side display
- **Justifications**: Explains reasoning for each selection
- **Enhanced Visualization**: Improved heatmap with better gene label visibility

#### 🔧 Heuristic Mode (no API key)  
- **Basic Selection**: Contrast ranking by DEG counts and informativeness
- **Standard Thresholds**: Statistically sound fallback parameters
- **Template Narratives**: Structured summaries with key findings
- **Transparent Logic**: Clear explanations of selection criteria
- **Improved Layout**: Same side-by-side display and visualization enhancements

### Research Context Examples

#### Oncology Research
```
MYCN amplification effects in neuroblastoma. Focus on tumor 
progression, metastasis markers, and potential therapeutic 
targets including ALK pathway components.
```

#### Immunology Studies
```
T cell activation and immune checkpoint pathways. Interested 
in PD-1/PD-L1 responses, CAR-T cell effectiveness, and 
interferon signaling cascades.
```

#### Developmental Biology
```
Neural development and differentiation. Studying transcriptional 
programs controlling neuronal maturation, synaptic formation, 
and axon guidance mechanisms.
```

#### Drug Discovery
```
Compound mechanism of action and drug resistance. Focus on 
kinase inhibitor responses, metabolic reprogramming, and 
resistance pathway activation.
```

## 🛠️ Installation & Setup

### Method 1: Using UV (Recommended)

```bash
# Clone or navigate to UORCA project
cd UORCA/main_workflow/reporting

# Verify UV is available (should be installed in your project)
uv --version

# Check AI readiness
uv run python verify_ai_ready.py

# Launch UORCA Explorer
python launch_uorca_explorer.py
```

### Method 2: Manual Setup

```bash
# Install dependencies manually if needed
pip install streamlit plotly pandas numpy

# For AI features
pip install openai

# Set API key (choose one):
export OPENAI_API_KEY=your_key_here
# OR create .env file with: OPENAI_API_KEY=your_key_here

# Launch directly
streamlit run uorca_explorer.py --server.port 8501
```

### Environment Configuration

**Option 1: Environment Variables**
```bash
# Optional: Set default results directory
export UORCA_DEFAULT_RESULTS_DIR=/path/to/your/UORCA_results

# Required for full AI features
export OPENAI_API_KEY=your_openai_api_key

# Optional: Use different OpenAI model
export OPENAI_MODEL=gpt-4o  # Default: gpt-4o-mini
```

**Option 2: .env File (Recommended)**
Create a `.env` file in the project root:
```
# Required for full AI features
OPENAI_API_KEY=your_openai_api_key

# Optional settings
UORCA_DEFAULT_RESULTS_DIR=/path/to/your/UORCA_results
OPENAI_MODEL=gpt-4o-mini
```

The application automatically loads `.env` files from the project root directory.

## 🔧 Troubleshooting

### AI Landing Page Not Available

**Issue**: Tab shows "AI Landing Page: Not Available"

**Solutions**:
```bash
# Check if running with uv
uv run streamlit run uorca_explorer.py --server.port 8501

# Verify OpenAI package
uv run python -c "import openai; print('OpenAI available')"

# Check installation
python verify_ai_ready.py
```

### Heuristic Mode Only

**Issue**: AI features working but not using OpenAI

**Solutions**:
```bash
# Option 1: Set API key as environment variable
export OPENAI_API_KEY=your_key_here

# Option 2: Create .env file in project root
echo "OPENAI_API_KEY=your_key_here" > .env

# Verify API key is accessible
echo $OPENAI_API_KEY

# Check quota and billing at platform.openai.com
```

### No Contrasts Selected

**Issue**: "No suitable contrasts found" error

**Solutions**:
- **Check data structure**: Ensure DEG.csv files exist in contrast subdirectories
- **Verify file format**: Confirm required columns (Gene, logFC, adj.P.Val)
- **Broaden context**: Use more general biological descriptions
- **Check file permissions**: Ensure read access to all result files

### Empty Gene Table

**Issue**: No genes appear in results table

**Solutions**:
- **Lower thresholds**: AI may be selecting overly strict parameters
- **Check data quality**: Verify DEG.csv files contain significant genes
- **Increase max genes**: Set higher limit in options
- **Examine source data**: Look at individual DEG files for content

### Performance Issues

**Issue**: Slow loading or generation

**Solutions**:
- **Reduce scope**: Lower max_contrasts and max_genes settings
- **Check data size**: Very large datasets may need more processing time
- **Optimize browser**: Close unnecessary tabs and applications
- **Internet connection**: AI features require stable internet for OpenAI API

## 📊 Understanding Results

### Contrast Selection Criteria

The AI evaluates contrasts based on:

- **Biological Relevance**: Alignment with research context
- **Statistical Power**: Number of differentially expressed genes
- **Experimental Design**: Clarity and interpretability
- **Scientific Interest**: Novelty and potential impact

### Gene Ranking Algorithm

```python
# Genes ranked by combined score
rank_score = frequency_across_contrasts × log10(max_absolute_logFC + 1)
```

- **Frequency**: How many contrasts show the gene as significant
- **Effect Size**: Maximum absolute log2 fold change observed
- **Minimum Frequency**: Genes must appear in multiple contrasts

### Statistical Thresholds

AI considers:

- **Dataset Size**: Stricter thresholds for smaller studies
- **Sample Numbers**: Adjustment based on statistical power
- **Research Context**: Biology-specific threshold optimization
- **Multiple Testing**: Appropriate FDR control

## 💡 Advanced Usage

### Programmatic Access

```python
# Access landing page data after generation
if hasattr(st.session_state, 'landing_data'):
    landing_data = st.session_state.landing_data
    
    # Get AI-selected genes
    top_genes = landing_data.top_genes
    
    # Get AI narrative
    interpretation = landing_data.narrative
    
    # Get selected contrasts with justifications
    for contrast in landing_data.selected_contrasts:
        print(f"{contrast.contrast_id}: {contrast.justification}")
```

### Custom Research Contexts

```python
# Multi-modal research prompt
research_context = """
Integrated analysis of cancer progression and immune response:

Primary Focus:
- Epithelial-mesenchymal transition (EMT) markers
- Immune infiltration patterns and checkpoint expression

Secondary Analysis:
- Angiogenesis and vascular remodeling
- Drug resistance mechanisms and metabolic reprogramming

Model Systems:
- Patient-derived xenografts and organoids
- Co-culture systems with immune cells

Therapeutic Relevance:
- Kinase inhibitor combinations
- Immunotherapy response prediction
"""
```

### Workflow Integration

```python
# Typical research workflow
def research_workflow():
    # 1. Start with AI Landing Page (side-by-side mode)
    ai_summary = generate_ai_landing_page(
        integrator=data_integrator,
        biological_prompt=research_context,
        max_genes=50
    )
    
    # 2. Review narrative and heatmap simultaneously
    narrative = ai_summary.narrative
    heatmap = ai_summary.heatmap_fig  # Enhanced with better gene labels
    
    # 3. Apply AI selections to detailed analysis
    selected_genes = ai_summary.top_genes
    selected_contrasts = ai_summary.selected_contrasts
    
    # 4. Generate detailed visualizations
    detailed_heatmap = create_detailed_heatmap(selected_genes, selected_contrasts)
    expression_plots = create_expression_analysis(selected_genes)
    
    # 5. Export for further analysis
    export_gene_lists(selected_genes)
    export_narrative(ai_summary.narrative)
```

### Integration with External Tools

```python
# Export for pathway analysis
gene_list = landing_data.top_genes
with open("genes_for_pathway_analysis.txt", "w") as f:
    f.write("\n".join(gene_list))

# Export for publication
narrative_sections = {
    "methods": f"Statistical thresholds: FDR < {landing_data.thresholds.fdr_cutoff}",
    "results": landing_data.narrative,
    "contrasts": [c.justification for c in landing_data.selected_contrasts]
}
```

## 📈 Performance Tips

### Optimal Settings

- **Research Context**: Be specific but not overly narrow
- **Layout Mode**: Use side-by-side for cross-referencing, stacked for detailed review
- **Max Genes**: 30-50 for targeted analysis (optimal for heatmap visibility)
- **Heatmap Height**: Adjust slider for comfortable gene label reading
- **Browser**: Use Chrome, Firefox, or Safari for best performance

### Data Organization

- **Consistent Naming**: Use clear, descriptive contrast names
- **Rich Metadata**: Include detailed contrast descriptions
- **Quality Control**: Ensure all DEG files have required columns
- **Documentation**: Maintain analysis_info.json files

### API Usage

- **Batch Requests**: Process multiple datasets sequentially
- **Rate Limits**: OpenAI allows generous limits for typical usage
- **Fallback Ready**: System automatically uses heuristics if API unavailable
- **Cost Optimization**: gpt-4o-mini is cost-effective for most analyses

## 🔮 Future Enhancements

### Planned Features

- **Multi-species Support**: Enhanced organism-specific analysis
- **Pathway Integration**: Automatic pathway enrichment analysis  
- **Publication Export**: Direct export to manuscript formats
- **Collaborative Features**: Shared analysis sessions
- **Custom Models**: Support for domain-specific AI models

### Extensibility
## 📁 Simplified File Structure

After cleanup, the main files are:

```
reporting/
├── uorca_explorer.py           # Main integrated app with AI landing page
├── launch_uorca_explorer.py    # Launch script with dependency checking
├── ResultsIntegration.py       # Data loading and visualization infrastructure
├── verify_ai_ready.py         # Simple verification script
└── README_AI_Landing_Page.md   # This documentation
```

## 🔮 Future Enhancements

The integrated AI system is designed for easy extension:

- **New AI Models**: Add support for other LLM providers
- **Custom Scoring**: Implement domain-specific relevance scoring
- **Additional Visualizations**: Extend plot types and interactions
- **Export Formats**: Add support for new output formats

## 📞 Support & Resources

### Getting Help

1. **Check Logs**: Look at browser console for detailed error messages
2. **Verify Setup**: Run `python verify_ai_ready.py` for diagnostic info
3. **Test with Sample Data**: Ensure basic functionality works
4. **API Status**: Check OpenAI API status and billing if using full AI mode

### Documentation Links

- **UORCA Main Documentation**: `../../README.md`
- **Streamlit Documentation**: https://docs.streamlit.io/
- **OpenAI API Documentation**: https://platform.openai.com/docs/
- **Plotly Documentation**: https://plotly.com/python/

### Common Solutions

| Problem | Solution |
|---------|----------|
| "Module not found" errors | Use `uv run` instead of plain `python` |
| API rate limits | Wait and retry, or use heuristic mode |
| Slow performance | Reduce max_genes and max_contrasts |
| Empty results | Check data format and column names |
| Browser issues | Try incognito mode or different browser |

---

**The AI Landing Page makes RNA-seq analysis accessible to researchers at all levels - start with intelligent automation, then dive deep with powerful exploration tools.**

🧬 **Ready to explore your data?** Launch UORCA Explorer and let AI guide your discovery!
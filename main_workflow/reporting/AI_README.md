# UORCA AI Assistant

The UORCA AI Assistant is an interactive, AI-powered interface for exploring your RNA-seq analysis results. It uses a Model Context Protocol (MCP) server to provide LLMs with access to your UORCA analysis data, enabling natural language queries about gene expression, contrasts, and experimental results.

## 🚀 Features

- **Natural Language Queries**: Ask questions about your data in plain English
- **Interactive Data Explorer**: Point-and-click interface for gene-contrast combinations
- **Quick Queries**: Pre-built templates for common analysis questions
- **Chat Interface**: Free-form conversation with your analysis results
- **Real-time Analysis**: Get log fold changes, gene lists, and statistical insights instantly

## 📋 Requirements

### Dependencies
All required dependencies are included in the main UORCA environment:
- `pydantic-ai[mcp]` - AI agent framework with MCP support
- `openai` - OpenAI API client
- `streamlit` - Web interface
- `toml` - Configuration file parsing
- `mcp` - Model Context Protocol

### API Key
You need an OpenAI API key to use the AI assistant:
```bash
export OPENAI_API_KEY="your-api-key-here"
```

Or add it to your `.env` file in the project root:
```
OPENAI_API_KEY=your-api-key-here
```

## 🔧 Setup

### 1. Configuration
The AI assistant is automatically configured when you have:
- Valid UORCA analysis results in your results directory
- OpenAI API key set as environment variable
- All dependencies installed (included in UORCA environment)

### 2. Server Configuration
The MCP server configuration is defined in `servers.toml`:
```toml
[servers.uorca_data]
enabled = true
path = "main_workflow/reporting/mcp_servers/uorca_data_server.py"
```

### 3. Testing Setup
Run the validation script to check your setup:
```bash
uv run python test_ai_setup.py
```

## 🎯 Usage

### Accessing the AI Assistant
1. Start UORCA Explorer normally
2. Navigate to the "🤖 AI Assistant" tab
3. Wait for initialization (first time may take a few seconds)
4. Start exploring your data!

### Quick Queries Tab
Pre-built queries for common questions:
- **Overview Queries**: Get counts of contrasts, genes, and datasets
- **Specific Queries**: Look up individual genes or find top regulated genes

Example queries:
- "How many contrasts do I have?"
- "Show me the most upregulated genes"
- "Find gene TP53 across all contrasts"

### Data Explorer Tab
Interactive interface for detailed exploration:
1. **Select a contrast** from the dropdown
2. **Choose a gene** either from the list or enter manually
3. **Analyze** to get detailed information including:
   - Log fold change values
   - Statistical significance
   - Biological interpretation

**Batch Analysis**: Enter multiple genes (one per line) to analyze them all at once in a selected contrast.

### Chat Tab
Free-form conversation with your data:
- Ask complex questions spanning multiple genes or contrasts
- Request comparisons between conditions
- Get explanations of statistical results
- Ask for biological interpretations

Example queries:
- "Compare the expression of TP53, MYC, and GAPDH across all my contrasts"
- "What are the most significant changes in the stress vs control comparison?"
- "Explain what a log fold change of 2.5 means biologically"

## 🛠 Technical Architecture

### MCP Server (`uorca_data_server.py`)
The core data access layer that provides four main tools:

1. **`list_contrasts()`** - Returns all available experimental contrasts
2. **`list_genes(limit)`** - Returns available gene symbols (with optional limit)
3. **`get_lfc(contrast, gene)`** - Gets log fold change for specific gene-contrast pair
4. **`get_analysis_info()`** - Returns metadata about analyses and datasets

### AI Agent Factory (`ai_agent_factory.py`)
Creates and configures the AI agent with:
- OpenAI GPT-4o-mini model (optimized for cost and performance)
- Low temperature (0.1) for consistent, factual responses
- Specialized system prompt for bioinformatics context
- Connection to the UORCA data MCP server

### Landing Page (`ai_landing_page.py`)
Streamlit interface components:
- Tab-based navigation
- Real-time query execution with status indicators
- Error handling and user guidance
- Session state management for chat history

## 🔍 Data Access

### Supported Data Types
- **DEG Files**: Differential expression results from RNA-seq analysis
- **CPM Files**: Counts per million normalized expression data
- **Analysis Metadata**: Study information, sample groups, experimental design

### Data Format Requirements
The AI assistant expects standard UORCA output structure:
```
results_directory/
├── analysis_id/
│   ├── RNAseqAnalysis/
│   │   ├── contrast_name/
│   │   │   └── DEG.csv
│   │   └── CPM.csv
│   └── metadata/
│       └── analysis_info.json
```

### Gene and Contrast Identifiers
- **Genes**: Standard gene symbols (e.g., TP53, MYC, GAPDH)
- **Contrasts**: Format "analysis_id:contrast_name" (e.g., "GSE12345:treatment_vs_control")

## 🐛 Troubleshooting

### Common Issues

**"AI Assistant not available"**
- Check that all dependencies are installed: `uv sync`
- Verify imports are working: `uv run python test_ai_setup.py`

**"OpenAI API Key Required"**
- Set the environment variable: `export OPENAI_API_KEY=your_key`
- Or add to `.env` file in project root
- Verify with: `echo $OPENAI_API_KEY`

**"Failed to initialize data access"**
- Check that results directory contains valid UORCA analyses
- Verify file permissions and directory structure
- Run data access test: see test script output

**"Query failed" or timeouts**
- Check internet connection for OpenAI API access
- Verify API key is valid and has credits
- Try simpler queries first to test connectivity

### Debug Mode
Enable detailed logging by setting environment variable:
```bash
export UORCA_AI_DEBUG=1
```

### Performance Tips
- First query may be slower due to initialization
- Complex queries with many genes may take longer
- Use gene limits for large datasets to improve response times

## 🔒 Security Considerations

- **API Key**: Never commit your OpenAI API key to version control
- **Data Privacy**: UORCA data is processed locally; only metadata queries are sent to OpenAI
- **Rate Limits**: The assistant respects OpenAI rate limits automatically

## 🆘 Support

If you encounter issues:

1. **Run the test script**: `uv run python test_ai_setup.py`
2. **Check logs**: Look for error messages in the Streamlit console
3. **Verify setup**: Ensure all requirements are met (API key, data directory, dependencies)
4. **Start simple**: Try basic queries before complex ones

## 📈 Future Enhancements

Planned features:
- Support for additional analysis types (pathway analysis, GO enrichment)
- Export capabilities for AI-generated insights
- Integration with visualization tools
- Custom query templates
- Multi-modal analysis support

---

**Note**: The AI assistant is designed to complement, not replace, detailed bioinformatics analysis. Always validate important findings using traditional statistical methods and domain expertise.
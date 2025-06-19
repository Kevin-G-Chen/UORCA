UORCA/repository_navigation.md
# UORCA Repository Navigation Guide for Agents

This document provides an overview of the repository structure of the UORCA project to assist AI agents in efficiently navigating the repository for their tasks.

## Overall notes

The repository uses `uv` for package management. This means that you will typically conduct tests and run scripts using `uv run <script_name>` or `uv run <command>`, rather than directly executing Python files. This is important to remember when navigating the repository and running code.

---

## High-Level Structure

- **main_workflow/**: Core directory containing all major code and scripts for the UORCA RNA-seq analysis pipeline.
- **sample_inputs/**: Contains example input files (CSV) for testing and running batch analyses.
- **scratch/**: Miscellaneous information not needed for code editing.
- **data/**: Raw data files not needed for code editing.
- **archive/**: Contains old code that is not relevant for current development or editing.

---

## main_workflow/ Directory Structure and Key Components

### 1. **Agents**

Contains specialized agents handling different parts of the workflow:

- **agents/extraction.py**
  Handles data extraction from GEO/SRA, including metadata fetching and FASTQ downloading. It updates context with necessary files for downstream analysis.

- **agents/metadata.py**
  Processes and cleans sample metadata, identifies biologically relevant columns, merges columns into a grouping variable, and generates contrasts for differential expression.

- **agents/analysis.py**
  Performs RNA-seq data analysis. Runs Kallisto quantification, prepares edgeR analysis sample mapping, and conducts differential expression analysis with edgeR/limma via R scripts.

- **agents/reporting.py** (legacy - reporting functionality has been moved to main_workflow/reporting/ with comprehensive modular architecture)

### 2. **Additional scripts**

- **additional_scripts/RNAseq.R**
  R script that performs the main edgeR/limma differential expression workflow, visualizations, and result export.

- **additional_scripts/DatasetIdentification.py**
  Identifies relevant GEO datasets based on a research query using AI-powered term extraction, GEO search, ESC, and clustering.

- **additional_scripts/ResultsIntegration.py**
  Integrates results from multiple RNA-seq analyses, identifies important genes, and creates interactive heatmaps and expression plots.

### 3. **Run helpers**

Scripts to help submit and manage SLURM batch jobs:

- **run_helpers/submit_datasets.py**
  Handles batch submission of multiple datasets with resource (storage) management and monitors job status.

- **run_helpers/run_single_dataset.sbatch.j2** and **run_helpers/run_dataset_array.sbatch.j2**
  Jinja2 job script templates for single and array dataset submission.

### 4. **Reporting**

The reporting directory contains a comprehensive modular Streamlit application (UORCA Explorer) and supporting tools for interactive analysis and visualization:

#### Core Integration and Analysis
  - **ResultsIntegration.py**: Main integration script for multi-dataset results with advanced features:
    - Cross-dataset gene prioritization with intelligent ranking
    - Interactive clustered heatmaps with hierarchical clustering
    - Expression plots with sample grouping and statistical filtering
    - Cached important gene identification for performance
    - HTML report generation with publication-quality visualizations

  - **single_analysis_plots.py**: Utilities for creating interactive QC and differential expression plots:
    - PCA plots replacing static MDS plots with sample grouping
    - Interactive volcano plots and MA plots with hover information
    - Clustered DEG heatmaps with sample group annotations
    - All plots feature hover data and biological context

#### UORCA Explorer - Modular Streamlit Application
  - **uorca_explorer.py**: Main modular Streamlit app with advanced caching and fragment isolation

  - **streamlit_tabs/**: Modular tab architecture for maintainability:
    - **data_selection_tab.py**: Interactive tables for dataset and contrast selection
    - **heatmap_tab.py**: Clustered heatmaps with advanced filtering and significance thresholds
    - **expression_plots_tab.py**: Violin plots with pagination and sample grouping
    - **analysis_plots_tab.py**: Quality control and differential expression plots from individual datasets
    - **datasets_info_tab.py**: Dataset metadata browsing with study details and filtering
    - **contrasts_info_tab.py**: Contrast information with descriptions and DEG counts
    - **ai_assistant_tab.py**: Complete AI-powered analysis with transparency features
    - **sidebar_controls.py**: Comprehensive parameter controls with auto gene selection
    - **helpers/**: Shared utilities including caching, logging, and session state management

#### AI-Powered Analysis Infrastructure
  - **ai_agent_factory.py**: Factory for creating AI agents that connect to UORCA MCP servers:
    - Configures pydantic-ai agents with OpenAI models
    - Manages MCP server connections for tool access
    - Handles structured output validation with schemas

  - **ai_gene_schema.py**: Pydantic models for structured AI analysis output:
    - Gene selection validation (1-50 genes)
    - Filter parameter validation
    - Biological interpretation requirements

  - **mcp_server_core.py**: MCP server providing structured tool access to UORCA data:
    - **get_most_common_genes()**: Find frequently differentially expressed genes
    - **get_gene_contrast_stats()**: Get statistics for specific genes across contrasts
    - **filter_genes_by_contrast_sets()**: Find genes specific to contrast subsets
    - **summarize_contrast()**: Provide contrast summaries with top DEGs
    - **Tool call logging**: Transparent tracking of AI tool usage with parameters and results

  - **contrast_relevance.py**: AI-powered contrast relevance assessment:
    - Automated scoring of contrast relevance to research questions (0-1 scale)
    - Intelligent contrast selection with diversity optimization
    - Categorization of contrasts (primary, control, comparative, supportive)
    - Batch processing with parallel API calls and result aggregation

#### Supporting Infrastructure
  - **example_output_files/**: Example files demonstrating output structure:
    - **analysis_info.json**: Analysis metadata with checkpoints and organism information
    - **contrasts.csv**: Contrast definitions with descriptions and expressions
    - **edger_analysis_samples.csv**: Sample mapping for expression analysis
    - **GSE111143_metadata.csv**: Example GEO metadata structure

  - **prompts/**: AI prompt templates for various analysis workflows:
    - **assess_and_select_contrasts.txt**: Prompts for intelligent contrast selection
    - **assess_contrast_relevance.txt**: Prompts for contrast relevance scoring

#### Key Features and Capabilities
  - **Performance Optimization**: Advanced caching of expensive operations, fragment isolation for efficient updates
  - **AI Transparency**: Real-time display of AI tool usage with parameters and results
  - **Modular Architecture**: Separate modules for each functionality enabling independent development
  - **Interactive Visualizations**: All plots feature hover information, biological context, and export capabilities
  - **Comprehensive Analysis**: From raw results to publication-ready insights with AI-powered interpretation


### 5. **Prompts**

- **prompts/**: Contains text prompt templates used by agents for different workflow steps, including:

  - `analysis.txt`: System prompt guiding the RNA-seq analysis agent.

  - `extraction.txt`: Prompt for the data extraction agent.

  - `metadata.txt`: Prompt for the metadata processing agent.

  - `master.txt`: Orchestration agent prompt.

  - `dataset_identification/`: Prompts for dataset identification submodule.

### 6. **Shared**

- **shared/**: Core utilities and context definitions shared across agents, including:

  - `__init__.py`: Context models, enums, and checkpoint system.

  - `entrez_utils.py`: Safe Entrez API utilities with rate limiting.

  - `workflow_logging.py`: Logging decorators and setup utilities.

---

## Notes on Other Directories

- **scratch/**: Contains miscellaneous information and scripts not required for coding or editing.

- **data/**: Includes raw data files used for testing or references but not required for code editing.

- **archive/**: Contains outdated or retired code not relevant for current analysis or development.

---

## Usage Summary for Agents

- When performing **data extraction**, focus on `agents/extraction.py`.

- For **metadata processing and contrast design**, use `agents/metadata.py`.

- For **RNA-seq quantification and differential expression analysis**, utilize `agents/analysis.py` along with the R script in `additional_scripts/RNAseq.R`.

- To **generate integrated reports or explore results**, use the scripts and APIs under `main_workflow/reporting/`.

- For managing **batch runs and resource-aware scheduling**, look in `run_helpers/`.

- For **AI assistance or dataset identification**, use `additional_scripts/DatasetIdentification.py` and the AI prompt definitions.

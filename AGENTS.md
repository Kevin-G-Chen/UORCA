# UORCA Repository Navigation Guide for Agents

This document provides an overview of the repository structure of the UORCA project to assist AI agents in efficiently navigating the repository for their tasks.

## Overall notes

The repository uses `uv` for package management and has been restructured as a proper Python package. UORCA can now be installed and run via a unified CLI interface with the `uorca` command. You will typically run commands using `uv run uorca <subcommand>` or directly via `uorca <subcommand>` if installed. The package supports containerized execution via Docker and Apptainer.

---

## High-Level Structure

- **uorca/**: Core Python package containing CLI interface and batch processing modules
- **main_workflow/**: Core directory containing all major code and scripts for the UORCA RNA-seq analysis pipeline
- **sample_inputs/**: Contains example input files (CSV) for testing and running batch analyses
- **logs/**: Runtime logs and AI tool call tracking
- **scratch/**: Miscellaneous information not needed for code editing
- **data/**: Raw data files not needed for code editing

---

## Package Structure (uorca/)

The main package provides a unified CLI interface and modular batch processing:

### 1. **CLI Interface**

- **uorca/cli.py**
  Main command-line interface providing three primary subcommands:
  - `uorca identify`: Dataset identification from GEO using AI-powered relevance assessment
  - `uorca run`: Batch RNA-seq analysis pipeline (supports both SLURM and local execution)
  - `uorca explore`: Interactive results explorer web application

- **uorca/identify.py**
  Wrapper for dataset identification functionality, integrating with the main workflow

- **uorca/explore.py**
  Launcher for the Streamlit-based UORCA Explorer web application

### 2. **Batch Processing Framework**

Modular batch processing system supporting different execution backends:

- **uorca/batch/base.py**
  Abstract base class defining the batch processing interface

- **uorca/batch/slurm.py**
  SLURM cluster batch processing implementation with resource management and job monitoring

- **uorca/batch/local.py**
  Local multiprocessing batch execution with automatic resource detection

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

### 2. **Dataset Identification**

AI-powered dataset discovery and relevance assessment:

- **dataset_identification/DatasetIdentification.py**
  Identifies relevant GEO datasets based on a research query using AI-powered term extraction, GEO search, clustering, and multi-round relevance scoring. Supports configurable thresholds and batch processing.

### 3. **Additional scripts**

- **additional_scripts/RNAseq.R**
  R script that performs the main edgeR/limma differential expression workflow, visualizations, and result export.

### 4. **Run helpers**

Scripts to help submit and manage SLURM batch jobs (legacy - functionality moved to uorca.batch):

- **run_helpers/submit_datasets.py**
  Handles batch submission of multiple datasets with resource (storage) management and monitors job status.

- **run_helpers/run_single_dataset.sbatch.j2** and **run_helpers/run_dataset_array.sbatch.j2**
  Jinja2 job script templates for single and array dataset submission.

### 5. **Reporting**

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

  - **tool_relevance_analyzer.py**: Advanced tool usage analysis and optimization for AI workflows

  - **config_loader.py**: Configuration management for the reporting system

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

### 6. **Prompts**

- **prompts/**: Contains text prompt templates used by agents for different workflow steps, including:

  - `analysis.txt`: System prompt guiding the RNA-seq analysis agent.

  - `extraction.txt`: Prompt for the data extraction agent.

  - `metadata.txt`: Prompt for the metadata processing agent.

  - `master.txt`: Orchestration agent prompt.

### 7. **Shared**

- **shared/**: Core utilities and context definitions shared across agents, including:

  - `__init__.py`: Context models, enums, and checkpoint system.

  - `entrez_utils.py`: Safe Entrez API utilities with rate limiting.

  - `workflow_logging.py`: Logging decorators and setup utilities.

### 8. **Master Workflow**

- **master.py**: Main orchestration script that coordinates the entire workflow from data extraction through analysis and reporting.

---

## Supporting Infrastructure

### 1. **Container Support**

- **run_uorca_explorer.sh**: Comprehensive container launcher script supporting both Docker and Apptainer execution with automatic engine detection, environment variable handling, and SSH tunneling instructions.

- **Dockerfile**: Container definition for reproducible UORCA deployments.

### 2. **Package Configuration**

- **pyproject.toml**: Python package configuration with dependencies, entry points, and workspace definitions. Defines the `uorca` CLI command.

- **uv.lock**: Dependency lock file for reproducible environments.

### 3. **Development and Testing**

- **logs/**: Runtime logging directory containing:
  - AI tool call tracking and transparency logs
  - Dataset identification logs
  - Streamlit application logs

---

## Notes on Other Directories

- **scratch/**: Contains miscellaneous information and scripts not required for coding or editing.

- **data/**: Includes raw data files used for testing or references but not required for code editing.

---

## CLI Usage Summary

The UORCA package provides a unified command-line interface:

### Dataset Identification
```bash
uorca identify -q "cancer stem cell differentiation" -o results.csv
```

### Batch Analysis
```bash
# SLURM cluster execution
uorca run slurm --csv datasets.csv --output_dir ../UORCA_results

# Local multiprocessing execution
uorca run local --csv datasets.csv --output_dir ../UORCA_results
```

### Results Exploration
```bash
uorca explore ../UORCA_results --port 8501
```

## Usage Summary for Agents

- When performing **dataset identification**, use the CLI command `uorca identify` or work directly with `main_workflow/dataset_identification/DatasetIdentification.py`.

- For **data extraction**, focus on `agents/extraction.py`.

- For **metadata processing and contrast design**, use `agents/metadata.py`.

- For **RNA-seq quantification and differential expression analysis**, utilize `agents/analysis.py` along with the R script in `additional_scripts/RNAseq.R`.

- To **generate integrated reports or explore results**, use the scripts and APIs under `main_workflow/reporting/` or run `uorca explore`.

- For managing **batch runs**, use the CLI commands `uorca run slurm` or `uorca run local`, or work directly with the batch processing modules in `uorca/batch/`.

- For **containerized execution**, use `run_uorca_explorer.sh` or the Docker/Apptainer containers directly.

- For **AI assistance or dataset identification**, use the CLI or the integrated AI analysis features in the UORCA Explorer.

- All components support the `uv` package manager workflow and can be run via `uv run` commands or direct CLI execution after installation.

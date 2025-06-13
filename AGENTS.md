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

- **agents/reporting.py** (not mentioned explicitly, but reporting functionality is in main_workflow/reporting/ - you can likely ignore this)

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

- **reporting/**: Contains tools and interfaces for report generation and data exploration, including:

  - **ResultsIntegration.py** (see above): main integration script for multi-dataset results.

  - **single_analysis_plots.py**: Utilities for plotting QC and differential expression plots for single analyses.

  - **ai_agent_factory.py**: Factory that creates AI agents to interact with UORCA results via MCP.

  - **ai_landing_page.py**: Streamlit-based AI assistant landing page for interactive data exploration.

  - **mcp_servers/uorca_data_server.py**: MCP server providing JSON/structured data tool access to UORCA analyses.

  - **mcp_utils.py**: Helper functions for MCP server setup.

  - **reporting/uorca_explorer.py**: Streamlit app for interactive exploration of UORCA RNA-seq results.

  - **reporting/MCP_examples/**: Contains example MCP server and client demonstrating how to work with MCP.
  - **reporting/example_output_files/** Contains example files from the output of the main workflow. Refer to these files to understand the sturcture of the output, and therefore how to approach interpreting these files for the purpose of the reporting app (e.g. JSON outputs)


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

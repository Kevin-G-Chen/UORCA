# Unified -Omics Reference Corpus of Analyses (UORCA)

A fully containerized, AI-powered workflow for automated RNA-seq analysis of public datasets from the Gene Expression Omnibus (GEO).

## Objective and purpose

UORCA addresses a critical bottleneck in genomics research: the time-consuming process of manually analyzing public RNA-seq datasets. While thousands of high-quality datasets exist in public repositories like GEO, researchers often lack the time, expertise, or computational resources to systematically analyze them for their research questions.

UORCA provides:
- **Automated end-to-end analysis**: From GEO accession to publication-ready results
- **AI-driven decision making**: Intelligent metadata processing and experimental design
- **Containerized deployment**: No complex dependency management required
- **Interactive exploration**: Web-based tools for results visualization and analysis
- **Scalable processing**: SLURM-based batch processing for multiple datasets

The workflow automatically handles:
- Data extraction from GEO/SRA
- Metadata cleaning and experimental design
- RNA-seq quantification and differential expression analysis
- Statistical analysis and visualization
- Report generation with publication-quality figures

## What specific value would UORCA provide?

1. **Time savings**: Reduces weeks of manual analysis to hours of automated processing
2. **Consistency**: Standardized analysis pipeline ensures reproducible results across datasets
3. **Accessibility**: No bioinformatics expertise required - just provide a GEO accession
4. **Scalability**: Process dozens of datasets simultaneously using HPC resources
5. **Quality**: AI-powered metadata processing handles complex experimental designs
6. **Integration**: Tools for meta-analysis across multiple datasets
7. **Exploration**: Interactive web interface for deep-dive analysis of results

## Quickstart

Get started with UORCA in just a few steps:

```bash
# Clone the repository
git clone https://github.com/Kevin-G-Chen/UORCA.git
cd UORCA

# Set up the uv environment
uv venv
uv pip install -e .

# 2. Find relevant datasets for your research
uv run uorca identify -q "cancer stem cell differentiation" -o datasets.csv

# 3. Process datasets (choose your environment)
uv run uorca run slurm --csv datasets.csv --output_dir ../UORCA_results  # For HPC clusters
# OR
uv run uorca run local --csv datasets.csv --output_dir ../UORCA_results  # For local machines

# 4. Explore results interactively
uv run uorca explore ../UORCA_results
```

## End Result

For each processed dataset, UORCA generates:
- Comprehensive HTML reports with statistical summaries and visualizations
- Differential expression tables for all experimental contrasts
- Publication-quality plots (PCA, volcano plots, heatmaps, etc.)
- Processed count matrices and metadata files
- Interactive results that can be explored via the UORCA Explorer web app

For multiple datasets, additional integration tools provide:
- Meta-analysis across datasets identifying consistent patterns
- Interactive multi-dataset visualizations
- Integrated gene prioritization and pathway analysis

## Who will benefit from UORCA?

- **Biologists** seeking to leverage public data without bioinformatics expertise
- **Bioinformaticians** looking to standardize and accelerate their analysis workflows
- **Research groups** wanting to systematically analyze multiple related datasets
- **Meta-analysis researchers** integrating findings across studies
- **Educational institutions** teaching RNA-seq analysis with real datasets

## Prerequisites

UORCA runs entirely within an Apptainer/Singularity container, eliminating the need for complex dependency management. The only requirements are:

- Access to an HPC system with SLURM and Apptainer/Singularity
- The UORCA container file (`.sif`)
- OpenAI API key for AI-powered analysis (set as environment variable)

## Navigating this repository

Everything you need for running UORCA lives under the `main_workflow/` folder. The rest of the repo contains documentation, helper scripts, and containerization files.

### Container-based workflow

All UORCA analyses run within an Apptainer container that includes all necessary dependencies:
- Python environment with scientific computing libraries
- R with Bioconductor packages for RNA-seq analysis
- Command-line tools (Kallisto, SRA Toolkit, etc.)
- Web frameworks for interactive visualization

## Unified Command Line Interface

UORCA provides a unified CLI that simplifies access to all three main workflows:

```bash
# Install UORCA package (enables CLI)
uv pip install -e .

# Dataset identification - find relevant GEO datasets
uv run uorca identify -q "cancer stem cell differentiation" -o results.csv

# Dataset analysis - run batch RNA-seq pipeline
uv run uorca run slurm --csv datasets.csv --output_dir ../UORCA_results

# Results exploration - launch interactive web app
uv run uorca explore ../UORCA_results --port 8501
```

**Key advantages of the unified CLI:**
- **Shorter commands**: `uv run uorca identify` vs `uv run main_workflow/dataset_identification/DatasetIdentification.py`
- **Simplified deployment**: `uorca explore` runs directly without container dependency
- **Batch-oriented processing**: `uorca run` focuses on efficient multi-dataset processing
- **Consistent interface**: All commands follow the same pattern with standardized help
- **Easy discovery**: `uv run uorca --help` and `uv run uorca COMMAND --help` for documentation

### Detailed CLI Usage

**Dataset Identification:**
```bash
# Basic usage
uv run uorca identify -q "neuroblastoma tumor vs normal" -o my_results

# Advanced options
uv run uorca identify -q "stem cell differentiation" \
    --threshold 7.5 \
    --max-datasets 1000 \
    --scoring-rounds 5 \
    --verbose
```

**Dataset Analysis:**
```bash
# SLURM batch processing
uv run uorca run slurm --csv datasets.csv --output_dir ../UORCA_results --cleanup

# Local parallel processing
uv run uorca run local --csv datasets.csv \
    --output_dir /path/to/results \
    --resource_dir ./data/kallisto_indices \
    --max_workers 4 --cleanup
```

**Results Exploration:**
```bash
# Explore specific results directory
uv run uorca explore /path/to/results

# Custom settings with different port
uv run uorca explore /path/to/results --port 8502 --headless

# For remote access via SSH tunnel
uv run uorca explore ../UORCA_results --port 8501
```

#### Single dataset analysis

```bash
# Using the container directly
./scratch/sbatch_script/Run_MasterAgent.sh

# Or submit via SLURM with custom parameters
sbatch --output=logs/my_analysis.out \
       --error=logs/my_analysis.err \
       Run_MasterAgent.sh
```

#### Multiple dataset analysis

UORCA provides efficient batch processing for multiple datasets:
```bash
# SLURM batch processing
uv run uorca run slurm --csv datasets.csv --output_dir ../UORCA_results --cleanup

# Local parallel processing
uv run uorca run local --csv datasets.csv --output_dir ../UORCA_results --max_workers 4 --cleanup
```

For legacy batch processing with advanced resource management:
```bash
# Process multiple datasets with storage-aware scheduling
uv run main_workflow/run_helpers/submit_datasets.py \
    --csv_file datasets.csv \
    --output_dir ../UORCA_results \
    --max_parallel 10 \
    --max_storage_gb 500 \
    --cleanup
```

### Agentic workflow architecture

UORCA employs a multi-agent AI system where specialized agents handle different aspects of the analysis:

1. **Master agent (`master.py`)**
   - Orchestrates the entire workflow
   - Handles CLI inputs and coordinates between agents
   - Manages checkpointing and error recovery

2. **Extraction agent (`agents/extraction.py`)**
   - Downloads GEO Series metadata and sample information
   - Fetches FASTQ files from SRA using optimized parallel downloading
   - Automatically determines organism from taxonomic metadata

3. **Metadata agent (`agents/metadata.py`)**
   - Cleans and standardizes sample metadata
   - Intelligently merges biological replicates and conditions
   - Designs experimental contrasts for differential expression

4. **Analysis agent (`agents/analysis.py`)**
   - Runs Kallisto quantification with automatic parallelization
   - Performs edgeR/limma differential expression analysis
   - Generates statistical summaries and quality control metrics

5. **Reporting agent (`agents/reporting.py`)**
   - Creates publication-quality visualizations
   - Builds comprehensive HTML reports using Sphinx
   - Organizes results for downstream exploration

Each agent has access to specialized tools and can make autonomous decisions while maintaining full traceability of the analysis process.

## Additional scripts and tools

### Dataset identification

Automatically discover relevant GEO datasets for your research question:

```bash
# Using unified CLI (recommended)
uv run uorca identify -q "cancer stem cell differentiation" -o relevant_datasets.csv

# Or using direct script execution
uv run main_workflow/dataset_identification/DatasetIdentification.py \
    -q "cancer stem cell differentiation" \
    -t 7.0 \
    -o relevant_datasets.csv
```

**Features:**
- AI-powered search term extraction from research queries
- Automated relevance scoring of GEO datasets
- SRA metadata validation for RNA-seq experiments
- Direct output formatting for batch processing

### Multi-dataset processing
#### Multiple dataset analysis

Process datasets using the batch-oriented CLI:

```bash
# SLURM batch processing (for HPC clusters)
uv run uorca run slurm --csv datasets.csv --output_dir ../UORCA_results --cleanup

# Local parallel processing (for workstations)
uv run uorca run local --csv datasets.csv --output_dir ../UORCA_results --max_workers 4 --cleanup

# Legacy batch processing with storage-aware scheduling
uv run main_workflow/run_helpers/submit_datasets.py \
    --csv_file multi_dataset_input.csv \
    --output_dir ../UORCA_results \
    --resource_dir ./data/kallisto_indices \
    --max_parallel 10 \
    --max_storage_gb 500 \
    --cleanup
```

**Key features:**
- Intelligent job scheduling based on dataset size
- Storage-aware resource management to prevent disk overflow
- Automatic cleanup of intermediate files
- Comprehensive progress tracking and error reporting
- Results ready for immediate exploration in UORCA Explorer

**CSV file format:**
```csv
Accession,DatasetSizeBytes
GSE123456,15000000000
GSE789012,8500000000
```

### Results integration

Integrate findings across multiple datasets:

```bash
# Create integrated analysis across all results
uv run main_workflow/additional_scripts/ResultsIntegration.py \
    --results_dir ../UORCA_results \
    --output_dir ../UORCA_results/integrated_results \
    --pvalue_threshold 0.05 \
    --lfc_threshold 1.0

# Then explore integrated results using the CLI
uv run uorca explore ../UORCA_results/integrated_results
```

**Features:**
- Cross-dataset gene prioritization with intelligent ranking
- Interactive multi-dataset heatmaps with clustering
- Meta-analysis statistical summaries
- Publication-ready integrated visualizations
- Seamless integration with UORCA Explorer for interactive analysis

## UORCA Explorer: Interactive Results Exploration

UORCA Explorer is a containerized Streamlit web application for interactive exploration of analysis results, featuring a modern modular architecture for enhanced performance and maintainability.

### Container-based usage

**Running the Explorer:**

```bash
# Using the CLI (recommended - runs directly without container)
uv run uorca explore ../UORCA_results --port 8501

# Or using containerized script (for container environments)
./run_uorca_explorer.sh [results_directory] [port]

# Example with custom settings
./run_uorca_explorer.sh ../UORCA_results 8501
```

**Remote access:**

1. Start the Explorer on your HPC system:
   ```bash
   # Using the CLI (direct Python execution)
   uv run uorca explore ../UORCA_results --port 8501

   # Or using container script (containerized execution)
   ./run_uorca_explorer.sh ../UORCA_results 8501
   ```

2. Create SSH tunnel from your laptop:
   ```bash
   ssh -L 8000:127.0.0.1:8501 your_username@hpc_server
   ```

3. Open `http://127.0.0.1:8000` in your browser

### Key Features

**Interactive Analysis:**
- **Smart Gene Selection**: Auto-identified DEGs with intelligent ranking or custom gene lists with validation
- **Dataset & Contrast Selection**: Interactive tables for selecting subsets of datasets and contrasts for focused analysis
- **Dynamic Visualization**: Real-time heatmaps, expression plots, and statistical summaries with advanced caching
- **Export Capabilities**: Download results as CSV or interactive HTML reports

**Advanced Exploration:**
- **Multi-dataset Integration**: Compare patterns across multiple studies with clustered heatmaps
- **Statistical Filtering**: Adjustable significance thresholds and fold-change cutoffs with separate heatmap filters
- **Quality Control**: Dataset-specific diagnostic plots including PCA, volcano plots, and MA plots
- **Analysis Plots**: Interactive alternatives to static R plots with hover information and sample grouping

**AI-Powered Analysis with Transparency:**
- **One-Click Complete Analysis**: Streamlined workflow combining contrast relevance and gene analysis
- **Contrast Relevance Assessment**: AI scoring and intelligent selection of contrasts relevant to research questions
- **Tool Call Transparency**: Real-time display of AI tool usage with parameters and results for full transparency
- **Automated Gene Discovery**: AI-driven identification of key differential expression patterns across selected contrasts
- **Intelligent Threshold Selection**: AI automatically chooses appropriate statistical thresholds based on data characteristics
- **Pattern Recognition**: AI identifies both shared and context-specific gene expression signatures
- **Biological Interpretation**: AI provides structured reasoning and biological context for findings

**Performance & Architecture:**
- **Modular Design**: Separate tab modules for better maintainability and development workflow
- **Advanced Caching**: Expensive operations like gene identification cached across sessions
- **Fragment Isolation**: Efficient tab switching without recomputation using Streamlit fragments
- **Automatic Pagination**: Large gene sets handled efficiently with 30 genes per page
- **Responsive Design**: Professional styling optimized for analysis workflows

### Available Tabs

**Data Selection & Exploration:**
- **☑️ Select Data & Contrasts**: Interactive tables for choosing datasets and contrasts with DEG counts
- **🌡️ Explore DEG Heatmap**: Clustered heatmaps of log2 fold changes with hover information and filtering
- **📈 Plot Gene Expression**: Violin plots showing expression distributions across sample groups
- **🧑‍🔬 Analyze Experiments**: Quality control and differential expression plots from individual datasets

**Information & AI Analysis:**
- **📋 View Dataset Info**: Browse and filter dataset metadata with study details and descriptions
- **🔍 View Contrast Info**: Browse contrast details with descriptions and DEG counts
- **🤖 AI Assistant**: Complete AI-powered analysis with contrast relevance assessment and gene discovery

### Troubleshooting

**Container Issues:**
- Ensure the container file exists at the specified path
- Verify Apptainer/Singularity module is loaded
- Check that results directory contains properly formatted analysis outputs

**Access Issues:**
- Confirm SSH tunnel is active for remote access
- Verify firewall settings allow the specified port
- Check that no other services are using the same port

**Performance Tips:**
- Use auto-selected DEGs for optimal performance with cached gene identification
- Enable cleanup mode for long-running batch jobs
- Monitor storage usage when processing large datasets
- Use the dedicated "Select Data & Contrasts" tab for efficient dataset browsing
- AI analysis benefits from OpenAI API key configuration for full functionality

## File organization

UORCA generates a structured output hierarchy:

```
UORCA_results/
├── GSE123456/                    # Individual dataset results
│   ├── metadata/                 # Sample information and experimental design
│   ├── quantification/           # Kallisto abundance files
│   ├── analysis/                 # Differential expression results
│   ├── figures/                  # Publication-quality plots
│   └── report/                   # HTML report and summary
├── integrated_results/           # Multi-dataset integration
├── logs/                         # Analysis logs and progress tracking
└── job_status/                   # SLURM job management files
```

## Getting started

1. **Obtain the container**: Contact the UORCA development team for access to the latest container image

2. **Set up your environment**:
   ```bash
   export OPENAI_API_KEY="your_api_key_here"  # Required for AI Assistant functionality
   ```

3. **Run a test analysis**:
   ```bash
   # Modify the script to point to your container location
   ./scratch/sbatch_script/Run_MasterAgent.sh
   ```

4. **Explore results interactively**:
   ```bash
   ./run_uorca_explorer.sh ../UORCA_results
   ```

5. **Use the AI Assistant** (if OpenAI API key is configured):
   - Navigate to the "🤖 AI Assistant" tab
   - Enter your research question
   - Click "🚀 Run Complete AI Analysis" for automated contrast relevance assessment and gene discovery
   - View transparent AI tool usage and biological interpretations

For questions, issues, or feature requests, please contact the UORCA development team or submit an issue through the project repository.

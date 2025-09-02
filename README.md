# Unified -Omics Reference Corpus of Analyses (UORCA)

A fully containerized, AI-powered workflow for automated RNA-seq analysis of public datasets from the Gene Expression Omnibus (GEO).

## Preparation

The following steps will be needed to set up UORCA:

```bash
# 1. Clone the repository
git clone https://github.com/Kevin-G-Chen/UORCA.git
cd UORCA

# 2. Set up the uv environment
uv venv
uv pip install -e .

# 3. Download Kallisto indices (required for RNA-seq quantification)
./download_kallisto_indices.sh        # Downloads human, mouse, dog, monkey, zebrafish
# OR download specific species:
./download_kallisto_indices.sh human  # Download only human index

# 4. (Optional) Pull the container for containerized execution
# For Singularity/Apptainer (recommended for HPC/SLURM):
singularity pull uorca_0.1.0.sif docker://kevingchen/uorca:0.1.0
# OR Apptainer (newer name, same command):
apptainer pull uorca_0.1.0.sif docker://kevingchen/uorca:0.1.0

# For Docker (local development only):
docker pull kevingchen/uorca:0.1.0
```

## 🔑 Environment Setup

UORCA requires API credentials for accessing biological databases and AI services.

### Required Variables
- **`ENTREZ_EMAIL`** - Your email address (required by NCBI guidelines)
- **`OPENAI_API_KEY`** - OpenAI API key (required for AI-powered dataset identification)

### Optional (Recommended)
- **`ENTREZ_API_KEY`** - NCBI API key (enables >3x faster processing - free from [NCBI](https://www.ncbi.nlm.nih.gov/account/settings/))

### Setup Instructions
1. Copy the template: `cp .env.example .env`
2. Edit `.env` with your actual values:
```bash
ENTREZ_EMAIL=your.email@institution.edu        # Required: Any valid email
OPENAI_API_KEY=sk-proj-your-key-here          # Required: Get from https://platform.openai.com/api-keys
ENTREZ_API_KEY=your-ncbi-key                  # Optional: For faster processing
```

### Troubleshooting
- **Missing variables**: UORCA will show clear error messages with setup instructions
- **Invalid OpenAI key**: Check your key format (should start with `sk-proj-` or `sk-`)
- **Rate limiting**: Add `ENTREZ_API_KEY` for faster NCBI queries (3→10 requests/second)

## 🚀 The UORCA Workflow: Identify → Run → Explore

UORCA follows a streamlined three-step workflow for comprehensive RNA-seq analysis:

### Step 1: Identify - Find relevant datasets
Discover GEO datasets that match your research question using AI-powered search:

```bash
uv run uorca identify -q "cancer stem cell differentiation" -o identification_results

# Options:
# -q: Your research question (required)
# -o: Output directory name (default: adds timestamp)
# -m: Max results to fetch for each search term (default: 500)
# -r: Number of ranking iterations (default: 3)
# -t: Relevance threshold 1-10 (default: 7.0)
# --model: GPT model to use (default: gpt-5-mini)
```

**Output**: Directory containing `selected_datasets.csv` with relevant GEO accessions and metadata about your search.

### Step 2: Run - Process datasets through the pipeline
Execute the complete RNA-seq analysis pipeline on identified datasets:

```bash
# For HPC clusters with SLURM
uv run uorca run slurm --input identification_results/ --output_dir ../UORCA_results

# For local machines
uv run uorca run local --input identification_results/ --output_dir ../UORCA_results --max_workers 4

# Options:
# --input: Directory from identify step OR CSV file
# --output_dir: Where to save results
# --cleanup: Remove intermediate files after processing
# --max_workers: Number of parallel jobs (local only)
```
**Note**: Inputting the directory rather than the CSV file is recommended, as this will ensure the research question is considered in the automated analyses. See `sample_inputs` for an example - feel free to test using this directory!

**Output**: Complete analysis results including differential expression and visualisations.

### Step 3: Explore - Interact with your results
Launch the interactive web interface to explore and analyze results:

```bash
uv run uorca explore ../UORCA_results

# Options:
# --port: Web server port (default: 8501)
# --headless: Run without opening browser

# For remote access via SSH:
# 1. On HPC: uv run uorca explore ../UORCA_results --port 8501
# 2. On laptop: ssh -L 8000:127.0.0.1:8501 username@hpc_server
# 3. Open browser: http://127.0.0.1:8000
```

**Features**: Interactive heatmaps, gene expression plots, and cross-dataset integration.

## What UORCA Does

UORCA automates the entire RNA-seq analysis workflow:

1. **Dataset Discovery**: AI-powered identification of relevant GEO datasets based on your research question
2. **Data Processing**: Automated download, quality control, and RNA-seq quantification using Kallisto
3. **Statistical Analysis**: Differential expression analysis with automatic experimental design
4. **Interactive Exploration**: Web-based interface for visualizing and analyzing results across multiple datasets
5. **AI-Powered Insights**: Intelligent analysis assistant for biological interpretation

## Prerequisites

- **Python 3.10+** with `uv` package manager
- **For HPC/SLURM**: 
  - SLURM job scheduler
  - Apptainer or Singularity (preferred over Docker for security)
  - Note: Singularity/Apptainer runs containers without root privileges, making it ideal for multi-user HPC environments
- **For local**: Sufficient storage and compute resources
- **API Keys**: OpenAI and NCBI credentials (see Environment Setup)

### Container Runtime Notes

**Why Singularity/Apptainer for HPC?**
- **Security**: Runs without root privileges (unlike Docker)
- **SLURM Integration**: Native support via `--container` flag
- **Performance**: Bare-metal performance with minimal overhead
- **Portability**: SIF files are single files, easy to share across nodes

The command `singularity pull` automatically converts Docker images to Singularity's SIF format, allowing you to use Docker Hub images while maintaining HPC security requirements.

## Output Structure

UORCA generates organized results for each dataset:

```
UORCA_results/
├── GSE123456/                           # Individual dataset results
│   ├── metadata/                        # Sample and experimental information
│   │   ├── GSE123456_metadata.csv      # Sample metadata
│   │   ├── analysis_info.json          # Analysis parameters
│   │   ├── contrasts.csv               # Experimental contrasts
│   │   └── edger_analysis_samples.csv  # Samples used in analysis
│   ├── RNAseqAnalysis/                  # Differential expression results
│   │   ├── CPM.csv                     # Counts per million
│   │   ├── DGE_norm.RDS                # Normalized expression object
│   │   ├── MDS.png                     # Multidimensional scaling plot
│   │   ├── filtering_density.png       # Expression filtering diagnostics
│   │   ├── normalization_boxplots.png  # Normalization quality control
│   │   ├── sa_plot.png                 # Sample relationship plot
│   │   ├── voom_mean_variance.png     # Mean-variance trend
│   │   └── Contrast1_vs_Contrast2/    # Per-contrast results
│   │       ├── DEG.csv                # Differentially expressed genes
│   │       ├── volcano_plot.png       # Volcano plot
│   │       ├── ma_plot.png            # MA plot
│   │       └── heatmap_top50.png     # Top 50 DEGs heatmap
│   └── logs/                           # Processing logs
│       ├── analysis_tool_logs_*.json  # Tool execution logs
│       └── *.log                       # Timestamped process logs
├── GSE789012/                          # Another dataset...
├── job_status/                         # SLURM job tracking
│   └── GSE*_status.json               # Individual job status files
└── logs/                               # Batch processing logs
    └── run_GSE*.out/.err              # SLURM output/error logs
```

## Getting Help

- **Command help**: `uv run uorca --help` or `uv run uorca COMMAND --help`
- **Issues**: Submit via [GitHub Issues](https://github.com/Kevin-G-Chen/UORCA/issues)

## Citation

If you use UORCA in your research, please cite our work (manuscript in preparation).

## License

This project is under development. License information will be added upon public release.

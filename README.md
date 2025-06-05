# Unified -Omics Reference Corpus of Analyses (UORCA)

![image](https://github.com/user-attachments/assets/1c9f9551-0ea2-4c4a-a41c-6348e9e3ebc1)

The UORCA is an agentic workflow that performs and collates analyses on multiple datasets relevant to a given research question.

## Objective and purpose

Over the past years, many biological and clinical datasets have been produced and made publicly available. These datasets span many different tissues, diseases, and are derived from different technologies. As such, identifying and analysing relevant datasets would prove valuable in providing information into any biomedical research question we may have.

However, doing so in reality is much more complicated. Some of the specific complications are as follows:
- analysing datasets traditionally requires manual efforts, such as:
  - identifying relevant datasets (does this dataset actually help with my research question?)
  - parsing metadata (which columns contain scientifically meaningful variables?)
  - designing accurate code (how do I set up a design matrix to capture the nature of different samples?)
- it is difficult to know how valuable any findings are:
  - are these DEGs specific to this disease subtype, or are they known to play critical roles in other contexts too?
  - how are these transcriptomics data related to known proteomics or genomics findings?
  - if we rely on literature review, can we be sure that authors have relegated possibly crucial information to supplementary data?

Combined, the analysis of these public datasets at present represents a manual and time-consuming process. This limits both the rate at which we can glean information out of public repositories, but also the amount of knowledge that can be gained.

To circumvent these challenges, we would need a system that can:
- automatically identify datasets
- automate analysis of these datasets (including data extraction and results generation)
- place these results in the context of other available information, e.g. from other public datasets

We can leverage the power of LLMs to achieve the above. LLMs have been used in, for example, literature review (evidencing  LLMs' ability to collate information), as well as in assisting analysis of datasets (e.g. [Bioinformatics Agent](https://www.biorxiv.org/content/10.1101/2024.05.22.595240v2) or [scChat](https://www.biorxiv.org/content/10.1101/2024.10.01.616063v2.full)). The existence of such tools indicates that similar or identical workflows are not valuable to the scientific community.

__Instead, the main area where we can differentiate, gain value, and provide meaningful scientific contributions will be integrating the findings across datasets. Doing so will allow us to derive richer conclusions from datasets.__

## What specific value would UORCA provide?

UORCA will generate "insights" by collating information from different datasets. These insights might include:
- identifying genes which are differentially expressed in disease states, including those which affect disease progression/severity
- more generally - identifying common biological trends between different datasets on same/similar diseases. Alternatively - identifying biological factors underpinning different contexts
- how genetic, transcriptomic, and proteomic data relate to eachother (e.g. is the transcriptomic information correlated to the proteomic information)
- what this might look like in practice:
  - in a specific contrast, what genes are uniquely differentially expressed?
  - what genes are differentially expressed across all/most (similar) contrasts?

These findings are nice, and would generate testable hypotheses (e.g. if we induce this mutation, findings from UORCA predict that we should induce a disease state similar to X, but distinct to Y). However, I would need to
- convince researchers that any generated insights/testable hypotheses are actually worth pursuing further
- produce findings that can be more explicitly relevant for the clinic (e.g. are these DEGs meant to therefore be drug targets, or have I included DEGs for fun?)

To convince researchers that UORCA can produce valuable findings, I need to demonstrate that any "insights" reflect underlying biology, and that I can capture novel biological ideas. For this purpose, ideas of experiments I have include (dependent on project progression):
- identifying novel clinical datasets, and demonstrating that results from UORCA can be used for diagnostic purposes (e.g. that a mutation can be predictive of disease progression). I would hope this could translate into better diagnoses for diseases.
- performing docking simulations to demonstrate the possibility of repurposing a drug for treatment (e.g. if I find that gene activity for some disease is comparable to another disease, can I support this with drug-target binding predictions?)

If I can successfully establish this, then I believe there is a chance that UORCA can be used to foster collaborations to incorporate wet lab validation.

## Envisaged end result

A web server where the user
1. Inputs a research question (regardless of whether this is a web-server or not, this should be the entry point)
2. OPTIONALLY includes their own data. The intent here would be to assist with "how does my data fit in the broader picture?"

The user should then see:
- key biological findings (e.g. important genes, proteins, transcripts, how they contribute to disease state)
- LLM inputs/outputs, e.g. datasets analysed, chosen parameters, raw results, generated code

## Would will benefit from UORCA?

1. Researchers
- Reserachers will see the current "state of play," and also be able to see where their works fits in the broader picture

2. Clinicians
- This is dependent on validating findings from UORCA - i.e. the method I use to validate generated "insights".
- If UORCA can be shown to produce meaningful results, the "testable hypotheses" generated from UORCA can foster collaborations with wet lab researchers, potentially paving the way to influencing diagnostic/treatment protocols through experimental validation of any findings
- These findings can then improve patient care

3. Patients
- With better diagnostic and treatment options, we will get better health outcomes for disease patients.
- This is dependent on clinical adoption of any findings (i.e. benefits for clinicians.)

## Internal collaboration opportunities
- SampleExplorer: this would facilitate identification of relevant samples and datasets.
- GeneInsight: I will get DEGs as part of the analysis of individual RNAseq datasets. GeneInsight could be used generate insights

## External collaboration opportunities

Given that this workflow should be agnostic to disease/gene etc., I anticipate that any researcher that would benefit from pulling from public data should see benefit here. Therefore, there are collaboration opportunities to help researchers identify datasets that are relevant to their research question, and then to analyse these datasets.

Dependencies and Installation

### Python Environment

UORCA uses `uv` for Python package management:

```bash
# Install uv
curl -LsSf https://astral.sh/uv/install.sh | sh

# Create and activate a new environment
uv venv
source .venv/bin/activate

# Install required Python packages
uv pip install openai pandas numpy pydantic pydantic-ai matplotlib tqdm scipy scikit-learn hdbscan plotly biopython unidecode python-dotenv GEOparse
```

### R Packages

R 4.1.2 or newer is recommended. Install the required R packages:

```R
# Install required packages
install.packages(c("gplots"))

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("edgeR", "tximport", "limma", "ComplexHeatmap"))
```

### Command-line Tools

#### NCBI Entrez Direct E-utilities

```bash
# Install Entrez Direct
sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
echo 'export PATH=${PATH}:$HOME/edirect' >> ~/.bashrc
source ~/.bashrc
```

#### NCBI SRA Toolkit

```bash
# Download and extract SRA Toolkit
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar -xzf sratoolkit.current-ubuntu64.tar.gz
echo 'export PATH=${PATH}:/path/to/sratoolkit.x.x.x-ubuntu64/bin' >> ~/.bashrc
source ~/.bashrc

# Configure SRA Toolkit
vdb-config --interactive
```

#### Kallisto

```bash
# Install Kallisto via conda (recommended)
conda install -c bioconda kallisto

# Or download pre-built binary
wget https://github.com/pachterlab/kallisto/releases/download/v0.46.1/kallisto_linux-v0.46.1.tar.gz
tar -xzf kallisto_linux-v0.46.1.tar.gz
echo 'export PATH=${PATH}:/path/to/kallisto' >> ~/.bashrc
source ~/.bashrc
```

#### Pigz (Parallel gzip)

```bash
# On Ubuntu/Debian
apt-get install pigz

# On CentOS/RedHat
yum install pigz

# On macOS with Homebrew
brew install pigz
```
## Navigating this repository

Everything you need for running UORCA lives under the `main_workflow/` folder.  The rest of the repo (this README, CI files, docs, etc.) is mainly for project bookkeeping— day-to-day work will be in `main_workflow/`.


### Agentic workflow

1. **Master agent (`master.py`)**
   Reads your CLI inputs (accession, output folder, resources) and then in turn invokes:
   - `extract()` → Extraction agent (automatically determines organism from GEO taxonomic metadata)
   - `analyse()` → Analysis (and metadata) agents

2. **Extraction agent (`agents/extraction.py`)**
   Tools:
   - `fetch_geo_metadata(accession)`
     • Downloads GEO Series metadata, builds a GSM ↔ SRX ↔ SRR table
   - `download_fastqs(threads, max_spots)`
     • Uses SRA-Toolkit (`prefetch` + `fasterq-dump`) to get compressed FASTQ files

3. **Metadata agent (`agents/metadata.py`)**
   Tools:
   - `process_metadata()`
     • Cleans column names/values, drops uninformative columns
   - `merge_analysis_columns(columns…)`
     • Merges one or more biological columns into a single grouping factor
   - `extract_unique_values()`
     • Lists the distinct groups for downstream contrast design

4. **Analysis agent (`agents/analysis.py`)**
   Tools:
   - `list_files(directory, pattern)`
     • Utility for finding files (indices, FASTQs, txts)
   - `run_kallisto_quantification(kallisto_index)`
     • Runs Kallisto on paired FASTQs in parallel
   - `prepare_edgeR_analysis()`
     • Matches abundance files to metadata, builds sample mapping CSV
   - `process_metadata_with_agent()`
     • Proxies into the Metadata agent to pick/merge columns & build contrasts
   - `run_edger_limma_analysis(tx2gene_path)`
     • Calls the R script (`RNAseq.R`), runs edgeR/limma, generates plots and DEG tables

5. **Reporting agent (`agents/reporting.py`)**
   Tools:
   - `identify_png_files()`
     • Finds all `.png` plots in your output
   - `generate_rst_from_pngs()`
     • Copies images into a local `images/` folder and writes an RST with `.. figure::` directives
   - `build_report()`
     • Initializes a Sphinx project, injects your RST/images, runs `sphinx-build` → HTML

---

By drilling down into **`main_workflow/`**, you’ll see the full LLM-driven pipeline.  All the business logic—data extraction, metadata wrangling, differential-expression analysis, and final reporting—is orchestrated there by a small team of specialised “agents,” each with its own clear responsibility.

## Additional scripts and tools

### Dataset identification

The `DatasetIdentification.py` script automates the process of finding relevant GEO datasets for your biological research query and fetching their associated SRA metadata.

#### Features:
- Uses OpenAI to extract relevant search terms from your research query
- Searches GEO database using NCBI Entrez API
- Scores datasets for relevance to your research question
- Fetches SRA metadata for relevant datasets
- Flags valid RNA-seq runs based on library properties

#### Usage:
```bash
python main_workflow/additional_scripts/DatasetIdentification.py \
    -q "Your research query here" \
    -n 3 \                  # Number of relevance scoring repeats
    -t 7.0 \                # Minimum relevance score threshold
    -o output_file.csv \    # Output file path
    --generate-multi-csv    # Generate a CSV formatted for multi-dataset analysis
```

#### Key parameters:
- `-q, --query`: Your biological research question (required)
- `-n, --num-queries`: Number of times to repeat relevance scoring (default: 3)
- `-t, --relevance-threshold`: Minimum score to fetch SRA metadata (default: 7.0)
- `-o, --output`: Path for output CSV (default: final_combined.csv)
- `--generate-multi-csv`: Create an additional CSV file for multi-dataset analysis

### Multi-dataset processing

The `submit_datasets.py` script allows you to process multiple datasets in parallel using SLURM.

#### Usage:
```bash
uv run main_workflow/additional_scripts/submit_datasets.py \
    --csv_file multi_dataset_input.csv \
    --output_dir ../UORCA_results/output_directory \
    --resource_dir ./data/kallisto_indices \
    --max_parallel 10
```

#### Parameters:
- `--csv_file`: Path to CSV file containing datasets to process (required)
- `--output_dir`: Directory where all results will be stored (default: ../UORCA_results)
- `--resource_dir`: Directory with Kallisto indices (default: ./data/kallisto_indices/)
- `--max_parallel`: Maximum number of parallel SLURM jobs (default: 10)

#### CSV file format:
The CSV should contain at least two columns:
1. `Accession`: GEO accession numbers (e.g., GSE123456)
2. `organism`: Organism name (e.g., human, mouse) - this is maintained for multi-dataset processing but will be automatically determined for individual analyses

Additional columns can be included, but will be ignored by the script. The script will automatically generate a new CSV file with the following columns.

Example CSV content:
```
Accession,organism,RelevanceScore,Valid
GSE123456,human,8.5,Yes
GSE789012,mouse,9.2,Yes
```

You can use the `--generate-multi-csv` option in the DatasetIdentification.py script to automatically generate this CSV file.

### Results Integration

The `ResultsIntegration.py` script integrates results from multiple RNA-seq analyses, identifying important genes and generating interactive visualizations that highlight patterns across datasets.

#### Features:
- Intelligently selects important genes using two complementary strategies:
  1. Most frequently differentially expressed genes across all contrasts
  2. Contrast-specific genes with large fold changes that appear in few contrasts
- Creates interactive heatmaps showing log fold changes across contrasts
- Generates expression plots (violin/box) showing CPM values across samples
- Produces a comprehensive HTML report with interactive visualizations

#### Usage:
```bash
uv run main_workflow/additional_scripts/ResultsIntegration.py \
    --results_dir ../UORCA_results \
    --output_dir ../UORCA_results/integrated_results \
    --pvalue_threshold 0.05 \
    --lfc_threshold 1.0 \
    --top_frequent 20 \
    --top_unique 10 \
    --plot_type violin
```

#### Key parameters:
- `--results_dir`: Directory containing multiple analysis results (required)
- `--output_dir`: Directory for integrated results (default: results_dir/integrated_results_TIMESTAMP)
- `--pvalue_threshold`: P-value threshold for DEGs (default: 0.05)
- `--lfc_threshold`: Log fold change threshold for DEGs (default: 1.0)
- `--top_frequent`: Number of most frequently occurring genes to include (default: 20)
- `--top_unique`: Number of unique genes with large fold changes to include per contrast (default: 10)
- `--max_genes`: Maximum number of genes to include in visualizations (default: 100)
- `--plot_type`: Type of expression plot to generate (default: violin)
- `--gene_list`: Path to text file with specific genes to plot (optional)
- `--analysis_summary`: Generate a separate summary of all analyses and contrasts (optional)

## UORCA Explorer: Interactive Results Exploration

In addition to the HTML produced by the integration script, you can also explore results using UORCA Explorer. This is a complementary Streamlit web application designed for interactive exploration of RNA-seq analysis results generated by UORCA. It enables researchers to dive deeper into the data by providing:

- Selection of any subset of genes, datasets, or contrasts for focused analysis
- Interactive heatmaps visualizing differential expression patterns with dynamic significance filtering
- Expression plots showing gene-level counts across samples, grouped by conditions
- Dynamic filtering and search capabilities to quickly identify significant genes or contrasts
- Real-time RNA-seq analysis plots (volcano, MA plots, heatmaps) for individual datasets
- Comprehensive dataset and contrast information browsing
- Export capabilities for both interactive HTML reports and raw CSV data
- Optimized performance with fragment isolation and advanced caching

### Installation

**Prerequisites:**
- Python 3.8+
- UORCA analysis results directory

**Setup with uv (recommended):**

[uv](https://github.com/astral-sh/uv) is a fast, reliable Python package installer and resolver.

```bash
# Install uv if you don't have it yet
curl -LsSf https://astral.sh/uv/install.sh | sh

# Create a virtual environment
uv venv

# Activate the environment
source .venv/bin/activate  # On Linux/macOS
.venv\Scripts\activate     # On Windows

# Install required packages
uv pip install streamlit plotly pandas scipy
```

**Alternative setup with pip:**

```bash
# Create a virtual environment
python -m venv .venv

# Activate the environment
source .venv/bin/activate  # On Linux/macOS
.venv\Scripts\activate     # On Windows

# Install required packages
pip install streamlit plotly pandas scipy
```

### Usage

**Running the app locally:**

```bash
# Navigate to the directory containing uorca_explorer.py
cd UORCA/main_workflow/reporting

# Run the Streamlit app
streamlit run uorca_explorer.py
```

This will start the app on the default port (8501) and open it in your browser.

**Remote access via SSH tunnel:**

If the app is running on a remote server (e.g., HPC cluster), you can access it using an SSH tunnel:

1. Start the app on the remote server:
   ```bash
   # On remote server
   streamlit run uorca_explorer.py --server.port 8501
   ```

2. Create an SSH tunnel from your local machine:
   ```bash
   # On your local machine
   ssh -L 8000:127.0.0.1:8501 user@remote.server
   ```

3. Open `http://127.0.0.1:8000` in your browser

**Running the app in a persistent session:**

For long-running sessions on a server, use tmux or screen:

```bash
# Start a new tmux session
tmux new -s uorca-explorer

# Run the app
streamlit run uorca_explorer.py

# Detach from the session (Ctrl+B, then D)
# To reattach later:
tmux attach -t uorca-explorer
```

### Features

**Gene Selection:**
- **Auto-selected DEGs** (default): Automatically identify important differentially expressed genes using two strategies:
  - *Frequent DEGs*: Genes significant across multiple contrasts
  - *Contrast-specific DEGs*: High fold-change genes unique to few contrasts
- **Custom input**: Enter your own gene list (one gene per line) with validation and feedback
- Dynamic significance filtering with adjustable p-value and fold-change thresholds
- Up to 200 genes can be selected (automatically ranked by fold-change if more are available)

**Dataset and Contrast Selection:**
- **Dedicated "Selections" Tab**: Full-width tables for easy browsing of datasets and contrasts
- **Interactive Data Tables**: Sortable, searchable tables with checkboxes for selection
- **Synchronized Selection**: Checkboxes stay in sync across Selections, Dataset Info, and Contrast Info tabs
- **Smart Defaults**: First 5 datasets and all their contrasts selected by default
- **Automatic Updates**: Plots update automatically when selections change
- **Rich Metadata Display**: Full dataset titles, descriptions, accession numbers, and DEG counts visible without truncation

**Visualization Options:**
- Interactive heatmap for differential expression with dynamic significance filtering
- Violin plots for gene expression across samples with increased vertical space
- Customizable plot settings and pagination for large gene sets
- Real-time analysis plots (volcano, MA, heatmaps) from individual datasets
- Tab descriptions explaining the purpose of each visualization type
- **Display Settings**: Font size controls, grid line options, legend positioning
- **Professional Styling**: Consistent fonts (Inter), optimized spacing, and clean design
- **Performance Optimized**: Fragment isolation prevents unnecessary recomputation while maintaining responsive UI

**Advanced Options:**
- Adjustable significance thresholds (p-value and log fold change)
- Gene selection parameters for automatic identification
- **Per-tab Display Settings**: Customize font sizes, grid opacity, legend position
- **Violin Plot Options**: Lock y-axis across genes, toggle raw points, adjust transparency
- **Heatmap Customization**: Grid line control, font sizing, simplified contrast labels

### Performance Features

- **Optimized tab isolation**: Uses Streamlit fragments to isolate different tabs, so switching between tabs won't trigger expensive recomputations
- **Efficient caching**: The ResultsIntegrator object is cached using `@st.cache_resource`, keeping your data loaded throughout the session
- **Fragment isolation**: Expensive plot generation is isolated to prevent costly recomputation on UI changes
- **Automatic pagination**: Large gene sets are automatically paginated (30 genes per page) for optimal performance
- **Scalable selection interface**: Table-based selection handles hundreds of datasets/contrasts without UI performance issues
- **Professional styling**: Consistent fonts (Inter), optimized grid lines, and responsive design

### Troubleshooting

**Common Issues:**
- **No data loaded**: Ensure you've entered the correct path to your UORCA results directory
- **Missing datasets**: Check that your results contain properly formatted CPM.csv files
- **Slow loading**: Large datasets may take time to load initially (data is cached after first load using advanced caching)
- **Browser freezes**: The app now automatically paginates large gene selections and uses fragment isolation to prevent freezes

**Performance Tips:**
- Use the "Auto-selected DEGs" feature or custom gene input for efficient gene selection
- Use the "Selections" tab for choosing datasets and contrasts - it's optimized for large numbers of items  
- Plots update automatically when you change selections - fragment isolation keeps the interface responsive
- The app automatically paginates large gene sets and limits visualizations to prevent performance issues
- Take advantage of synchronized checkboxes across tabs for consistent selection management
- In Debug mode, enable **Show MCP logs only** to focus logging output on Multi-Core Pipeline modules

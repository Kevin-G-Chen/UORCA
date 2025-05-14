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

These findings are nice, and would generate testable hypotheses (e.g. if we induce this mutation, findings from UORCA predict that we should induce a disease state similar to X, but distinct to Y). However, I would need to
- convince researchers that any generated insights/testable hypotheses are actually worth pursuing further
- produce findings that can be more explicitly relevant for the clinic (e.g. are these DEGs meant to therefore be drug targets, or have I included DEGs for fun?)

To convince researchers that UORCA can produce valuable findings, I need to demonstrate that any "insights" reflect underlying biology, and that I can capture novel biological ideas. For this purpose, ideas of experiments I have include (dependent on project progression):
- identifying novel clinical datasets, and demonstrating that results from UORCA can be used for diagnostic purposes (e.g. that a mutation can be predictive of disease progression). I would hope this could translate into better diagnoses for diseases.
- performing docking simulations to demonstrate the possibility of repurposing a drug for treatment (e.g. if I find that gene activity for some disease is comparable to another disease, can I support this with drug-target binding predictions?)

If I can successfully establish this, then I believe there is a chance that UORCA can be used to foster collaborations to incorporate wet lab validation.

## Envisaged end result

A web server where the user
1. Inputs a research question
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

## Navigating this repository

Everything you need for running UORCA lives under the `main_workflow/` folder.  The rest of the repo (this README, CI files, docs, etc.) is mainly for project bookkeeping— day-to-day work will be in `main_workflow/`.


### Agentic workflow

1. **Master agent (`master.py`)**
   Reads your CLI inputs (accession, output folder, organism, resources) and then in turn invokes:
   - `extract()` → Extraction agent
   - `analyse()` → Analysis (and metadata) agents
   - `report()` → Reporting agent

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

Additional scripts and tools

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
2. `organism`: Organism name (e.g., human, mouse)

Additional columns can be included, but will be ignored by the script. The script will automatically generate a new CSV file with the following columns.

Example CSV content:
```
Accession,organism,RelevanceScore,Valid
GSE123456,human,8.5,Yes
GSE789012,mouse,9.2,Yes
```

You can use the `--generate-multi-csv` option in the DatasetIdentification.py script to automatically generate this CSV file.

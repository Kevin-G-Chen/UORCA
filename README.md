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
- LLM/KG mining: while I can't see using this as a sole validation method (LLM with LLM), this framework COULD help clarify any relationships between biological entities that I make. (I'd probably need to see more of it before concretely working out what would happen in reality)
- IMPPROVE: as a variant prioritization method, if I find anything about variants (i.e. a reported insight that variant X is of importance to disease Y), I could look towards performing a similar spike-in experiment as part of a validation approach.

__(note everything below here is from before)__

## Applications of work

Keeping in line with the goal of improving child health, the **generic** goals of UORCA are two-fold - to consolidate what is known in the existing literature, and to pinpoint key directions in a novel dataset, given the consolidation of the existing literature.

At time of writing this segement of the README, I am working on developing an automated pipeline that performs and RNAseq-based analysis. Once this "module" is complete, UORCA will be able to comprehensively evaluate what is known for a given research question. As a specific application I might be exploring, this will likely be "biomarkers associated with a specific cancer". Specific findings of interest would include:
- Whether these biomarkers are generically present across the literature (i.e. similar datasets)
- If not, what are factors affecting the presence/absence of these biomarkers
- Molecular/phenotypic implications (e.g. altered pathways -> disease state?)
- Drug/treatment strategies - what is proposed, do these align with the biomarkers we know from above?


## Navigating this repository

The repository is now structured into several distinct sections:

- **notebooks/**: Contains historical Jupyter notebooks used for exploratory work and prototyping.
- **archive/**: Holds previous versions and benchmarking scripts (e.g. data extraction and analysis benchmarks).
- **main/**: Will contain the core application code including web routes, services, and model definitions.
- **script_development/**: This is the primary folder for current development. All new analysis modules, agent workflows, and integration of PydanticAI (as seen in DataAnalysisAgent.py) are actively being developed here.
- **Other folders**: Such as SingleDatasetAnalysis and experiments provide additional supporting scripts or prototypes.

Overall, while multiple directories contain useful legacy or auxiliary code, the **script_development/** directory is where the active workflow work is being done.

## Intended workflow

The workflow will comprise a single master agent, which takes in some input (likely a research query - to be decided), and will delegate tasks to sub-agents. The four proposed sub-agents are as follows:
- Dataset identification agent: responsible for identifying relevant datasets
- Data extraction agent: responsible for extracting relevant data from these datasets
- Data analysis agent: responsible for performing the analysis
- Reporting agent: responsible for generating a report based on the analysis

Current progress is focussed on the DataAnalysisAgent

### Dataset analysis agent

The **DataAnalysisAgent** is a central component of the UORCA workflow, implemented using PydanticAI to drive an agentic, modular RNAseq analysis pipeline. Its design leverages dependency injection (through the `RNAseqData` model) to encapsulate paths, metadata, and runtime parameters, ensuring that each analytical step is executed consistently and reproducibly.

Key features include:

- **Modular Tool Architecture:**
  The agent is structured as a series of asynchronous tools (decorated with `@rnaseq_agent.tool`), each responsible for a specific task. These tools include:
  - **File Discovery & Management:** Tools like `list_fastq_files` and `find_files` recursively locate required FASTQ files.
  - **Data Cleaning & Metadata Processing:** Functions such as `clean_string` and `process_metadata` load, clean, and normalize metadata to ensure consistent grouping and merging of sample data.
  - **Contrast Design:** The `design_contrasts` tool automates experimental contrast creation based on merged metadata columns.
  - **Quantification & Differential Expression:**
    - `run_kallisto_quantification` executes Kallisto for RNAseq abundance estimation.
    - `prepare_edgeR_analysis` and `run_edger_limma_analysis` handle downstream differential expression analyses using established tools like edgeR and limma.
  - **Gene Set Enrichment Analysis (GSEA):** Through the `run_gsea_analysis` tool, the agent also performs pathway enrichment analyses, integrating results into the overall workflow (NOT YET TESTED)

- **Integrated Logging and Error Handling:**
  Each tool logs its progress and any encountered errors using standardized log functions (e.g., `log_tool_header` and `log_tool_result`). This ensures full traceability of the workflow execution and facilitates troubleshooting.

- **Extensibility & Future Enhancements:**
  The structure is designed to easily accommodate additional analytical steps or modifications. Future updates may include enhanced data visualization, integration with external databases, or additional omics data considerations.

Overall, the DataAnalysisAgent exemplifies a modern, agent-driven approach to complex RNAseq analyses, ensuring that as new datasets are identified, they can be seamlessly processed and integrated into the UORCA framework.

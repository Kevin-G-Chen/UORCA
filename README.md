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

We can leverage the power of LLMs to achieve the above. LLMs have been used in, for example, literature review (evidencing  LLMs' ability to collate information), as well as in assisting analysis of datasets (e.g. [Bioinformatics Agent](https://www.biorxiv.org/content/10.1101/2024.05.22.595240v2) or [scChat](https://www.biorxiv.org/content/10.1101/2024.10.01.616063v2.full)). These feature chat interfaces - not a limitation (and indeed, can be beneficial for data interaction), however the existence of such tools indicates that similar or identical workflows are not valuable to the scientific community.

__Instead, the main area where we can differentiate, gain value, and provide meaningful scientific contributions will be integrating the findings across datasets. Doing so will allow us to derive richer conclusions from datasets.__

## What specific value would UORCA provide?

UORCA will generate "insights" by collating information from different datasets. These insights might include:
- identifying genes which are differentially expressed in disease states, including those which affect disease progression/severity
- more generally - identifying common biological trends between different datasets on same/similar diseases. Alternatively - identifying biological factors underpinning different contexts
- how genetic, transcriptomic, and proteomic data relate to eachother (e.g. is the transcriptomic information correlated to the proteomic information)

These findings are nice, and would generate testable hypotheses (e.g. if we induce this mutation, findings from UORCA predict that we should induce a disease state similar to X, but distinct to Y). However, I would need to
a) convince researchers that any generated insights/testable hypotheses are actually worth pursuing further
b) produce findings that can be more explicitly relevant for the clinic (e.g. are these DEGs meant to therefore be drug targets, or have I included DEGs for fun?)

To convince researchers that UORCA can produce valuable findings, I don't see performing a lab experiment being a viable approach (once I validate I would argue this might be an appropriate follow up to any "insights," defined as per above). For this purpose, ideas of experiments I have include (dependent on project progression):
- identifying novel clinical datasets, and demonstrating that results from UORCA can be used for diagnostic purposes (e.g. that a mutation can be predictive of disease progression). I would hope this could translate into better diagnoses for diseases.
- performing docking simulations to demonstrate the possibility of repurposing a drug for treatment (e.g. if I find that gene activity for some disease is comparable to another disease, can I support this with drug-target binding predictions?)

## Envisaged end result

A web server where the user
1. Inputs a research question
2. OPTIONALLY includes their own data. The intent here would be to assist with "how does my data fit in the broader picture?"

The user should then see:
- key biological points (e.g. important genes, proteins, transcripts)
- LLM inputs/outputs, e.g. datasets analysed, chosen parameters, raw results, generated code

## Applications of work

Keeping in line with the goal of improving child health, the **generic** goals of UORCA are two-fold - to consolidate what is known in the existing literature, and to pinpoint key directions in a novel dataset, given the consolidation of the existing literature.

At time of writing this segement of the README, I am working on developing an automated pipeline that performs and RNAseq-based analysis. Once this "module" is complete, UORCA will be able to comprehensively evaluate what is known for a given research question. As a specific application I might be exploring, this will likely be "biomarkers associated with a specific cancer". Specific findings of interest would include:
- Whether these biomarkers are generically present across the literature (i.e. similar datasets)
- If not, what are factors affecting the presence/absence of these biomarkers
- Molecular/phenotypic implications (e.g. altered pathways -> disease state?)
- Drug/treatment strategies - what is proposed, do these align with the biomarkers we know from above?


## Navigating this repository

- [Notebooks](notebooks) contains the Jupyter notebooks in which I am developing the code. These are structured into smaller "sub-modules," indicated by the numbers. I will include READMEs in each directory to indicate what each file contains.
- [Data](data) will contain the results from my code. I will try to standardize the naming convention to match my code sub-modules at some point.

## Pipeline to developing the UORCA

(Work in progress)

### Dataset identification

A variety of different data sources will be considered, with some candidates being:
- ArrayExpress - functional genomics data.
- SRA (sequence read archive) - genomic information.
- ProteomeXchange/proteome central - proteomic data
- Metabolomics workbench - metabolomics data.
- dbGap - genotypes and phenotypes.
- Text - PubMed/biorXiv

The above is not exhaustive, but gives an idea of where I am currently looking.

__My focus for the time being has been with NCBI GEO__, that is utilizing transcriptomic data.

Currently, the pipeline is as follows:
1. Given a research query, use an LLM to extract keywords and identify related terms.
2. Use NCBI API to query NCBI GEO for the above keywords and related terms, identifying candidate datasets
3. Extract the title and summary associated with these datasets
4. Use an LLM to score candidate datasets, using the title and summary, based on their relevance to the research query

### Data extraction

This section currently operates as follows:
1. Download raw FASTQ files using a bash script. This entails linking the GEO datasets to sample IDs, and the SRA Run IDs linked to each sample.
2. I only download a subset of each FASTQ files for space reasons. Note that the FASTQ files are not hosted on this repository.

### Data processing - quantification

At the moment, the pipeline operates as follows:
1. Automating the Kallisto quantification, including generation of the code/identification of appropriate indices.
2. I currently focus only on human samples, however I do intend to integrate other species in - mainly those where a [pre-built index](https://github.com/pachterlab/kallisto-transcriptome-indices) is already provided.
3. After execution of the code, checking mechanisms (evaluations) are employed to check whether the quantification proceeded as expected, and if the code represents a valid method of quantification (e.g. were any warning messages biologically relevant? Were there error messages? Were parameters correctly used?)
4. Depending on the evaluation, the generated code is corrected.

(Note - I will probably rework this into developing a code "template" and having the LLM fill in the blanks)

### Using quantification data

1. Template scripts are created for:
   a. Creating the DGEList object, including filtering and normalisation with visual outputs, as well as metadata cleaning. I do hope to add PCA plots to this as well.
   b. Seeing output of cleaned metadata
   c. Generation of contrast/design matrix, and execution of DEG analysis
2. I call an LLM to determine appropriate inputs for script a and c.

I am also beginning to implement more logging and token output information - I have so far only integrated these into the first DEG analysis script, but I plan on adding this to the remaining scripts as well.

### Building the corpus

Work to begin soon. The general idea is to a) Conduct an analysis, b) Extract findings from the analysis, and c) Repeat for further datasets, however use findings (b) to contextualise these results further. This is how we can begin to see trends in datasets - do we see the same DEGs across multiple datasets? Do we notice that DEGs are only apparent in XYZ conditions?

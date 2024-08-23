# Unified -Omics Reference Corpus of Analyses (UORCA)

![image](https://github.com/user-attachments/assets/1c9f9551-0ea2-4c4a-a41c-6348e9e3ebc1)


## Objective and purpose

Biological and clinical data is rich and often holds many insights. The intepretation of these data can be complicated and time-consuming, requiring deep expertise of the subject matter and careful literature review to determine what is both novel and interesting. The advent of large language models (LLMs) provides the opportunity to expedite this otherwise resource-intensive process, in turn accelerating the rate at which we can process data and enact change. 

As such, this project aims to leverage the power of LLMs to improve the interpretation of biological data by developing a Unified -Omics Reference Corpus of Analyses, a set of pre-analysed data and results which can be utilized to both better understand the set of existing knowledge, and also uncover interesting results in novel data.

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

### Using quantification data

I will hopefully then apply the standard RNAseq pipeline, though specifics are to be discussed. 

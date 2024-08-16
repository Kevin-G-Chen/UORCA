# Unified -Omics Reference Corpus of Analyses (UORCA)
## Objective and purpose

Biological and clinical data is rich and often holds many insights. The intepretation of these data can be complicated and time-consuming, requiring deep expertise of the subject matter and careful literature review to determine what is both novel and interesting. The advent of large language models (LLMs) provides the opportunity to expedite this otherwise resource-intensive process, in turn accelerating the rate at which we can process data and enact change. 

As such, this project aims to leverage the power of LLMs to improve the interpretation of biological data by developing a Unified -Omics Reference Corpus, a set of pre-analysed data and results which can be utilized to both better understand the set of existing knowledge, and also uncover interesting results in novel data.

## Pipeline to developing the UORC

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

The plan at this stage:
1. Download raw FASTQ files, using a bash script. This entails linking the GEO datasets to sample IDs, and the SRA Run IDs linked to each sample
2. For now at least, subsetting FASTQ files - just because I can't be working FASTQ files that are too large

### Kallisto quantification

Continued from above...

3. Hopefully automating the Kallisto quantification, including generation of the code/identification of appropriate indices. For the moment, I really would like to integrate multiple species together (this is something that will be a consideration for later on). Alternatively, I could just only use human samples at the moment, and add in other species integration later.
4. Related to above - adding checking mechanisms (did this proceed as intended? If there was an error, can we use this error message to improve the generated code?)

### Using quantification data

I will hopefully then apply the standard RNAseq pipeline, though specifics are to be discussed. 

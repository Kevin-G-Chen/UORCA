# Overview

This file is intended for me to plan any experiments that I intend to be carrying out. Here, I will specify:
- the purpose of the experiment. This should be done in the context of the overall goal
- the way the experiment will be executed (i.e. methodology)
- the results to be measured
- the way the results will be interpreted

# Benchmarking performance of different LLMs for different stages of analysis pipeline for bulk RNAseq datasets

## PURPOSE

One of the fundamental components of the UORCA workflow is the analysis of data. While the main intended benefit of UORCA is the ability to collate and unify analyses from different datasets and modalities, the crux of this is a valid and reliable analysis approach.

This experiment aims to
- compare the performance of different LLMs in facilitating the automated analysis of publicly available datasets
- determine the robustness of the agentic workflow in producing a scientifically accurate workflow
- compare the performance of the agentic workflow compared to other alternatives

Upon completion of this experiment, I will have a better idea of:
- which LLMs are best at facilitating bioinformatic analyses
- if the agentic workflow is able to handle a variety of bulk RNAseq datasets
- (a specific concern) how easily more complicated experimental designs can be constructed
- the performance of my workflow as compared to alternative options:
  - manual analysis (?)
  - CrewAI
  - (This needs to be fleshed out...)

## EXPERIMENTAL DESIGN

### Inputs/necessary preparation steps

- A completed (almost) end-to-end pipeline
  - I'll be starting from the data extraction, rather than the dataset identification
- Benchmark dataset of datasets
  - i.e. the datasets that I will be using to assess performance

__Notes__: In my current iteration, I do not use AI at all in the data extraction step. However, a rework that was discussed is to instead approach it as "Here's a schematic and all the steps involved. Generate code to execute these steps".

In any case - irrespective of how the workflow is developed, I intend to use this as a way of OBJECTIVELY assessing the performance of the workflow (e.g. through quantification).

### Method

1. Run the pipeline - however, hardcode the input dataset as those included in the benchmark datasets
2. Assess performance at each stage:
  - Extraction of FASTQ files. Measurement: were all FASTQ files CORRECTLY extracted? (Y/N)
  - Extraction of metadata. Measurement: were all metadata sheets CORRECTLY extracted? (Y/N)
  - Kallisto quantification. Measurement: were the chosen parameters valid? Either: Yes, Partially (but doesn't affect results meaningfully), Partially (compromising quantification), No (e.g. completely wrong, or didn't run)
  - Metadata parsing - i.e. contrast identification and generation of contrast expressions (code for design matrix). Measurement: (how many contrast successfully identified? )


### Considerations

- If I am using this to "optimize" my workflow, I need to be careful of overfitting to the benchmark dataset.
  - I am expecting to need to make changes (e.g. at time of writing this, I don't know how long read data will be handled)
-

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
  - (edit - I am going to shift this into a different experiment)

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
  - (Note that performance should also incorporate: time taken, token usage/cost)
  - Extraction of FASTQ files. Measurement: what percentage of FASTQ files were CORRECTLY downloaded?
  - Extraction of metadata. Measurement: were all metadata sheets CORRECTLY extracted? (Y/N)
  - Kallisto quantification. Measurement: were the chosen parameters valid? (Have been given the suggestion of )
  - Metadata parsing - i.e. contrast identification and generation of contrast expressions (code for design matrix). Measurement: (how many contrast successfully identified? )
  - DEG identification. Measurement: does the generated code produce DEGs without error? Either: Yes, Partially (but doesn't affect results meaningfully), Partially (compromising results), No (e.g. completely wrong, or didn't run)

__Note__ - I think the critical part here is specifically how I choose to measure each part (i.e. what specifically am I benchmarking?)

## ANTICIPATED RESULTS/DISCUSSION POINTS

### Results

Key measurements:
- quantification of results
- (ok maybe not "key" but A measurement...) results for each "type" of dataset...?
  - the consideration to make here is do I want to comment on the types of datasets characterised by poor/strong performance?

e.g. "We assessed the performance of UORCA on different stages when using different LLMs. We found that (xyz) model had the greatest overall performance, outperforming all other tested models on all steps of the pipeline (some quotable number here...)."

e.g. "We therefore sought to use models X/Y/Z for A/B/C to compare optimal UORCA performance against other methods" (e.g. when performing follow up experiments)

### Discussion

- Do we see better performance in the more powerful models?
  - if so, this would indicate that as models improve, UORCA will continue to stay relevant

## CONSIDERATIONS

- If I am using this to "optimize" my workflow, I need to be careful of overfitting to the benchmark dataset.
  - I am expecting to need to make changes (e.g. at time of writing this, I don't know how long read data will be handled)
  - A way around this would be to include a hold-out dataset (training/testing/validation...)
- How to incorporate a statistical test?
  - I suspect the validity of the test I use will be dependent on the scoring metric I use (i.e. quantifying (?) the performance at each stage)
- What are the controls?
  - positive control - perhaps I manually determine correct parameters...?
  - negative control - ...? (I can't think of what could be applicable here...)
- The success of this is also dependent on my benchmark dataset. I have tried to capture a variety of datasets... but have I really done a good job of that?
  - I'm unsure how I'll determine how "good" my benchmark dataset is...

# Comparing performance of UORCA to other methods

(to be fleshed out)
- overall purpose is to assess (after considering results from above) does UORCA perform better than other methods for data analysis?
- my initial idea was to use the same benchmark dataset that I use for the previous experiment. However, I can kind of already see an issue here - "if you're using the benchmark dataset, which was used to optimise the performance, to comapre UORCA to other methods, how can this be fair?"

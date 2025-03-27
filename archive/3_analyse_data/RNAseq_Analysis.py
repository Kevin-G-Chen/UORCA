
# This notebook is for testing the RNAseq analysis of Kallisto quantification data.

# %% Initial preparations

##### PREPARATION

# Load modules

from openai import OpenAI
import sys
import openai # I need this and above
import os
from tqdm import tqdm
import time
import re
import numpy as np
import pandas as pd
from dotenv import load_dotenv
from pydantic import BaseModel, Field
from typing import List, Dict, Literal, Optional
import subprocess
import glob
import asyncio
import json
import base64 # image interpretation
import requests # image interpretation
import shlex # suggested for command-line strings
from datetime import datetime
from unidecode import unidecode
import tempfile

# Load .env file
load_dotenv('../../.env')
openai_api_key = os.getenv('OPENAI_API_KEY')

# %% Test API

client = OpenAI(
  api_key=openai_api_key,  # this is also the default, it can be omitted
)

chat_completion = client.chat.completions.create(
    messages=[
        {
            "role": "user",
            "content": "You are an expert of the universe, and know everything about the world. Be highly descriptive and share your knowledge with me. Share absolutely EVERYTHING that you know. For example - not just what is the universe. Also consider everything within the universe. Any topic you can think of. Again - share EVERYTHING  you know. Be as long as descriptive as you need to describe everything.",
        }
    ],
    model="gpt-4o-mini",
)

result = chat_completion.choices[0].message.content
print(result)

# %% File extraction

####### DEFINE FUNCTION FOR IDENTIFYING FILES ########

def get_files(directory, suffix):
    """
    Recursively lists all files in a given directory and its subdirectories that end with the specified suffix,
    returning their absolute paths.

    Parameters:
    directory (str): The path to the directory to search in.
    suffix (str): The file suffix to look for (e.g., 'abundance.tsv').

    Returns:
    list: A list of absolute file paths that match the given suffix.
    """
    matched_files = []
    try:
        # Walk through directory and subdirectories
        for root, _, files in os.walk(directory):
            for f in files:
                if f.endswith(suffix):
                    matched_files.append(os.path.join(root, f))

        return matched_files
    except FileNotFoundError:
        print(f"Directory '{directory}' not found.")
        return []
    except Exception as e:
        print(f"An error occurred: {e}")
        return []
#### TEST CASE ####

# Identify the metadata files for GSE133702
metadata_files = get_files(directory = "../2_extract_data/GSE133702_data",
    suffix = ".csv")
# Next identify the Kallisto quantification files
abundance_files = get_files(directory = "Benchmark_Kallisto/GSE133702_Kallisto",
    suffix = "abundance.tsv")
# Now identify the linking file
link_file = get_files(directory = "../2_extract_data/GSE133702_data",
    suffix = "sra_ids.txt")

## The files have been identified correctly. This won't be robust to superseries, but will suffice for the moment. My goal is now to load the files, i.e. the metadata and link files, link these together along with the Kallisto paths.

# Load the metadata file
metadata = pd.read_csv(metadata_files[0])
#... and the SRA ID file
SRA_IDs = pd.read_table(link_file[0])

# %% Functions to clean the metadata
##### DEVELOPING FUNCTIOSN TO PROCESS THE SAMPLE METADATA #####

##### Removing constant columns
def remove_constant_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Remove columns from the DataFrame where all values are identical.

    Parameters
    ----------
    df : pd.DataFrame
        The input DataFrame.

    Returns
    -------
    pd.DataFrame
        DataFrame with constant columns removed.
    """
    # Identify columns where all values are the same
    constant_columns = [col for col in df.columns if df[col].nunique() <= 1]

    # Drop these columns
    df_cleaned = df.drop(columns=constant_columns)

    print(f"Removed {len(constant_columns)} constant columns: {constant_columns}")

    return df_cleaned

##### Clean column names
import pandas as pd
import re

def clean_column_names(df: pd.DataFrame) -> pd.DataFrame:
    """
    Clean column names by converting to lowercase, replacing spaces with underscores,
    and removing non-alphanumeric characters.

    Parameters
    ----------
    df : pd.DataFrame
        The input DataFrame.

    Returns
    -------
    pd.DataFrame
        DataFrame with cleaned column names.
    """
    original_columns = df.columns.tolist()
    cleaned_columns = []

    for col in original_columns:
        # Convert to lowercase
        col_clean = col.lower()
        # Replace spaces and hyphens with underscores
        col_clean = re.sub(r'[ \-]+', '_', col_clean)
        # Remove non-alphanumeric characters except underscores
        col_clean = re.sub(r'[^\w]', '', col_clean)
        cleaned_columns.append(col_clean)

    # Assign cleaned column names
    df.columns = cleaned_columns

    print(f"Cleaned column names from:\n{original_columns}\nto:\n{cleaned_columns}")

    return df

#### Clean column values
import pandas as pd

def clean_column_values(df: pd.DataFrame) -> pd.DataFrame:
    """
    Clean column values to ensure compatibility with downstream analyses like edgeR/limma.
    This includes:
    - Stripping leading/trailing whitespace
    - Replacing spaces with underscores
    - Removing special characters
    - Ensuring consistent capitalization

    Parameters
    ----------
    df : pd.DataFrame
        The input DataFrame.

    Returns
    -------
    pd.DataFrame
        DataFrame with cleaned column values.
    """
    for col in df.columns:
        # Only clean columns with object dtype (typically strings)
        if df[col].dtype == 'object':
            # Strip whitespace, replace spaces with underscores, remove special chars
            df[col] = df[col].astype(str).str.strip()
            df[col] = df[col].str.replace(r'[ \-]+', '_', regex=True)
            df[col] = df[col].str.replace(r'[^\w]', '', regex=True)
            # Optionally, convert to lowercase or keep as is
            # df[col] = df[col].str.lower()

    print("Cleaned column values for object-type columns.")

    return df

# EXAMPLE USAGE
metadata_cleaned = remove_constant_columns(metadata)
metadata_cleaned = clean_column_names(metadata_cleaned)
metadata_cleaned = clean_column_values(metadata_cleaned)
metadata_cleaned

# %% Function to combine everything
##### CREATING FUNCTION TO CONSOLIDATE INTO SINGLE DATAFRAME ####

import pandas as pd
from pathlib import Path
import re
from typing import List

def combine_datasets(abundance_files: List[str],
                    metadata_df: pd.DataFrame,
                    sra_ids_df: pd.DataFrame) -> pd.DataFrame:
    """
    Combine abundance files, metadata, and SRA IDs into a single DataFrame.

    Parameters
    ----------
    abundance_files : List[str]
        List of file paths to Kallisto abundance.tsv files.

    metadata_df : pd.DataFrame
        DataFrame containing sample metadata. Must contain a column named 'geo_accession'
        which corresponds to 'sample_ID' in sra_ids_df.

    sra_ids_df : pd.DataFrame
        DataFrame containing SRA IDs. Must contain columns named 'sample_ID' and 'SRA_ID'.

    Returns
    -------
    combined_df : pd.DataFrame
        Combined DataFrame with one SRR ID per row, including metadata and abundance file path.
    """

    # Step 1: Create DataFrame from abundance_files
    abundance_df = pd.DataFrame({'abundance_file': abundance_files})

    # Extract SRR_ID from abundance_file path using regex
    # Assuming SRR IDs are in the format SRR followed by digits
    abundance_df['SRA_ID'] = abundance_df['abundance_file'].apply(
        lambda x: re.search(r'(SRR\d+)', x).group(1) if re.search(r'(SRR\d+)', x) else None
    )

    # Drop rows where SRR_ID could not be extracted
    abundance_df.dropna(subset=['SRA_ID'], inplace=True)

    # Step 2: Merge abundance_df with sra_ids_df on 'SRA_ID'
    merged_df = pd.merge(abundance_df, sra_ids_df, on='SRA_ID', how='left', validate='many_to_one')

    # Check for any SRR_IDs that did not find a match in sra_ids_df
    missing_sra = merged_df[merged_df['sample_ID'].isna()]
    if not missing_sra.empty:
        print("Warning: The following SRR_IDs from abundance_files were not found in sra_ids_df:")
        print(missing_sra['SRA_ID'].unique())
        # Optionally, you can choose to drop these rows or handle them as needed
        merged_df.dropna(subset=['sample_ID'], inplace=True)

    # Step 3: Merge with metadata_df on 'sample_ID' == 'geo_accession'
    combined_df = pd.merge(merged_df, metadata_df, left_on='sample_ID', right_on='geo_accession', how='left', validate='many_to_one')

    # Check for any sample_IDs that did not find a match in metadata_df
    missing_meta = combined_df[combined_df['title'].isna()]  # Assuming 'title' is a key column in metadata_df
    if not missing_meta.empty:
        print("Warning: The following sample_IDs from sra_ids_df were not found in metadata_df:")
        print(missing_meta['sample_ID'].unique())
        # Optionally, handle these rows as needed
        # For example, you might want to drop them
        combined_df.dropna(subset=['title'], inplace=True)

    # Optional: Reorder columns for clarity
    # Move 'abundance_file' to the front
    cols = ['abundance_file'] + [col for col in combined_df.columns if col != 'abundance_file']
    combined_df = combined_df[cols]

    # Optional: Reset index
    combined_df.reset_index(drop=True, inplace=True)

    return combined_df

# USAGE

combined_df = combine_datasets(abundance_files, metadata_cleaned, SRA_IDs)
combined_df

# %% Function to extract study summary

def get_study_metadata(accession):
    """
    Get both title and summary for a GEO study accession.

    Parameters:
    -----------
    accession : str
        The GEO accession number (e.g., "GSE133702")

    Returns:
    --------
    dict
        Dictionary containing 'title' and 'summary' of the study
    """
    # Commands to get title and summary
    title_command = (
        f'esearch -db gds -query "{accession}[ACCN]" | '
        'efetch -format docsum | '
        'xtract -pattern DocumentSummarySet -block DocumentSummary '
        f'-if Accession -equals {accession} -element title'
    )

    summary_command = (
        f'esearch -db gds -query "{accession}[ACCN]" | '
        'efetch -format docsum | '
        'xtract -pattern DocumentSummarySet -block DocumentSummary '
        f'-if Accession -equals {accession} -element summary'
    )

    # Execute commands
    title_result = subprocess.run(title_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    summary_result = subprocess.run(summary_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    # Check if commands were successful
    if title_result.returncode != 0:
        raise Exception(f"Error getting title: {title_result.stderr}")
    if summary_result.returncode != 0:
        raise Exception(f"Error getting summary: {summary_result.stderr}")

    return {
        'title': title_result.stdout.strip(),
        'summary': summary_result.stdout.strip()
    }

study_metadata = get_study_metadata("GSE133702")
print("Title:", study_metadata['title'])
print("\nSummary:", study_metadata['summary'])
# %% Now to determine appropriate contrasts to make...
###### FUNCTION TO DETERMINE CONTRASTS ######
client = OpenAI(
    api_key = openai_api_key,
)
metadata_json = combined_df.to_json(orient = 'records',
    lines=False,
    indent = 2)

class Contrast(BaseModel):
    name: str = Field(..., description = "Name of contrast to perform")
    values: list[str] = Field(..., description = "Values involved in analysis of the contrast")
    description: str = Field(..., description = "Description of the contrast")
    justification: str = Field(..., description = "Justification of why the contrast is of interest to analyse")

class AllAnalysisContrasts(BaseModel):
    contrasts: list[Contrast]


def IdentifyContrasts():
    prompt = f"""

### IDENTITY AND PURPOSE

You are an expert in bioinformatics. You advise on the most scientifically valuable experiments that can be performed, and have a deep awareness of DEG analysis tools, such as limma and edgeR.

Your task is to study the provided information, and determine what contrasts would be interesting to study.

### STEPS

1. You will be given input sample metadata. The crux of the decision making should be based on this.
2. You will be provided information about the dataset summary. Use this to inform about the scientific purpose of the dataset.
43 Having considered and digested the input information, carefully decide what the most valuable contrasts to analyse will be. Keep in mind the following guidelines:
- The values you specify should be derived ONLY from the merged column
- The contrasts you analyse should have scientific value, and not simply be "control experiments"
- The contrasts should be focussed and have a clear defined purpose

### OUTPUT

- Assign a name for each contrast
- State the values required to correctly analyse each contrast. These values must EXACTLY match the value in the metadata - do not make any changes, even if you think there are typos
- Describe what the contrast is investigating
- Justify why the contrast is scientifically valuable

### INPUTS

Study Title: {study_metadata['title']}

Study Summary: {study_metadata['summary']}

Sample Metadata: {metadata_json}

"""
    chat_completion = client.beta.chat.completions.parse(
        messages=[
            {
                "role": "user",
                "content": prompt,
            }
        ],
        model="gpt-4o",
        response_format = AllAnalysisContrasts
        )
    result = chat_completion.choices[0].message.parsed
    print(f"Generated tokens: ", chat_completion.usage.completion_tokens)
    print(f"Prompt tokens: ", chat_completion.usage.prompt_tokens)
    print(f"Total tokens: ", chat_completion.usage.total_tokens)
    return(result)

contrasts_data = IdentifyContrasts()
contrasts_data

# %% Generating expressions to analyse the contrasts
class Expressions(BaseModel):
    name: str = Field(..., description = "Name of contrast to perform")
    expressions: str = Field(..., description = "Expressions representing contrasts")

class ContrastMatrix(BaseModel):
    contrasts: list[Expressions]

def GenerateContrastExpressions():
    prompt = f"""

    ### IDENTITY AND PURPOSE

    You are an expert in bioinformatics. You advise on the most scientifically valuable experiments that can be performed, and have a deep awareness of DEG analysis tools, such as limma and edgeR.

    Your task is to study the provided information, and determine the epxressions to use to construct the contrast matrix.

    ### STEPS

    1. You will be given input information about the contrasts to use. Make note of the description of the contrast, as well as the values
    2. For each suggested contrast, state a simple name to represent it (e.g. TreatmentInKO). The fewer characters the better, however it should still be informative.
    3. For each suggested contrast, use an expression to represent it. The expression must only use values, exactly as written, indicated in the information about contrasts. Note that this expression MUST be compatible with the makeContrasts function


    ### OUTPUT

    - State a simple name for each contrast
    - State an appropriate expression for each contrast

    ### INPUTS

    Contrast information: {contrasts_data}


    """
    chat_completion = client.beta.chat.completions.parse(
        messages=[
            {
                "role": "user",
                "content": prompt,
            }
        ],
        model="gpt-4o-mini",
        response_format = ContrastMatrix
        )
    result = chat_completion.choices[0].message.parsed
    print(f"Generated tokens: ", chat_completion.usage.completion_tokens)
    print(f"Prompt tokens: ", chat_completion.usage.prompt_tokens)
    print(f"Total tokens: ", chat_completion.usage.total_tokens)
    return(result)

exprs = GenerateContrastExpressions()
exprs

# %% R script execution preparation

##### COLLECTION OF FUNCTIONS/SCRIPTS TO PREPARE FOR R EXEUCTION #####

# To prepare for the R execution, I will need a) a template R script, b) input parameters to insert into the template script
# In my original iteration, I had a single column for the analysis group. However, I don't think this will effectively scale up (e.g. mixed effects/confounding). I do need a single column for filtering purposes, so I will create a prompt to extract all this information.
metadata_json = combined_df.to_json(orient = 'records',
    lines=False,
    indent = 2)
study_metadata = get_study_metadata("GSE133702")
t2g_files = get_files(directory = "../../",
    suffix = "t2g.txt")
template_R_script_for_prompt = f"""
    library(tximport)
    library(tidyverse)
    library(edgeR)

    # Read tx2gene
    tx2gene <- read_tsv("[TX2GENE_PATH_HERE]", col_names = FALSE) %>%
      dplyr::select(1, 3) %>%
      drop_na()

    # Define abundance files
    files <- kallisto_files

    # Import data using tximport
    kallisto <- tximport(files = files,
                        type = "kallisto",
                        tx2gene = tx2gene,
                        ignoreAfterBar = TRUE,
                        countsFromAbundance = "lengthScaledTPM")

    # Read metadata
    meta <- read.csv(metadata_path, row.names = 1)

    # Create DGEList
    DGE <- DGEList(counts = kallisto$counts,
                  samples = meta)

    keep.exprs <- filterByExpr(DGE, group = DGE$samples$[COLUMN_FOR_FILTERING_HERE])
    DGE.filtered <- DGE[keep.exprs, keep.lib.sizes = FALSE]
    print(dim(DGE.filtered))
    # Normalize
    DGE.final <- calcNormFactors(DGE.filtered)

"""

prompt = f"""

### IDENTITY AND PURPOSE

You are an expert bioinformatician, and have a deep understanding of RNAseq analysis pipelines, including the code necessary to execute such pipelines as well as the biological implications of different chosen parameters. You will be given information about a dataset, including the input metadata, and a template R script. You will be asked to generate parameters to facilitate an R-based RNAseq analysis.

### STEPS

1. Carefully digest and interpret the provided sample metadata. Ensure you
2. Carefully digest and interpret the information associated with the dataset. This information will help contextualise what the dataset was studying.
3. Take your time to study and analyse the provided R script. This R script contains placeholder values, which will be the values that you will need to provide. Therefore, ensure you carefully consider what the purpose of each placeholder parameter is - what is the biological significance of each, and what is the required format for each entry?
4. Having absorbed all of the above information, indicate what the most appropriate value, or values, is for each placeholder. Keep in mind the following guidelines for each parameter:
    - tx2gene path: here, ensure you select a SINGLE t2g.txt file that is most appropriate. The distinguishing factor is the species - so ensure that the species is accurate based on any information that is provided
    - filter_analysis_group: this is a SINGLE COLUMN that will be used to filter lowly expressed genes. The chosen column must EXACTLY match a column name in the metadata, and should be a column that most accurately represents different "experimental groups"

### OUTPUT

Your response will include:
    - the EXACT parameter that should be inputted for the tx2gene path
    - the EXACT parameter that should be inputted for the filter_analysis_group
    - the justification for each of the above

### INPUTS

Template R script: {template_R_script_for_prompt}
Information about study: {study_metadata}
Possible tx2gene paths: {t2g_files}
Sample metadata: {metadata_json}
"""

class RNAseq_Preparation_Parameters(BaseModel):
    tx2gene: str = Field(..., description = "Path to the appropriate tx2gene file")
    filter_column: str = Field(..., description = "Name of the column to use for filtering")
    tx2gene_justification: str = Field(..., description = "Justification for why the tx2gene file is correct for this dataset")
    filter_column_justification: str = Field(..., description = "Justification for why the chosen column is appropriate for use in filtering")


chat_completion = client.beta.chat.completions.parse(
    messages=[
        {
            "role": "user",
            "content": prompt,
        }
    ],
    model="gpt-4o",
    response_format = RNAseq_Preparation_Parameters
    )
result = chat_completion.choices[0].message.parsed
print(f"Generated tokens: ", chat_completion.usage.completion_tokens)
print(f"Prompt tokens: ", chat_completion.usage.prompt_tokens)
print(f"Total tokens: ", chat_completion.usage.total_tokens)
print(result)
# %% Use identified parameters to perform RNAseq

# Need to specify the different parameters - including: tx2gene (from above), column (from above), Kallisto files (need code), and metadata (need code)

t2g_path = result.tx2gene
kallisto_paths = combined_df['abundance_file']
filter_column = result.filter_column

with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.csv', encoding='utf-8') as tmp_meta:
    metadata_path = tmp_meta.name
    combined_df.to_csv(metadata_path, index=False)

with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.R', encoding='utf-8') as tmp_r_script:
    r_script_path = tmp_r_script.name

    # Start constructing the R script
    r_script = f"""
    library(tximport)
    library(tidyverse)
    library(edgeR)

    # Read tx2gene
    tx2gene <- read_tsv("{t2g_path}", col_names = FALSE) %>%
      dplyr::select(1, 3) %>%
      drop_na()

    # Define abundance files
    files <- c({', '.join([f'"{file}"' for file in kallisto_paths])})

    # Import data using tximport
    kallisto <- tximport(files = files,
                        type = "kallisto",
                        tx2gene = tx2gene,
                        ignoreAfterBar = TRUE,
                        countsFromAbundance = "lengthScaledTPM")

    # Read metadata
    meta <- read.csv("{metadata_path}", row.names = 1)

    # Create DGEList
    DGE <- DGEList(counts = kallisto$counts,
                  samples = meta)

    keep.exprs <- filterByExpr(DGE, group = DGE$samples${filter_column})
    DGE.filtered <- DGE[keep.exprs, keep.lib.sizes = FALSE]
    print(dim(DGE.filtered))
    # Normalize
    DGE.final <- calcNormFactors(DGE.filtered)
    str(DGE.final)
    """

    tmp_r_script.write(r_script)

print("Generated R Script:\n")
with open(r_script_path, 'r') as f:
    print(f.read())

# Execute the R script
try:
    print("Executing R script...")
    subprocess.run(["Rscript", r_script_path], check=True)
    print("R script executed successfully.")
except subprocess.CalledProcessError as e:
    print("An error occurred while executing the R script.")
    print("Error message:")
    print(e.stderr if e.stderr else e)
finally:
    # Clean up temporary files if desired
    os.remove(metadata_path)
# %% Construction of design and contrast matrix

# When I did this before I only assumed a single columnm, which allowed for a simple experimental design (and a simple construction of the design matrix). However, this will not scale very well for other experimental designs.

# And so some thought needs to be put on the best way to approach this. Earlier I generated "contrasts to analyse" - I think I will instead combine this as part of a larger prompt where I simulataneously determine the contrasts to analyse and how to analyse them. My justification here is that by combining them together, it will be easier to determine what contrasts are feasible to analyse (in comparison, when doing them separately, I might be trying to construct a design matrix for a comparison that isn't feasible. Or something like that...)


class Contrast(BaseModel):
    name: str = Field(..., description = "Name of contrast to perform")
    expression: list[str] = Field(..., description = "The exact expression used to represent the contrast in the makeContrasts function")
    description: str = Field(..., description = "Description of the contrast")
    justification: str = Field(..., description = "Justification of why the contrast is of interest to analyse")

class ContrastAnalyses(BaseModel):
    contrasts: list[Contrast]
    model_matrix: str = Field(..., description = "Code used to design the model matrix")
    model_matrix_justification: str = Field(..., description = "Explanation of the construction of model matrix")

template_R_script_for_prompt = """

# Create the model matrix
design <- model.matrix([INSERT_CODE_FOR_TEMPLATE_MATRIX_HERE], data = DGE.final$samples)
print(design)

# Create contrast matrix
contrast.matrix <- makeContrasts(
{CONTRAST_EXPRESSIONS_GO_HERE},
levels = colnames(design)
    )


v <- voom(DGE.final,
          design)
vfit <- lmFit(v,
              design)

vfit <- contrasts.fit(vfit,
contrast.matrix)

efit <- eBayes(vfit)

contrasts <- colnames(contrast.matrix)

LFC.summary <- sapply(contrasts, function(x){{
lfc.list <- list()
top <- topTable(efit,
                    coef = x,
                    number = Inf) %>%
    list()
    lfc.list <- append(lfc.list, top)
    }})

str(LFC.summary)

"""

prompt = f"""
## IDENTITY AND PURPOSE

You are an expert bioinformatician and specialise in conducting RNAseq analyses. You provide guidance on how to produce the most accurate and informative analyses for a wide variety of datasets, particularly as they apply to pipelines incorporating packages such as DESeq2/limma/edgeR.

You will be provided information about datasets, the samples within these datasets, and some template code. Your task is to use the provided information to determine what are the most appropriate inputs to complete the template code. Follow the provided steps carefully to produce the best possible outcome.

## STEPS
The first set of steps relates to understanding the context of the answers you will provide

1. Carefully digest and interpret the dataset summary. You should use this information to work out what where the data is derived from, and to make determinations about what is scientifically interesting to study
2. Carefully digest and interpret the sample metadata. Use this information to determine the SPECIFIC variables that can be studied, and to determine what is actually possible to analyse.
- Note this this, in combination with the dataset summary, should be the only information you use. Do not make any assumptions about the data that is not explicitly stated
3. Carefully consider the provided template R code. Note that your answers pertaining to code will need to fit into the template code.

Having considered how your answers must conform to the requirements, you will be now generating responses.

4. Determine what contrasts to analyse based on the sample metadata and dataset summary. Take into account the following recommendations:
    - The contrasts should be scientifically interesting
    - The contrasts should be feasibly studied (i.e. able to be implemented using code, taking into consideration the structure of the metadata)

5. For each contrast, provide a name for the contrast, the code that should be used to represent the contrast (i.e. the expression), and a justification for why the contrast should be studied
6. After constructing these contrasts, consider what appropriate model matrix code should be implemented to enable study of the identified contrasts. Note that the model matrix must be possible based on the provided metadata, with no additional processing. Therefore, only columns present in the metadata should be used - additional columns cannot be created.
7. If an appropriate model matrix cannot be constructed, revise the contrasts until a model matrix can be constructed.

## OUTPUT

Ensure your output conforms to the ContrastAnalyses class:

class Contrast(BaseModel):
    name: str = Field(..., description = "Name of contrast to perform")
    expression: list[str] = Field(..., description = "The exact expression used to represent the contrast in the makeContrasts function")
    description: str = Field(..., description = "Description of the contrast")
    justification: str = Field(..., description = "Justification of why the contrast is of interest to analyse.")

class ContrastAnalyses(BaseModel):
    contrasts: list[Contrast]
    model_matrix: str = Field(..., description = "Code used to design the model matrix")
    model_matrix_justification: str = Field(..., description = "Explanation of the construction of model matrix")

Note that for model_matrix, the output should emulate EXACTLY what will go in the template - do not copy any of the surrounding information

## INPUTS

Template R script: {template_R_script_for_prompt}
Information about study: {study_metadata}
Sample metadata: {metadata_json}

"""

chat_completion = client.chat.completions.create(
    messages=[
        {
            "role": "user",
            "content": prompt,
        }
    ],
    model="o1-mini",
)

result = chat_completion.choices[0].message.content
print(result)
prompt2 = f"""

You will be given a prompt that was passed to an LLM, and the output produced by the LLM. Your task is to convert the output produced by the LLM into the required structured output.

Original prompt: {prompt}
Output: {result}
"""

chat_completion = client.beta.chat.completions.parse(
    messages=[
        {
            "role": "user",
            "content": prompt2,
        }
    ],
    model="gpt-4o",
    response_format = ContrastAnalyses
    )
result = chat_completion.choices[0].message.parsed
print(f"Generated tokens: ", chat_completion.usage.completion_tokens)
print(f"Prompt tokens: ", chat_completion.usage.prompt_tokens)
print(f"Total tokens: ", chat_completion.usage.total_tokens)
print(result)
# %% Use results to execute code
t2g_path
# %% test

import tempfile
import os
import subprocess
import pandas as pd
from typing import List

def run_dge_analysis_with_contrasts(
    analysis,          # A ContrastAnalyses instance
    t2g_path,          # Path to tx2gene
    combined_df,       # Pandas DataFrame of sample metadata
    filter_column,     # Column name to group by in filterByExpr
    r_library_paths: List[str] = None,  # Optional: specify .libPaths() if needed
):
    """
    Generate and run an R script using the specified ContrastAnalyses object
    to create a design matrix and contrasts on-the-fly.

    :param analysis: ContrastAnalyses object with model_matrix and contrasts
    :param t2g_path: Path to the tx2gene file
    :param combined_df: Pandas DataFrame of sample metadata
    :param filter_column: Column name from metadata to be used in filterByExpr
    :param r_library_paths: Optionally a list of R library paths to prepend in .libPaths()
    """

    # 1) Convert the model_matrix string (e.g. "~0 + title")
    #    into something we can directly insert into R code.
    model_matrix_code = analysis.model_matrix  # e.g. "~0 + title"

    # 2) Build the "makeContrasts(...)" argument for each contrast.
    #
    #    Each Contrast has:
    #      - name (string)
    #      - expression (list[str])  --> Usually a single expression string in your example,
    #                                    but let's join them if there are multiple lines.
    #    We'll produce something like:
    #      `HighAltitude_vs_Basal` = ( ... ) - ( ... ),
    #      `Day7_vs_Day3_HighAltitude_Kyrgyz` = ( ... ) - ( ... ),
    #      ...
    #
    #    Note: the backticks around {c.name} allow for safe R variable naming
    #          even if the user-supplied name has special characters.
    contrast_statements = []
    for c in analysis.contrasts:
        # Join multiple expression lines into one (if you have multi-line).
        expression_str = " ".join(c.expression)
        # Build something like:  `ContrastName` = (expr1) - (expr2)
        stmt = f"`{c.name}` = {expression_str}"
        contrast_statements.append(stmt)

    # Join them with commas and newlines for readability
    contrast_code = ",\n".join(contrast_statements)

    # 3) Create temp files for metadata CSV and the R script
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.csv', encoding='utf-8') as tmp_meta:
        metadata_path = tmp_meta.name
        combined_df.to_csv(metadata_path, index=False)

    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.R', encoding='utf-8') as tmp_r_script:
        r_script_path = tmp_r_script.name

        # If you need custom library paths, build an R command like:
        #   .libPaths(c("/some/custom/path", .libPaths()))
        # This is optional. Only do this if you have special library locations.
        libpaths_code = ""
        if r_library_paths:
            # e.g. r_library_paths = ["/path/to/R/libs1", "/path/to/R/libs2"]
            joined_paths = '", "'.join(r_library_paths)
            libpaths_code = f'.libPaths(c("{joined_paths}", .libPaths()))'

        # 4) Construct the final R script:
        r_script = f"""
{libpaths_code}

library(tximport)
library(tidyverse)
library(edgeR)
library(limma)

# Read tx2gene
tx2gene <- read_tsv("{t2g_path}", col_names = FALSE) %>%
  dplyr::select(1, 3) %>%
  drop_na()

# Define abundance files
files <- c({', '.join(f'"{file}"' for file in combined_df['abundance_file'])})

# Import data using tximport
kallisto <- tximport(files = files,
                    type = "kallisto",
                    tx2gene = tx2gene,
                    ignoreAfterBar = TRUE,
                    countsFromAbundance = "lengthScaledTPM")

# Read metadata
meta <- read.csv("{metadata_path}", row.names = 1)

# Create DGEList
DGE <- DGEList(counts = kallisto$counts,
               samples = meta)

keep.exprs <- filterByExpr(DGE, group = DGE$samples${filter_column})
DGE.filtered <- DGE[keep.exprs, keep.lib.sizes = FALSE]
print(dim(DGE.filtered))

# Normalize
DGE.final <- calcNormFactors(DGE.filtered)
str(DGE.final)

########################################################
## The user-specified design matrix from ContrastAnalyses
########################################################
design <- model.matrix({model_matrix_code}, data = DGE.final$samples)
colnames(design) <- str_remove_all(colnames(design), "title")
print(design)
########################################################
## Create contrast matrix
########################################################
contrast.matrix <- makeContrasts(
{contrast_code},
levels = colnames(design)
)
print("Hello")

########################################################
## Fit models and compute statistics
########################################################
v <- voom(DGE.final, design)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrast.matrix)
efit <- eBayes(vfit)

contrasts <- colnames(contrast.matrix)

LFC.summary <- sapply(contrasts, function(x) {{
    top <- topTable(efit, coef = x, number = Inf)
    return(list(top))
}})

str(LFC.summary)

# Write out the full results table for each contrast
# (optional if you want to save them to disk)
# for (cnt in contrasts) {{
#   out.path <- paste0("{os.path.dirname(r_script_path)}", "/", cnt, "_results.csv")
#   write.csv(topTable(efit, coef = cnt, number = Inf), out.path, row.names = FALSE)
# }}

"""

        # Write the R script and close
        tmp_r_script.write(r_script)

    # 5) Print or log the generated R script for debugging
    print("Generated R Script:\n")
    with open(r_script_path, 'r') as f:
        print(f.read())

    # 6) Execute the R script in a subprocess
    try:
        print("Executing R script...")
        subprocess.run(["Rscript", r_script_path], check=True)
        print("R script executed successfully.")
    except subprocess.CalledProcessError as e:
        print("An error occurred while executing the R script.")
        print("Error message:")
        print(e.stderr if e.stderr else e)
    finally:
        # Clean up temporary files if you like
        os.remove(metadata_path)
        os.remove(r_script_path)

# %% usage

result_temp = ContrastAnalyses(contrasts=[Contrast(name='HighAltitude_vs_Basal', expression=['(High_altitude_Day_3_4111m_Kyrgyz_Sample_7 + High_altitude_Day_3_4111m_Kyrgyz_Sample_8 + High_altitude_Day_7_4111m_Kyrgyz_Sample_11 + High_altitude_Day_7_4111m_Kyrgyz_Sample_12 + High_Altitude_Day_14_4111m_Kyrgyz_Sample_15 + High_altitude_Day_14_4111m_Kyrgyz_Sample_16 + High_Altitude_Day_21_4111m_Kyrgyz_Sample_19 + High_Altitude_Day_21_4111m_Kyrgyz_Sample_20 + High_Altitude_day_21_4111m_Indian_Sample_17 + High_Altitude_day_21_4111m_Indian_Sample_18 + High_Altitude_day_21_4111m_Indian_Sample_18) / 11 - (Basal__800m_Kyrgyz_Sample_3 + Basal__800m_Kyrgyz_Sample_4) / 2'], description='Compare gene expression levels between combined high altitude time points (Day 3, Day 7, Day 14, Day 21) and basal level (800m) samples.', justification='To identify genes that are differentially expressed due to altitude-induced stress compared to basal conditions at sea level, thus assessing the biological variations occurring due to high altitude exposure across time.'), Contrast(name='Day7_vs_Day3_HighAltitude_Kyrgyz', expression=['(High_altitude_Day_7_4111m_Kyrgyz_Sample_11 + High_altitude_Day_7_4111m_Kyrgyz_Sample_12) / 2 - (High_altitude_Day_3_4111m_Kyrgyz_Sample_7 + High_altitude_Day_3_4111m_Kyrgyz_Sample_8) / 2'], description='Compare gene expression between Day 7 and Day 3 high altitude samples in Kyrgyz males.', justification='To investigate specific temporal gene expression changes that occur as a response to continued high altitude exposure, providing insights into adaptive responses over time in Kyrgyz males.'), Contrast(name='Indian_vs_Kyrgyz_HighAltitude', expression=['(High_altitude_Day_7_4111m_Indian_Sample_10 + High_altitude_day_14_4111m_Indian_Sample_14 + High_Altitude_day_21_4111m_Indian_Sample_17) / 3 - (High_altitude_Day_7_4111m_Kyrgyz_Sample_11 + High_altitude_Day_14_4111m_Kyrgyz_Sample_15 + High_Altitude_Day_21_4111m_Kyrgyz_Sample_19) / 3'], description='Compare gene expression at high altitudes (4111m) between Indian and Kyrgyz males over different days.', justification='To explore the inter-ethnic genetic response to high altitude, identifying variations that could contribute to differences in high-altitude adaptation between these two ethnic groups.')], model_matrix='~0 + title', model_matrix_justification="The 'title' variable within the sample metadata is leveraged as it uniquely categorizes and labels each sample. This design allows contrasts based on the experimental setup and specific features like altitude, sample day, and nationality to be systematically compared.")

result_temp = ContrastAnalyses(
    contrasts=[
        Contrast(
            name='HighAltitude_vs_Basal',
            expression=[
                # Summation of Day 3, 7, 14, 21 (Kyrgyz + Indian) / 10 minus Basal / 2
                "((High_altitude_Day_3_4111m_Kyrgyz_Sample_7 + "
                "High_altitude_Day_3_4111m_Kyrgyz_Sample_8 + "
                "High_altitude_Day_7_4111m_Kyrgyz_Sample_11 + "
                "High_altitude_Day_7_4111m_Kyrgyz_Sample_12 + "
                "High_Altitude_Day_14_4111m_Kyrgyz_Sample_15 + "
                "High_Altitude_Day_14_4111m_Kyrgyz_Sample_16 + "
                "High_Altitude_Day_21_4111m_Kyrgyz_Sample_19 + "
                "High_Altitude_Day_21_4111m_Kyrgyz_Sample_20 + "
                "High_Altitude_day_21_4111m_Indian_Sample_17 + "
                "High_Altitude_day_21_4111m_Indian_Sample_18) / 10) - "
                "((Basal__800m_Kyrgyz_Sample_3 + Basal__800m_Kyrgyz_Sample_4) / 2)"
            ],
            description=(
                "Compare gene expression levels between combined high-altitude "
                "time points (Day 3, Day 7, Day 14, Day 21) and basal (800m) samples."
            ),
            justification=(
                "To identify genes differentially expressed due to altitude-induced "
                "stress compared to basal conditions, assessing biological variations "
                "across time."
            )
        ),
        Contrast(
            name='Day7_vs_Day3_HighAltitude_Kyrgyz',
            expression=[
                # Compare Day 7 vs. Day 3 Kyrgyz
                "((High_altitude_Day_7_4111m_Kyrgyz_Sample_11 + "
                "High_altitude_Day_7_4111m_Kyrgyz_Sample_12) / 2) - "
                "((High_altitude_Day_3_4111m_Kyrgyz_Sample_7 + "
                "High_altitude_Day_3_4111m_Kyrgyz_Sample_8) / 2)"
            ],
            description=(
                "Compare gene expression between Day 7 and Day 3 high-altitude "
                "samples in Kyrgyz males."
            ),
            justification=(
                "To investigate specific temporal gene-expression changes in Kyrgyz "
                "males as a response to continued high-altitude exposure."
            )
        ),
        Contrast(
            name='Indian_vs_Kyrgyz_HighAltitude',
            expression=[
                # Compare Indian vs. Kyrgyz across (some) high-altitude time points
                "((High_altitude_Day_7_4111m_Indian_Sample_10 + "
                "High_altitude_Day_14_4111m_Indian_Sample_14 + "
                "High_Altitude_day_21_4111m_Indian_Sample_17) / 3) - "
                "((High_altitude_Day_7_4111m_Kyrgyz_Sample_11 + "
                "High_Altitude_Day_14_4111m_Kyrgyz_Sample_15 + "
                "High_Altitude_Day_21_4111m_Kyrgyz_Sample_19) / 3)"
            ],
            description=(
                "Compare gene expression at high altitudes (4111m) between Indian "
                "and Kyrgyz males over different days."
            ),
            justification=(
                "To explore inter-ethnic differences in high-altitude adaptation "
                "by comparing genetic responses in Indian vs. Kyrgyz males."
            )
        ),
    ],
    model_matrix='~0 + title',
    model_matrix_justification=(
        "The 'title' variable in the sample metadata uniquely labels each sample. "
        "This design supports contrasts based on altitude, day, and nationality."
    )
)

run_dge_analysis_with_contrasts(
    analysis=result_temp,
    t2g_path=t2g_path,
    combined_df=combined_df,
    filter_column=filter_column
)

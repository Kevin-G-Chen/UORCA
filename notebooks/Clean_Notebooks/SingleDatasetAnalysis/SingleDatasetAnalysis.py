#!/usr/bin/env python
# coding: utf-8

# # Purpose
# 
# This notebook will contain ONLY the code to complete the end to end workflow. The idea is to make sure I can execute everything from start to finish, and to understand where the challenges will lie.
# 
# I am also finding that a bit of the code with structured outputs in the past is now outdated. I will modernise it.
# 
# I will mark some of the hardcoded stuff with [HARDCODED]. (Ooh, fancy formatting. I hope that doesn't break anything...)

# # Part 0 - Loading modules/.env/preparing OpenAI client
# 
# The goal of this section is to load all the necessary modules, as well as prepare the OpenAI client.
# 
# Hmmm... I suspect in the end I would want all the classes/functions/prompts defined somewhere...
# 

# In[122]:


# Load modules

from openai import OpenAI
import os
import json
from tqdm import tqdm
import time
from Bio import Entrez
import numpy as np
import pandas as pd
from dotenv import load_dotenv
import instructor
from pydantic import BaseModel, Field
from typing import  List, Dict, Literal, Optional
import asyncio
from concurrent.futures import ThreadPoolExecutor
from collections import Counter
import statistics
import subprocess
from unidecode import unidecode
import re
import pandas as pd
import subprocess
import tempfile
import os
import json


# Prepare .env file

load_dotenv('../../.env') # [HARDCODED]

Entrez.email = os.getenv('ENTREZ_EMAIL')
Entrez.api_key = os.getenv('ENTREZ_API_KEY')
openai_api_key = os.getenv('OPENAI_API_KEY')

# Prepare OpenAI API Client

client = OpenAI(
  api_key=openai_api_key,  # this is also the default, it can be omitted
)


# # Part 1 - Dataset identification
# 
# The goal from this section is to identify the datasets that will be analysed. 
# 
# Keep in mind - for this notebook, I will only be analysing one dataset.
# 
# ## Specific approach
# - A user query is given. This is customizable, and will be what the user supplies themselves (i.e. the question they are interested in)
# - Based on this user query, appropriate search terms are identified using AI. Perhaps an option to specify how many iterations are performed?
# - Datasets are extracted from this. For the moment, I only extract the first 20 datasets (hardcoded). I suppose this could be a parameter?
# - Based on information from the extracted datasets, the relevance of the datasets to the research question is determined. This is performed three times. It'd be good to specify how many iterations are performed.

# In[34]:


# Prepare the initial search queries

# Define some initial variables for demonstration purposes [HARDCODED]

user_query = "Identify datasets and samples which are relevant to exploring immunotherapies for lung cancer"
num_queries = 3

# 

# Prepare functions for term extraction from research query

# Prepare output structures
## For term extraction from user query
class ExtractedTerms(BaseModel):
    extracted_terms: List[str] = Field(description="List of terms extracted from the query")
    expanded_terms: List[str] = Field(description="List of related terms generated from the extracted terms")

## For determining dataset relevance
class Assessment(BaseModel):
    ID: str
    RelevanceScore: int = Field(description="Score from 0 to 10, indicating relevance")
    Justification: str = Field(description="A brief explanation for the score")

class Assessments(BaseModel):
    assessments: List[Assessment]

# Define function for term extraction 
def extract_terms(user_query: str) -> List[str]:
    prompt = f"""

## IDENTITY AND PURPOSE
You are an expert in literature searches of biological ideas. Your task is to identify biological term(s) from a query, and generate related terms for the purposes of generating a search query. 

## STEPS

- First, extract the biological term(s) from the input query. These should be specific and fall into one of the following categories:
1. Genes - Examples: BRCA1, TP53
2. Treatments/Methods - Examples: chemotherapy, CRISPR
3. Tissues/Cells - Examples: lung, hepatocytes
4. Diseases - Examples: Alzheimer's disease, lung cancer.

Do not fabricate items if no relevant term exists. Avoid general terms such as "disease" or "variant."

- Second, for each extracted biological term, generate two related terms. Make a considered effort to keep these terms in the same category as the original term. These are examples of an identified term, and possible relevant terms:
1. Genes: BRCA1 - Examples: BRCA2, oncogene
2. Treatments: Chemotherapy - Examples: radiotherapy, monoclonal antibody
3. Tissues/Cells: Lung - Examples: respiratory, alveoli
4. Diseases: Alzheimer's disease - Examples: dementia, amyloid plaque

## OUTPUT

Provide two lists:
1. Extracted terms: The primary terms identified directly from the query.
2. Expanded terms: The related terms generated from the extracted terms.
Do not include categories or justifications.

## INPUT
User query: {user_query}"""
    
    extracted_terms = client.beta.chat.completions.parse(
        model="gpt-4o-mini",
        temperature=1,
        response_format=ExtractedTerms,
        messages=[
            {"role": "user", 
             "content": prompt}
        ]
    )
    extracted_terms = extracted_terms.choices[0].message.parsed

#    print(f"Raw extracted terms: {extracted_terms.extracted_terms}")
#    print(f"Raw expanded terms: {extracted_terms.expanded_terms}")
    
    all_terms = extracted_terms.extracted_terms + extracted_terms.expanded_terms
    terms_with_filter = [term + ' AND "gse"[Filter]' for term in all_terms]
    return terms_with_filter

# Extension - define function to perform term extraction multiple times
async def extract_terms_multiple(user_query: str, num_queries: int = 3) -> List[str]:
    async def single_extract():
        return extract_terms(user_query)
    
    tasks = [single_extract() for _ in range(num_queries)]
    results = await asyncio.gather(*tasks)
    
    # Flatten the list of lists and remove duplicates
    all_terms = list(set([term for sublist in results for term in sublist]))
    return all_terms

# Define function for performing search
def perform_search(term):
    search_handle = Entrez.esearch(db="gds", term=term, retmode="xml", retmax = 50) # CHANGE
    search_results = Entrez.read(search_handle)
    search_handle.close()
    return search_results

# Define function for extracting information from above search results
def extract_geo_info_batch(geo_ids):
    """
    Retrieve GEO information for a batch of GEO IDs.
    """
    ids_str = ",".join(geo_ids)
    handle = Entrez.esummary(db="gds", id=ids_str, retmode="xml")
    output = Entrez.read(handle)
    handle.close()

    data = []
    for geo_id, geo_data in zip(geo_ids, output):
        if isinstance(geo_data, dict):
            data.append({
                'ID': geo_id,
                'Title': geo_data.get('title', 'No title available'),
                'Summary': geo_data.get('summary', 'No summary available'),
                'Accession': geo_data.get('Accession', 'No accession available'),
                'Species': geo_data.get('taxon', 'No taxon available'),
                'Date': geo_data.get('PDAT', 'Date made public unknown')
            })
        else:
            data.append({'ID': geo_id, 'Title': 'Error', 'Summary': 'Unable to fetch data', 'Accession': 'Error'})

    return data

def create_geo_dataframe(geo_ids, batch_size=10):
    """Create a DataFrame from GEO search results using batch processing."""
    data = []
    for i in tqdm(range(0, len(geo_ids), batch_size), desc="Processing GEO IDs in batches"):
        batch_ids = geo_ids[i:i + batch_size]
        data.extend(extract_geo_info_batch(batch_ids))
        time.sleep(0.2)  # Be nice to NCBI servers
    return pd.DataFrame(data)


# Define function for determining relevance of datasets
def assess_relevance_batch(df, query, batch_size=10):
    results = []
    total_batches = (len(df) + batch_size - 1) // batch_size
    for i in tqdm(range(0, len(df), batch_size), desc="Determining dataset relevance", total=total_batches):
        batch = df.iloc[i:i+batch_size]
        prompt = f"""
## IDENTITY AND PURPOSE

You are a highly knowledgeable biologist tasked with identifying relevant datasets for a given research query. Your goal is to assess NCBI GEO datasets based on their titles and summaries, and determine their relevance to the research question at hand.

## STEPS

1. For each dataset, carefully analyze the provided title and summary.
2. Extract ALL biological concepts represented in the dataset, including but not limited to:
   - Genes and variants investigated (e.g., p53, BRCA1)
   - Species studied (e.g., Homo sapiens, Escherichia coli)
   - Sample sources (e.g., organoid cultures, human samples)
   - Diseases or phenotypes studied (e.g., Alzheimer's disease, lung cancer)
   - Cell types or tissues examined (e.g., lung tissue, neural progenitor cells)
   - Experimental techniques or methodologies used (e.g., RNA-seq, ChIP-seq)
3. Extract ALL biological concepts represented in the research query using the same categories.
4. Assign a relevance score from 0 to 10 in increments of 1, based solely on the provided information. 
   - Do not fabricate or assume information not explicitly stated about the dataset. 
   - If confirmed information about the dataset is POSSIBLY useful for the research question, view this favourably for determining dataset relevance. 
   - Note that the gene, disease, and cell type/tissue being studied are the most important in determining relevance. The other factors are considered minor aspects.
   - Use the following scoring guide:

   0: No relevance. All biological concepts (genes, species, samples, diseases, cell types, methods) are completely unrelated to the research query.
   1: Minimal relevance. One minor aspect is loosely related, but the overall focus is different.
   2: Low relevance. One major aspect aligns with the query, but other key elements differ significantly.
   3: Somewhat low relevance. Two aspects align, but critical elements are still mismatched.
   4: Moderate relevance. Multiple aspects align, but there are still significant differences in focus or approach.
   5: Moderately relevant. Most major aspects align, but there are some notable differences that may limit direct applicability.
   6: Relevant. All major aspects align, but there might be differences in specific genes, cell types, or methodologies that somewhat reduce direct applicability.
   7: Highly relevant. Very close alignment in all major aspects, with only minor differences that don't significantly impact applicability.
   8: Very highly relevant. Near-perfect alignment in all major aspects, with at most one or two minor differences.
   9: Extremely relevant. Perfect alignment in all major aspects, with at most one negligible difference.
   10: Perfectly relevant. The dataset appears to be an exact match for the research query in all aspects.

5. Provide a brief justification (3-4 sentences) for the assigned score, highlighting key similarities and differences.

## OUTPUT
For each dataset, provide a JSON object with the ID, relevance score, and justification. 
- The relevance score should be a number, with no other information.
- The justification should be a 4-5 sentence explanation for the relevance score. 

## HANDLING LIMITED INFORMATION
If the dataset title or summary lacks sufficient detail:
- Focus on the information that is available
- Do not make assumptions about missing information
- Assign a lower score if critical information is absent
- Note the lack of information in the justification

Remember, it's better to assign a lower score due to lack of information than to assume relevance without evidence.

Given the following datasets and query, determine if each dataset is relevant.
        Query: {query}
        Datasets:
        """
        for _, row in batch.iterrows():
            prompt += f"""
            ID: {row['ID']}
            Title: {row['Title']}
            Summary: {row['Summary']}
            Species: {row['Species']}
            """
        
        try:
            response = client.beta.chat.completions.parse(
                model="gpt-4o-mini",
                temperature=0.3,
                response_format=Assessments,
                messages=[
                    {"role": "user", 
                     "content": prompt}
                ],
                max_tokens=10000
            )
            response = response.choices[0].message.parsed
            results.extend([assessment.dict() for assessment in response.assessments])
        except Exception as e:
            results.extend([{"ID": row['ID'], "Relevance": "Error", "Justification": str(e)} for _, row in batch.iterrows()])
        time.sleep(1)  # Be nice to the API
    return results

# Extension - define function to assess relevance multiple times
async def assess_relevance_batch_multiple(df, query, num_queries: int = 2, batch_size=20):
    async def single_assess():
        return assess_relevance_batch(df, query, batch_size)
    
    tasks = [single_assess() for _ in range(num_queries)]
    results = await asyncio.gather(*tasks)
    
    # Collate results
    collated_results = {}
    for i, result_set in enumerate(results):
        for assessment in result_set:
            id = assessment['ID']
            if id not in collated_results:
                collated_results[id] = {'scores': [], 'justifications': []}
            collated_results[id]['scores'].append(assessment['RelevanceScore'])
            collated_results[id]['justifications'].append(assessment['Justification'])
    
    # Determine final relevance and format output
    final_results = []
    for id, data in collated_results.items():
        mean_score = statistics.mean(data['scores'])
        
        result = {
            'ID': id,
            'RelevanceScore': round(mean_score, 1),
        }
        
        # Add individual scores and justifications
        for i in range(num_queries):
            result[f'IndividualScore{i+1}'] = data['scores'][i] if i < len(data['scores']) else None
            result[f'Justification{i+1}'] = data['justifications'][i] if i < len(data['justifications']) else None
        
        final_results.append(result)
    
    return final_results

async def main(user_query):
    # Extract terms
    search_terms = await extract_terms_multiple(user_query)
    print("Search terms:", search_terms)

    # Perform Entrez search and remove duplicates
    geo_ids = set()  # Use a set to automatically remove duplicates
    for term in search_terms:
        search_results = perform_search(term)
        geo_ids.update(search_results.get('IdList', []))  # Update the set with new IDs
    if not geo_ids:
        return pd.DataFrame({'Error': ["No results found for the extracted terms"]})

    # Convert set back to list
    geo_ids = list(geo_ids)[1:20] # for the moment only use a subset of the IDs [HARDCODED]

    # Create DataFrame with GEO information
    df = create_geo_dataframe(geo_ids)

    # Assess relevance
    relevance_results = await assess_relevance_batch_multiple(df, user_query, num_queries=num_queries) # Currently have this at 3
    relevance_df = pd.DataFrame(relevance_results)

    # Merge results
    df['ID'] = df['ID'].astype(str)
    relevance_df['ID'] = relevance_df['ID'].astype(str)
    result_df = df.merge(relevance_df, on='ID', how='left')

    # Dynamically create the desired order of columns
    base_columns = ['ID', 'Title', 'Summary', 'Species', 'Accession', 'Date', 'RelevanceScore']
    score_columns = [f'IndividualScore{i+1}' for i in range(num_queries)]
    justification_columns = [f'Justification{i+1}' for i in range(num_queries)]
    desired_order = base_columns + score_columns + justification_columns

    # Reorder columns
    result_df = result_df[desired_order]

    # Reset index
    result_df = result_df.reset_index(drop=True)

    return result_df


# In[36]:


# View results

dataset_relevance_df = await (main(user_query))
dataset_relevance_df


# # Part 2 - Data extraction
# 
# In this part, relevant data is extracted out of the NCBI GEO datasets.
# 
# For this notebook, I will extract data only out of the best scoring dataset (i.e. highest relevance).
# 
# This part relies on scripts that I generated in a different notebook, so I will make sure to call those.
# 
# ## Specific approach:
# - Determine appropriate input parameters for script (i.e. dataset ID, output directory name, number of "spots"/reads)
# - Download the metadata associated with the dataset ID (which is a NCBI GEO ID)
# - Download FASTQ files associated with the dataset ID

# In[56]:


# [HARDCODED] For the moment, we will begin by determining the single dataset that we should analyse

top_accession = dataset_relevance_df.sort_values(by="RelevanceScore", ascending=False).iloc[0]["Accession"]

# We will then use this to determine input parameters. I think I am happy leaving these hardcoded.
output_dir_name = top_accession + "_data"

n_spots = 80000 # [HARDCODED]

script_dir = "OtherScripts"  # Adjust if necessary

# Construct the path to the process_geo.sh script
script_path = os.path.join(script_dir, "process_geo.sh")

# Run the subprocess
result = subprocess.run([
    script_path,
    "--geo_accession", top_accession,
    "--output_dir", output_dir_name,
    "--num_spots", str(n_spots),
    "--force"
], check=True)  # check=True will raise an exception if the subprocess fails


# Post performance notes:
# - There was a period of time where fastq-dump would not work (network issues - I could not isolate whether it was an issue on NCBI's end or my end). However, in other time periods it seems to work very well.
# - I will be using the sra_ids.txt file to link the various IDs and whatnot. This is necessary for the next step.

# # Part 3 - Data analysis
# 
# Now that the data has been extracted, we will now want to perform the analysis. 
# 
# The specifics steps involved here are:
# ## (3.1 - Kallisto quantification)
# - View the documentation
# - Identify the file locations (FASTQ files and Kallisto index files). In my original iteration, I had the file locations hardcoded, so I do need to determine how to resolve this... my vision was the have the AI determine appropriate values (that way I can rely on just a single function... (Hm, I think I can get away with the FASTQ files by specifying the output directory as defined earlier, but I don't have the same luxury for the index files. Perhaps I just search from the home directory...)
# - Identify sample metadata. I believe this is the SRA metadata... I need to iron this out, because this might be a duplicate (and perhaps it might be more suitable for the previous part...)
# - Get the study summary
# - Use the above information to determine the appropriate Kallisto parameters
# - Using these determined parameters, perform the Kallisto quantification
# 
# ## 3.2 - DEG analysis
# - Read in sample metadata (as extracted from before)
# - Identify location of abundance files
# - Determine appropriate contrasts from metadata, and structure this appropriately (i.e. compatible with makeContrasts)
# - Perform the DEG analysis, with the input files/contrasts
# 
# Later - I would need an evaluation mechanism for.. pretty much every step. I'm hoping I can simply develop something which goes "hey check these steps" and can be flexible beyond that.

# In[57]:


# Start by getting the documentation. This is necessary to ensure the OpenAI API knows the versions etc. that are being dealt with.

def get_documentation(command):
    try:
        # Execute the kallisto command
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        
        # Capture the stdout
        stdout = result.stdout
        
        # Return the results
        return stdout
    except Exception as e:
        print(f"An error occurred: {e}")
        return None

kallisto_docs = get_documentation("kallisto quant --help")


# In[89]:


# Now determine the file locations. This gets me the FASTQ files, but is also necessary to get the SRA IDs, which I use to extract the SRA metadata.

def list_files(directory, suffix, exclude_hidden=True):
    """
    Recursively lists all files in a given directory and its subdirectories that end with the specified suffix,
    optionally excluding hidden files and directories, returning their absolute paths.

    Parameters:
    directory (str): The path to the directory to search in.
    suffix (str): The file suffix to look for (e.g., 'fastq.gz').
    exclude_hidden (bool): If True, hidden files and directories are excluded. Defaults to True.

    Returns:
    list: A list of absolute file paths that match the given suffix.
    """
    matched_files = []
    
    try:
        for root, dirs, files in os.walk(directory):
            if exclude_hidden:
                # Skip hidden directories
                dirs[:] = [d for d in dirs if not d.startswith('.')]
            
            for f in files:
                if exclude_hidden and f.startswith('.'):
                    continue
                if f.endswith(suffix):
                    matched_files.append(os.path.join(root, f))
                    
        return matched_files
    except FileNotFoundError:
        print(f"Directory '{directory}' not found.")
        return []
    except Exception as e:
        print(f"An error occurred: {e}")
        return []
# First extract the FASTQ files. I can use the directory defined above to get these.

# output_dir_name = top_accession + "_data" This is how I defined this earlier

fastq_directory = output_dir_name
fastq_suffix = ".fastq.gz" # [HARDCODED] Hoping that I can automate this
fastq_files = list_files(fastq_directory, fastq_suffix)

# Next is the Kallisto indices...

index_directory = "/home/myuser/work/" # [HARDCODED] .. but maybe this is ok...
index_suffix = ".idx" # [HARDCODED]
index_files = list_files(index_directory, index_suffix)

# Now we extract the study summary

def get_study_summary(accession):

    # Define the command as a string
    command = (
        f'esearch -db gds -query "{accession}[ACCN]" | '
        'efetch -format docsum | '
        'xtract -pattern DocumentSummarySet -block DocumentSummary '
        f'-if Accession -equals {accession} -element summary'
    )

    # Execute the command
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    # Check if the command was successful
    if result.returncode == 0:
        # Return the output
        return result.stdout.strip()
    else:
        # Raise an error with the stderr output
        raise Exception(f"Error: {result.stderr}")

study_summary = get_study_summary(top_accession)

# Now extract the SRA sample metadata...

def fetch_sra_metadata_shell(
    sra_ids_file,
    entrez_api_key=None,
    delay=None,
    verbose=True
):
    """
    Fetches metadata for a list of SRA IDs using NCBI's esearch and efetch command-line tools.

    Parameters:
    - sra_ids_file (str): Path to the file containing SRA IDs with headers (tab-separated).
    - entrez_api_key (str, optional): NCBI Entrez API key. If not provided, it will be read from the
      'ENTREZ_API_KEY' environment variable.
    - delay (float, optional): Delay in seconds between requests to respect rate limits.
      If not provided, it defaults to 0.5 seconds without an API key and 0.1 seconds with an API key.
    - verbose (bool, optional): If True, prints progress messages. Defaults to True.

    Returns:
    - pandas.DataFrame: Combined DataFrame containing metadata for all fetched SRA IDs.
    """
    # Set Entrez API key from parameter or environment variable
    if entrez_api_key is None:
        entrez_api_key = os.getenv('ENTREZ_API_KEY')
    
    if entrez_api_key:
        Entrez.api_key = entrez_api_key
        if verbose:
            print("API key detected.")
    else:
        if verbose:
            print("No API key detected; proceeding without it.")
    
    # Set default delay based on API key presence
    if delay is None:
        delay = 0.1 if entrez_api_key else 0.5
    
    # List to store each SRA ID's fetched data as a DataFrame
    data = []
    
    # Check if the SRA IDs file exists
    if not os.path.isfile(sra_ids_file):
        raise FileNotFoundError(f"SRA IDs file '{sra_ids_file}' does not exist.")
    
    # Open and read the SRA IDs file using csv.DictReader for tab-separated values
    with open(sra_ids_file, 'r', newline='') as ids_file:
        reader = csv.DictReader(ids_file, delimiter='\t')
        if 'SRA_ID' not in reader.fieldnames:
            raise ValueError("Input file must contain a 'SRA_ID' column.")
        
        for line_num, row in enumerate(reader, start=2):  # start=2 accounts for header
            sra_id = row.get('SRA_ID', '').strip()
            if not sra_id:
                if verbose:
                    print(f"Line {line_num}: Missing 'SRA_ID'. Skipping.")
                continue  # Skip if SRA_ID is missing
        
            if verbose:
                print(f"\nProcessing SRA ID: {sra_id}")
        
            # Construct the command
            command = f"esearch -db sra -query {sra_id} | efetch -format runinfo"
        
            if verbose:
                print(f"Executing command: {command}")
        
            try:
                # Execute the command
                result = subprocess.run(
                    command,
                    shell=True,
                    capture_output=True,
                    text=True,
                    check=True
                )
        
                # Check if output is not empty
                if not result.stdout.strip():
                    if verbose:
                        print(f"No data returned for SRA ID: {sra_id}.")
                    continue
        
                # Convert the CSV output to a DataFrame
                csv_data = StringIO(result.stdout)
                df = pd.read_csv(csv_data)
                data.append(df)
        
                if verbose:
                    print(f"Successfully fetched data for SRA ID: {sra_id}.")
        
            except subprocess.CalledProcessError as e:
                print(f"Error processing {sra_id} on line {line_num}: {e}")
                print(f"Command output: {e.output}")
                continue  # Skip to the next SRA ID if there’s an error
        
            # Respect API rate limits
            if verbose:
                print(f"Sleeping for {delay} seconds to respect rate limits.")
            time.sleep(delay)
    
    # Combine all DataFrames into one
    if data:
        combined_df = pd.concat(data, ignore_index=True)
        
        # Remove columns where all entries are NaN
        combined_df.dropna(axis=1, how='all', inplace=True)
    
        if verbose:
            print("\nData fetching complete.")
        return combined_df
    else:
        if verbose:
            print("No data was fetched.")
        return pd.DataFrame()  # Return empty DataFrame
sra_file = list_files(fastq_directory,
                      "sra_ids.txt")

sra_metadata = fetch_sra_metadata_shell(sra_file[0])


# In[104]:


# Above is getting a bit too chunky, so next section.

# This is where I'll do some Kallisto stuff.

class KallistoCommand(BaseModel):
    index: str = Field(..., description="Filename for the Kallisto index to be used for quantification")
    fastq1: str = Field(..., description="Filename for the first FASTQ file (Read 1) to be quantified")
    fastq2: Optional[str] = Field(description="Filename for the second FASTQ file (Read 2) to be quantified (optional for single-end reads)")
    output: str = Field(..., description="Directory to write output to")
    bootstraps: int = Field(..., description="Number of bootstrap samples")
    single: bool = Field(..., description="If the reads are single-end")
    fr_stranded: bool = Field(..., description="If the reads are strand-specific, with first read forward")
    rf_stranded: bool = Field(..., description="If the reads are strand-specific, with first read reverse")
    frag_length: Optional[int] = Field(description="Estimated average fragment length (required for single-end reads)")
    sd: Optional[int] = Field(description="Estimated standard deviation of fragment length (required for single-end reads)")
    justification: str = Field(..., description="Justification for each chosen parameter, including if the parameter was excluded")

class KallistoCommands(BaseModel):
    commands: List[KallistoCommand] = Field(description="List of Kallisto quantification commands for each sample")

def identify_kallisto_params():
    prompt = f"""

## IDENTITY AND PURPOSE

You are an expert in bioinformatic analyses. You will be provided with various pieces of information, and use this information to determine the appropriate parameters for a Kallisto analysis.

## STEPS

1. Carefully digest the contents of the provided Kallisto documentation. Note that any existing knowledge you have of Kallisto may not be correct, so follow the documentation closely.
2. Carefully consider the contents of the sample metadata. Not all information will be relevant, however there will be content that will be needed.
3. Carefully look through the dataset metadata. This may contain details that are useful.
4. After considering all of the above, determine which Kallisto parameters should be set. Do not make any assumptions that are not explicitly stated for any optional fields. If unsure, leave blank.
5. In determining parameters, make sure you only choose valid files (i.e. pick out of the options which are provided)
6. Ensure that the chosen parameters allow for a robust analysis that would satisfy the most critical peer reviewers.
7. You should prioritize scientific robustness over ease of computational burden.
8. Note the following guidelines for some specific parameters:
- the output directory should be named such that the sample being quantified can be identified from this output directory.

## OUTPUT

Your output should consist of each parameter, and either:
- the value to be included for the parameter
- if the parameter should not be included, you should state NA
- For ALL chosen parameters, describe the justification for including the particular value, or excluding it.

This should be applied to all parameters identified as per the provided Kallisto documentation.

## INPUT

Kallisto documentation: {kallisto_docs}

Dataset summary: {study_summary}

FASTQ files: {fastq_files}

Possible Kallisto indices: {index_files}

Sample metadata: {sra_metadata.to_json}

"""
    chat_completion = client.beta.chat.completions.parse(
        messages=[
            {
                "role": "user",
                "content": prompt,
            }
        ],
        model="gpt-4o-mini",
        response_format = KallistoCommands
        )
    result = chat_completion.choices[0].message.parsed
    print(f"Generated tokens: ", chat_completion.usage.completion_tokens)
    print(f"Prompt tokens: ", chat_completion.usage.prompt_tokens)
    print(f"Total tokens: ", chat_completion.usage.total_tokens)
    return(result)

kallisto_params = identify_kallisto_params()

for cmd in kallisto_params.commands:
    # Construct the Kallisto command string
    kallisto_cmd = f"kallisto quant -i {cmd.index} -o {cmd.output} -t 4"
    
    if cmd.bootstraps > 0:
        kallisto_cmd += f" --bootstrap-samples={cmd.bootstraps}"
    
    if cmd.single:
        kallisto_cmd += " --single"
        if cmd.frag_length:
            kallisto_cmd += f" -l {cmd.frag_length}"
        if cmd.sd:
            kallisto_cmd += f" -s {cmd.sd}"
    else:
        # Paired-end
        if cmd.fr_stranded:
            kallisto_cmd += " --fr-stranded"
        elif cmd.rf_stranded:
            kallisto_cmd += " --rf-stranded"
    
    # Append FASTQ files
    kallisto_cmd += f" {cmd.fastq1} {cmd.fastq2}"

def execute_kallisto_commands(kallisto_commands: KallistoCommands):
    for cmd in kallisto_commands.commands:
        # Construct the Kallisto command string
        kallisto_cmd = f"kallisto quant -i {cmd.index} -o {cmd.output} -t 4 --plaintext"
        
        if cmd.bootstraps > 0:
            kallisto_cmd += f" --bootstrap-samples={cmd.bootstraps}"
        
        if cmd.single:
            kallisto_cmd += " --single"
            if cmd.frag_length:
                kallisto_cmd += f" -l {cmd.frag_length}"
            if cmd.sd:
                kallisto_cmd += f" -s {cmd.sd}"
        else:
            # Paired-end
            if cmd.fr_stranded:
                kallisto_cmd += " --fr-stranded"
            elif cmd.rf_stranded:
                kallisto_cmd += " --rf-stranded"
        
        # Append FASTQ files
        if cmd.fastq2 and cmd.fastq2.lower() != 'na':
            kallisto_cmd += f" {cmd.fastq1} {cmd.fastq2}"
        else:
            kallisto_cmd += f" {cmd.fastq1}"
        
        print(f"Executing Kallisto command for {cmd.fastq1}:")
        print(kallisto_cmd)

        # Execute the command
        try:
            subprocess.run(kallisto_cmd, shell=True, check=True)
            print(f"Kallisto quantification completed for {cmd.fastq1}\n")
        except subprocess.CalledProcessError as e:
            print(f"Error executing Kallisto for {cmd.fastq1}: {e}\n")
        
        # Optionally, log the justification
        justification_path = os.path.join(cmd.output, "justification.txt")
        os.makedirs(cmd.output, exist_ok=True)
        with open(justification_path, "w") as f:
            f.write(cmd.justification)
        print(f"Justification saved to {justification_path}\n")

if __name__ == "__main__":
    kallisto_commands = identify_kallisto_params()
    execute_kallisto_commands(kallisto_commands)


# ## Part 3b - Analysing the quantification data

# In[118]:


# Start by reading the metadata... I hope this won't lead to any complications...

metadata_csv = list_files(directory = fastq_directory, 
                          suffix = ".csv")
df = pd.read_csv(metadata_csv[0])
df = df.loc[:, df.nunique() > 1]
metadata_json = df.to_json(orient='records', lines=False, indent=2) # parse to JSON

# Quick inspection - will be interesting, because I previously had an example reliant on merging. This one does not. Will see if I've hardcoded

class ColumnMerging(BaseModel):
    merge: bool = Field(..., description="Whether or not columns should be merged")
    cols: Optional[list[str]] = Field(..., description="List of columns to be merged")
    justification: str = Field(..., description = "Justification of columns being merged/why no columns needed to be merged")

def Identify_ColMerges():
    prompt = f"""

### IDENTITY AND PURPOSE

You are an expert in bioinformatics. You advise on the most scientifically valuable experiments that can be performed, and have a deep awareness of DEG analysis tools, such as limma and edgeR.

Your task is to study the provided metadata, and determine which columns to use in proceeding with the analysis.

### STEPS

Note that a future step of the analysis will involve design of a matrix as follows:
design <- model.matrix(data = DGE.final$samples,
                       ~0 + column)

Crucially, this only includes a single column. As such, if there are columns with DISTINCT SCIENTIFIC INFORMATION, these should be merged. Columns with similar information DO NOT need to be merged. Therefore, take a deep breath, and follow these steps to ensure that subsequent analyses are as robust as possible:

1. Assess the content of each column in the provided metadata
2. Determine which columns contain anything of biological relevance
3. Determine if any columns are redundant, and do not need to be considered (e.g. similar content). In this case, only consider the column with simpler values (i.e. fewer special characters)
4. Determine which columns contain information that would be scientifically valuable to analyse, i.e. could result in a meaningful biological finding.
5. If there are multiple columns that contain scientifically valuable information, identify these columns as needing to be merged.
6. If there is one one column containing scientifically valuable information, no columns need to be merged
7. If you would be merging two redundant columns, these do not need to be merged. As such, no merge should occur (i.e. set merge to FALSE). Note that in this case, merging will COMPLICATE the analysis. Instead, IGNORING one of these columns is the best way to proceed.
8. Be very aware that no merging can be perfectly viable. Do not force a suboptimal merge.

Take into consideration that, suppose the values in one column are
A
B
C

And another column are 
1
2
3

The merged output would be
A_1
B_2
C_3 

(or something comparable to that)

### OUTPUT

- Specify if any columns will need to be merged
- State the names of the columns to be merged
- Justify your choice

### INPUT METADATA

{metadata_json}

"""
    chat_completion = client.beta.chat.completions.parse(
        messages=[
            {
                "role": "user",
                "content": prompt,
            }
        ],
        model="gpt-4o-mini",
        response_format = ColumnMerging
        )
    result = chat_completion.choices[0].message.parsed
    print(f"Generated tokens: ", chat_completion.usage.completion_tokens)
    print(f"Prompt tokens: ", chat_completion.usage.prompt_tokens)
    print(f"Total tokens: ", chat_completion.usage.total_tokens)
    return(result)

col_merge_info = Identify_ColMerges()
col_merge_info


# In[126]:


# omg I actually give up. I think a completely different approach is needed. Now, I don't think this will change any of the results, however it is indicative of this not working to the standard I want.

def clean_string(s: str) -> str:
    """
    Clean a string by normalizing special characters, replacing spaces with underscores,
    and removing non-word characters.

    Args:
        s (str): The string to clean.

    Returns:
        str: The cleaned string.
    """
    if pd.isnull(s):
        return "NA"  # Handle missing values
    s = str(s)
    s = s.strip()  # Remove leading and trailing whitespaces
    s = unidecode(s)  # Normalize special characters to ASCII
    s = s.replace(" ", "_")  # Replace spaces with underscores
    s = re.sub(r'[^\w]', '', s)  # Remove non-word characters (retain letters, digits, underscores)
    return s

def process_column_merging(df: pd.DataFrame, column_merge_info: ColumnMerging) -> pd.DataFrame:
    """
    Process column merging based on ColumnMerging information.

    Args:
        df (pd.DataFrame): The sample metadata DataFrame.
        column_merge_info (ColumnMerging): Information about column merging.

    Returns:
        pd.DataFrame: The updated DataFrame with merged columns if applicable.
    """
    if column_merge_info.merge:
        # Ensure that at least two columns are provided for merging
        if not column_merge_info.cols or len(column_merge_info.cols) < 2:
            raise ValueError("At least two columns must be specified for merging when merge=True.")
        
        cols_to_merge = column_merge_info.cols

        # Generate new column name by combining base names of the columns to merge
        # For example, merging 'genotype:ch1' and 'treatment:ch1' becomes 'genotype_treatment_clean'
        base_names = [col.split(":")[0] for col in cols_to_merge]
        new_col_name = "merged_analysis_group"

        # Clean the values in the columns to be merged
        cleaned_columns = df[cols_to_merge].map(clean_string)

        # Merge the cleaned columns by concatenating their values with underscores
        df[new_col_name] = cleaned_columns.apply(lambda row: "_".join(row.values), axis=1)

        print(f"Merged columns {cols_to_merge} into '{new_col_name}'.")
    else:
        # When merging is not required, ensure exactly one column is specified
        if not column_merge_info.cols or len(column_merge_info.cols) != 1:
            raise ValueError("Exactly one column must be specified for cleaning when merge=False.")
        
        col_to_clean = column_merge_info.cols[0]

        # Generate a new column name by appending '_clean' to the original column name
        new_col_name = "merged_analysis_group"

        # Rename the column in the DataFrame
        df = df.rename(columns={col_to_clean: new_col_name})

        # Clean the values in the renamed column
        df[new_col_name] = df[new_col_name].apply(clean_string)

        print(f"Cleaned column '{col_to_clean}' into '{new_col_name}'.")

    return df

cleaned_metadata_df = process_column_merging(df, col_merge_info)
cleaned_metadata_json = cleaned_metadata_df.to_json(orient='records', lines=False, indent=2)


# In[128]:


# Now to identify the contrasts...

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
2. You will be given some input information about a "merged column" called "merged_analysis_group". You should focus on the values in this column. However, the information will also detail where the merged values are derived from, so you can use this information as well.
3. You will be provided information about the dataset summary. Use this to inform about the scientific purpose of the dataset.
4. Having considered and digested the input information, carefully decide what the most valuable contrasts to analyse will be. Keep in mind the following guidelines:
- The values you specify should be derived ONLY from the merged column
- The contrasts you analyse should have scientific value, and not simply be "control experiments"
- The contrasts should be focussed and have a clear defined purpose
- Here are some examples of how to structure the contrasts:
    - If the samples to be compared are, for example "Treatment X vs. Y in genotpye A samples", the output should be "X_A, Y_A" (where X_A refers to the EXACT value in the merged_analysis_group column)
    - If the samples to be compared are, for example "Treatment X vs. Y", the output should be "X_A, X_B, Y_A, Y_B". 
5. Once you have produced the output, double check that:
- You have considered the correct column
- The values you have stated are derived from the correct column


### OUTPUT

- Assign a name for each contrast
- State the values required to correctly analyse each contrast. These values must EXACTLY match the value in the merged_analysis_group column
- Describe what the contrast is investigating
- Justify why the contrast is scientifically valuable

### INPUTS

Sample metadata: {cleaned_metadata_json}
Information about merged columns: {col_merge_info}
Dataset summary: {study_summary}


"""
    chat_completion = client.beta.chat.completions.parse(
        messages=[
            {
                "role": "user",
                "content": prompt,
            }
        ],
        model="gpt-4o-mini",
        response_format = AllAnalysisContrasts
        )
    result = chat_completion.choices[0].message.parsed
    print(f"Generated tokens: ", chat_completion.usage.completion_tokens)
    print(f"Prompt tokens: ", chat_completion.usage.prompt_tokens)
    print(f"Total tokens: ", chat_completion.usage.total_tokens)
    return(result)

contrasts_data = IdentifyContrasts()


# In[193]:


# Generate expressions for these contrasts

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
3. For each suggested contrast, use an expression to represent it. The expression must only use values, exactly as written, indicated in the information about contrasts. Note that this expression MUST be compatible with the makeContrasts function. See below for some examples:
"GNASknockout - WT"
"(GNASknockout_A - GNASknockout_B) - (WT_A - WT_B)"


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


# In[197]:


# The contrast (singular in this case, I would agree this is reasonable) to analyse has been determined, now to generate DGEList object to perform the DEG analysis.

# Part 1 of this is to identify the input files that will be needed. This includes the metadata, linking samples to their abundance files, 
# the actual abundance files, transcript to gene text file

# First the Kallisto and tx2gene files.

abundance_directory = "." # [HARDCODED]
abundance_suffix = "abundance.tsv" # [HARDCODED]
abundance_files = list_files(abundance_directory, abundance_suffix) # just for my own sanity I didn't print the output, but I can see it was able to find all the files
tx2gene_files = list_files(directory = "/home/myuser/work", # Ah. [HARDCODED]. This won't work well...
                          suffix = "t2g.txt")

SRA_IDs = pd.read_csv(sra_file[0], sep = '\t')


# In[189]:


class IDMatching(BaseModel):
    SRA_ID: str = Field(..., description = "Name of SRA ID")
    Kallisto_path: str = Field(..., description = "Name of matching Kallisto path")

class AllIDMatches(BaseModel):
    AllMatches: list[IDMatching]

def Match_SRAIDs():
    prompt = f"""

You will be given inputs for SRA IDs, as well as the path to abundance files generated from Kallisto. Your task is to generate 1-to-1 matches between abundance files and the SRA IDs.

That is - for each SRA ID, identify the single path that is most likely to correspond to that SRA ID.

Your output should consist of ONLY the SRA ID, and their matching Kallisto path. The output should match the input text EXACTLY with no other formatting included.

### INPUTS

SRA IDs: {SRA_IDs['SRA_ID']}
Kallisto paths: {abundance_files}


"""
    chat_completion = client.beta.chat.completions.parse(
        messages=[
            {
                "role": "user",
                "content": prompt,
            }
        ],
        model="gpt-4o-mini",
        response_format = AllIDMatches
        )
    result = chat_completion.choices[0].message.parsed
    print(f"Generated tokens: ", chat_completion.usage.completion_tokens)
    print(f"Prompt tokens: ", chat_completion.usage.prompt_tokens)
    print(f"Total tokens: ", chat_completion.usage.total_tokens)
    return(result)

# Now a commands to link these together...

SRA_ID_links = Match_SRAIDs()
SRA_IDs = pd.read_csv(sra_file[0], sep = '\t')
json_data = SRA_ID_links.json()
all_id_matches = AllIDMatches.parse_raw(json_data)
sra_to_kallisto = {match.SRA_ID: match.Kallisto_path for match in all_id_matches.AllMatches}
SRA_IDs['Kallisto_path'] = SRA_IDs['SRA_ID'].map(sra_to_kallisto)


# In[199]:


# With this being done, now I prepare the reading of the files and the execution of the R script


# Paths and data (Assuming these are defined elsewhere in your code)
tx2gene_path = tx2gene_files[1] # [HARDCODED] ... I think I need this as an LLM prompt to determine which path to use... 
analysis_group = "merged_analysis_group" # [HARDCODED] ...

# Export metadata to a temporary CSV file for R to read
with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.csv', encoding='utf-8') as tmp_meta:
    metadata_path = tmp_meta.name
    linked_data.to_csv(metadata_path, index=False)

# Create a temporary R script file
with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.R', encoding='utf-8') as tmp_r_script:
    r_script_path = tmp_r_script.name

    # Start constructing the R script
    r_script = f"""
    library(tximport)
    library(tidyverse)
    library(edgeR)
    
    # Read tx2gene
    tx2gene <- read_tsv("{tx2gene_path}", col_names = FALSE) %>%
      dplyr::select(1, 3) %>%
      drop_na()
    
    # Define abundance files
    files <- c({', '.join([f'"{file}"' for file in abundance_files])})
    
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
    
    keep.exprs <- filterByExpr(DGE, group = DGE$samples${analysis_group})
    DGE.filtered <- DGE[keep.exprs, keep.lib.sizes = FALSE]
    print(dim(DGE.filtered))
    # Normalize
    DGE.final <- calcNormFactors(DGE.filtered)
    """

    # Add Design Matrix Code
    r_script += f"""
    library(stringr)
    
    # Create design matrix using the specified grouping variable
    design <- model.matrix(~0 + {analysis_group}, data = DGE.final$samples)
    
    # Clean column names by removing the grouping variable string
    colnames(design) <- str_remove_all(colnames(design), "{analysis_group}")
    
    print(design)
    """

    # Add Contrast Matrix Code
    # Extract contrast names and expressions from exprs
    contrast_entries = []
    for contrast in exprs.dict()['contrasts']:
        name = contrast['name']
        expression = contrast['expressions']
        # Escape double quotes in expressions
        expression = expression.replace('"', '\\"')
        contrast_entries.append(f'{name} = "{expression}"')

    contrast_matrix_str = ",\n  ".join(contrast_entries)

    # Use single quotes in message to avoid conflicts with double quotes in contrast_matrix_str
    r_script += f"""
    colnames(design)
    # Create contrast matrix
    contrast.matrix <- makeContrasts(
      {contrast_matrix_str},
      levels = colnames(design)
    )
    
    
    # Optionally, you can proceed with fitting the model and other downstream analysis
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

    saveRDS(LFC.summary, "LFCs.RDS")
    """

    # Write the complete R script to the temporary file
    tmp_r_script.write(r_script)

# Optional: Print the generated R script for debugging
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


# In[ ]:





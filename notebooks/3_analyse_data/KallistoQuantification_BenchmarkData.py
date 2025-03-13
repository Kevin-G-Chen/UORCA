################################################################################
# KALLISTO QUANTIFICATION OF BENCHMARK DATA
################################################################################

# %% Preparation

####################
# IMPORTS
####################
from openai import OpenAI
import openai # Probably don't need above... but this is for testing tools with structured outputs
import re
import os
from tqdm import tqdm
import time
import numpy as np
import pandas as pd
from dotenv import load_dotenv
from pydantic import BaseModel, Field
from typing import List, Dict, Literal, Optional
from pathlib import Path
import subprocess
import glob
import asyncio
import json
from Bio import Entrez
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
from concurrent.futures import ThreadPoolExecutor

# Set number of CPU cores to use
N_CORES = max(1, multiprocessing.cpu_count() - 1)  # Leave one core free
BOOTSTRAP_SAMPLES = 5  # Set consistent number of bootstrap samples

# Prepare OpenAI API

load_dotenv('../../.env')

Entrez.email = os.getenv('ENTREZ_EMAIL')
Entrez.api_key = os.getenv('ENTREZ_API_KEY')
openai_api_key = os.getenv('OPENAI_API_KEY')

####################
# CONFIGURATION
####################

# Set number of CPU cores to use
N_CORES = max(1, multiprocessing.cpu_count() - 1)  # Leave one core free
BOOTSTRAP_SAMPLES = 5  # Set consistent number of bootstrap samples

####################
# DOCUMENTATION FUNCTIONS
####################

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



# Get files
def get_files(directory, suffix):
    """
    Recursively lists all files in a given directory and its subdirectories that end with the specified suffix,
    returning their absolute paths.

    Parameters:
    directory (str): The path to the directory to search in.
    suffix (str): The file suffix to look for (e.g., 'fastq.gz').

    Returns:
    list: A list of absolute file paths that match the given suffix.
    """
    matched_files = []

    try:
        # Expand user home directory and make path absolute
        abs_dir = os.path.abspath(os.path.expanduser(directory))

        if not os.path.exists(abs_dir):
            return []

        # Walk through directory and subdirectories
        for root, _, files in os.walk(abs_dir):
            for f in files:
                if f.endswith(suffix):
                    full_path = os.path.join(root, f)
                    matched_files.append(full_path)

        return matched_files
    except FileNotFoundError:
        print(f"Directory '{directory}' not found.")
        return []
    except Exception as e:
        print(f"An error occurred: {e}")
        return []



# Function to extract study summary

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




# %% Prepare the input file
# The input file will contain the GSE accessions of the benchmark datasets
benchmark_file = Path("../../data/BenchmarkDatasets/BulkRNAseq.txt")

# Read benchmark accessions and convert to a pandas series
benchmark = pd.read_csv(benchmark_file, header=None, names=['accession']).squeeze("columns")

# Validate GSE accessions
gse_pattern = re.compile(r'^GSE\d+$')
invalid_accessions = benchmark[~benchmark.str.match(gse_pattern)]

if not invalid_accessions.empty:
    raise ValueError(f"Invalid GSE accessions found: {', '.join(invalid_accessions)}")

print(f"Successfully loaded {len(benchmark)} benchmark dataset accessions")
# %% Extract sample metadata
import os
import time
import subprocess
import pandas as pd
from io import StringIO

####################
# METADATA FUNCTIONS
####################
def extract_sample_metadata(geo_accession):
    """
    Extract sample metadata for a given GEO accession.

    Parameters:
    geo_accession (str): GEO accession number (e.g., 'GSE122380')

    Returns:
    combined_df (pd.DataFrame): Combined metadata for all SRA IDs
    """
    # Define base directories
    base_data_dir = os.path.expanduser("~/Documents/UORCA/notebooks/2_extract_data")
    geo_data_dir = os.path.join(base_data_dir, f"{geo_accession}_data")

    # Define the get_files function
    def get_files(directory, suffix):
        """
        Recursively lists all files in a given directory and its subdirectories that end with the specified suffix,
        returning their absolute paths.

        Parameters:
        directory (str): The path to the directory to search in.
        suffix (str): The file suffix to look for (e.g., 'sra_ids.txt').

        Returns:
        list: A list of absolute file paths that match the given suffix.
        """
        matched_files = []

        try:
            # Expand user home directory and make path absolute
            abs_dir = os.path.abspath(os.path.expanduser(directory))

            if not os.path.exists(abs_dir):
                print(f"Directory '{abs_dir}' does not exist.")
                return []

            # Walk through directory and subdirectories
            for root, _, files in os.walk(abs_dir):
                for f in files:
                    if f.endswith(suffix):
                        full_path = os.path.join(root, f)
                        matched_files.append(full_path)

            return matched_files
        except FileNotFoundError:
            print(f"Directory '{directory}' not found.")
            return []
        except Exception as e:
            print(f"An error occurred while searching for files: {e}")
            return []

    # Use get_files to find 'sra_ids.txt' in geo_data_dir
    sra_ids_files = get_files(geo_data_dir, "sra_ids.txt")
    if not sra_ids_files:
        print(f"sra_ids.txt not found in directory {geo_data_dir}")
        return None
    elif len(sra_ids_files) > 1:
        print(f"Multiple sra_ids.txt files found in directory {geo_data_dir}. Using the first one: {sra_ids_files[0]}")
    sra_ids_file = sra_ids_files[0]

    # Temporary Debugging: Print the path to sra_ids.txt
    print(f"\n[DEBUG] Path to sra_ids.txt: {sra_ids_file}")

    # Temporary Debugging: Print the contents of sra_ids.txt
    print(f"\n[DEBUG] Contents of sra_ids.txt:")
    try:
        with open(sra_ids_file, 'r') as f:
            sra_contents = f.read()
            print(sra_contents)
    except Exception as e:
        print(f"Error reading {sra_ids_file}: {e}")
        return None

    # Set NCBI API key from environment variable
    entrez_api_key = os.getenv('ENTREZ_API_KEY')
    if entrez_api_key:
        print("\nAPI key detected.")
    else:
        print("\nNo API key detected; proceeding without it.")

    # Delay (in seconds) between requests to respect rate limits
    delay = 0.5 if not entrez_api_key else 0.1

    # List to store each SRA ID's fetched data as a dictionary
    data = []

    # Read the SRA IDs from the file
    try:
        with open(sra_ids_file, 'r') as ids_file:
            # Use pandas to read the file and skip the header
            sra_df = pd.read_csv(ids_file, sep='\t')
            sra_ids = sra_df['SRA_ID'].dropna().unique()
    except Exception as e:
        print(f"Error processing {sra_ids_file}: {e}")
        return None

    for sra_id in sra_ids:
        sra_id = sra_id.strip()
        if not sra_id:
            continue  # Skip empty lines or NaNs

        print(f"\nProcessing SRA ID: {sra_id}")

        # Construct the full command as a single shell command
        command = f"esearch -db sra -query {sra_id} | efetch -format runinfo"

        # Print the command being attempted for debugging
        print(f"Executing command: {command}")

        try:
            # Execute the combined command as a single shell pipeline
            result = subprocess.run(
                command, shell=True, capture_output=True, text=True, check=True
            )

            # Convert the result to a DataFrame and append it to the list
            csv_data = StringIO(result.stdout)
            df = pd.read_csv(csv_data)
            data.append(df)

        except subprocess.CalledProcessError as e:
            print(f"Error processing {sra_id}: {e}")
            print(f"Command output: {e.output}")
            continue  # Skip to the next SRA ID if there’s an error

        # Respect API rate limits
        time.sleep(delay)

    # Combine all DataFrames into one
    if data:
        combined_df = pd.concat(data, ignore_index=True)

        # Remove columns where all entries are NaN
        combined_df.dropna(axis=1, how='all', inplace=True)

        # Display the resulting DataFrame
        print("\nData fetching complete.")
        print(combined_df)
    else:
        print("No data was fetched.")
        combined_df = None

    return combined_df


def identify_and_remove_constants(df):
    """
    Identifies columns with constant values and removes them from the dataframe.

    Parameters:
    df (pd.DataFrame): Input dataframe to analyze

    Returns:
    tuple: (simplified_df, constant_values)
        - simplified_df: DataFrame with constant columns removed
        - constant_values: Dictionary of constant columns and their values
    """
    # Find columns where all values are the same
    constant_columns = {}
    for column in df.columns:
        unique_values = df[column].nunique()
        if unique_values == 1:
            constant_columns[column] = df[column].iloc[0]

    # Create new dataframe without constant columns
    simplified_df = df.drop(columns=constant_columns.keys())

    print(f"Removed {len(constant_columns)} constant columns")
    print("\nConstant values found:")
    for col, val in constant_columns.items():
        print(f"{col}: {val}")

    return simplified_df, constant_columns

# Save both versions for comparison
# simplified_metadata.to_csv("simplified_metadata.csv")

####################
# KALLISTO INPUT PREPARATION
####################
def prepare_kallisto_inputs(geo_accession: str) -> dict:
    """
    Prepare all required inputs for Kallisto quantification from a GEO accession.

    Parameters
    ----------
    geo_accession : str
        The GEO accession number to process

    Returns
    -------
    dict
        A dictionary containing:
        - study_summary: Text summary of the study
        - metadata: DataFrame of sample metadata with separate constant columns
        - fastq_files: List of available FASTQ files
        - index_files: List of available Kallisto indices
        - kallisto_docs: Kallisto documentation text
    """
    try:
        # Get study summary
        study_summary = get_study_summary(geo_accession)

        # Extract and process sample metadata
        metadata = extract_sample_metadata(geo_accession)
        if metadata is None:
            raise ValueError(f"Failed to extract metadata for {geo_accession}")

        # Separate constant columns
        simplified_metadata, constant_metadata = identify_and_remove_constants(metadata)

        # Get available FASTQ files
        data_dir = os.path.expanduser(f"~/Documents/UORCA/notebooks/2_extract_data/{geo_accession}_data")
        fastq_files = get_files(data_dir, ".fastq.gz")

        # Get available Kallisto indices
        index_dir = os.path.expanduser("~/Documents/UORCA/data/kallisto_indices/")
        index_files = get_files(index_dir, ".idx")

        # Get Kallisto documentation
        kallisto_docs = get_documentation("kallisto quant --help")

        return {
            "study_summary": study_summary,
            "metadata": simplified_metadata,
            "constant_metadata": constant_metadata,
            "fastq_files": fastq_files,
            "index_files": index_files,
            "kallisto_docs": kallisto_docs
        }

    except Exception as e:
        print(f"Error preparing Kallisto inputs: {str(e)}")
        raise
# %% Determine Kallisto parameters

client = OpenAI(
  api_key=openai_api_key,
)

####################
# KALLISTO PARAMETER MODELS
####################
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
    lr: bool = Field(..., description = "If reads are long read, e.g. from PacBio or Oxford Nanopore Technologies")
    justification: str = Field(..., description="Justification for each chosen parameter, including if the parameter was excluded")

class KallistoCommands(BaseModel):
    commands: List[KallistoCommand] = Field(description="List of Kallisto quantification commands for each sample")

def identify_kallisto_params(geo_accession: str):
    # Prepare all required inputs
    inputs = prepare_kallisto_inputs(geo_accession)

    prompt = f"""

## IDENTITY AND PURPOSE

You are an expert in bioinformatic analyses. You will be provided with various pieces of information, and use this information to determine the appropriate parameters for a Kallisto analysis.

## STEPS

1. Carefully digest the contents of the provided Kallisto documentation. Note that any existing knowledge you have of Kallisto may not be correct, so follow the documentation closely.
2. Carefully consider the contents of the sample metadata. Not all information will be relevant, however there will be content that will be needed.
- Note that the sample metadata is separated into a data frame, but there are some values which are constant over all values. These both comprise the sample metadata.
3. Carefully consider the dataset summary, which can contain information about experimental design etc.
4. When determining parameters corresponding to files, take care to ensure you only choose valid files (i.e. choose from options provided ONLY)
- That is: only choose one of the indicated FASTQ files, and one of the provided Kallisto index files.
5. Note that bootstrap samples should be set to exactly 5 for all samples.
6. Note that some parameters will be automatically set in certain conditions. If this is the case, then you do not need to manually set these parameters (and can leave these blank).
- However, also note that in some cases, parameters must be set. For example, for using single end reads, then standard deviation (sd) and fragment length MUST be supplied.
7. Repeat the above for ALL samples

## OUTPUT

Your output should consist of each parameter, and either:
- the value to be included for the parameter
- if the parameter should not be included, you should state NA
- For ALL chosen parameters, describe the justification for including the particular value, or excluding it.
- this should be done for all samples

This should be applied to all parameters identified as per the provided Kallisto documentation.

## INPUT

Kallisto documentation: {inputs['kallisto_docs']}
Kallisto documentation extra: Use the -lr parameter if using long reads
Dataset summary: {inputs['study_summary']}
FASTQ files: {inputs['fastq_files']}
Possible Kallisto indices: {inputs['index_files']}
Sample metadata:
- Note that the following values are constant for ALL samples, and should be considered as metadata:
    {inputs['constant_metadata']}
- The following is a representation of the remaining columns:
    {inputs['metadata'].to_string()}

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
# %% kallisto execution

####################
# KALLISTO EXECUTION FUNCTIONS
####################
def generate_output_name(fastq_path: str) -> str:
    """Generate meaningful output directory name from FASTQ path."""
    return Path(fastq_path).stem.split('_')[0]  # Gets SRR ID

from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import subprocess

def execute_single_kallisto(cmd: KallistoCommand, output_dir: Path) -> dict:
    """
    Execute a single Kallisto command and return results.

    Parameters
    ----------
    cmd : KallistoCommand
        The command to execute
    output_dir : Path
        Directory to write output to

    Returns
    -------
    dict
        Execution results containing status and any error messages
    """
    try:
        # Create the output directory
        output_dir.mkdir(parents=True, exist_ok=True)

        # Construct the Kallisto command string
        kallisto_cmd = (
            f"kallisto quant -i {cmd.index} -o {output_dir} "
            f"-t 4 --plaintext"
            # f"-t {max(1, N_CORES // 4)} --plaintext"
        )

        if cmd.bootstraps > 0:
            kallisto_cmd += f" --bootstrap-samples={cmd.bootstraps}"

        if cmd.single:
            kallisto_cmd += " --single"
            if cmd.frag_length:
                kallisto_cmd += f" -l {cmd.frag_length}"
            else:
                raise ValueError("frag_length (-l) must be set when using --single.")
            if cmd.sd:
                kallisto_cmd += f" -s {cmd.sd}"
            else:
                raise ValueError("sd (-s) must be set when using --single.")
        else:
            if cmd.fr_stranded:
                kallisto_cmd += " --fr-stranded"
            elif cmd.rf_stranded:
                kallisto_cmd += " --rf-stranded"

        # Append FASTQ files
        if cmd.fastq2 and cmd.fastq2.lower() != 'na':
            kallisto_cmd += f" {cmd.fastq1} {cmd.fastq2}"
        else:
            kallisto_cmd += f" {cmd.fastq1}"

        # Save command justification
        with open(output_dir / "justification.txt", "w") as f:
            f.write(cmd.justification)

        # Execute the command
        print(f"Executing command: {kallisto_cmd}")  # Added logging
        result = subprocess.run(
            kallisto_cmd,
            shell=True,
            check=True,
            capture_output=True,
            text=True
        )

        # Save command output
        with open(output_dir / "command_output.log", "w") as f:
            f.write(
                f"Command:\n{kallisto_cmd}\n\n"
                f"STDOUT:\n{result.stdout}\n\n"
                f"STDERR:\n{result.stderr}"
            )

        return {
            "status": "success",
            "output_dir": str(output_dir),
            "fastq": cmd.fastq1
        }

    except subprocess.CalledProcessError as e:
        # Capture subprocess-specific errors
        error_msg = (
            f"Subprocess error:\n"
            f"Return code: {e.returncode}\n"
            f"STDOUT: {e.stdout}\n"
            f"STDERR: {e.stderr}\n"
        )
        print(error_msg)

        # Save error to log file
        try:
            with open(output_dir / "error.log", "w") as f:
                f.write(error_msg)
        except Exception as log_error:
            print(f"Additionally failed to write error log: {log_error}")

        return {
            "status": "error",
            "error": error_msg,
            "output_dir": str(output_dir),
            "fastq": cmd.fastq1
        }

    except Exception as e:
        # Capture all other exceptions
        error_msg = f"General error: {str(e)}"
        print(error_msg)
        return {
            "status": "error",
            "error": error_msg,
            "output_dir": str(output_dir),
            "fastq": cmd.fastq1
        }

def execute_kallisto_commands(
    kallisto_commands: KallistoCommands,
    geo_accession: str,
    base_dir: Path | str,
) -> None:
    """
    Execute Kallisto quantification commands with proper directory structure.
    """
    # Convert base_dir to Path if it isn't already
    base_dir = Path(base_dir)

    # Create the base directory if it doesn't exist
    base_dir.mkdir(parents=True, exist_ok=True)

    # Create dataset-specific directory
    dataset_dir = base_dir / f"{geo_accession}_Kallisto"
    dataset_dir.mkdir(exist_ok=True)

    # Prepare arguments for parallel execution
    execution_args = []
    for cmd in kallisto_commands.commands:
        output_dir = dataset_dir / Path(cmd.output).name
        execution_args.append((cmd, output_dir))

    # Execute commands in parallel using threads
    with ThreadPoolExecutor(max_workers=max(1, N_CORES // 4)) as executor:
        # Submit all jobs
        future_to_cmd = {
            executor.submit(execute_single_kallisto, cmd, out_dir): (cmd, out_dir)
            for cmd, out_dir in execution_args
        }

        # Process results as they complete
        for future in tqdm(
            as_completed(future_to_cmd),
            total=len(future_to_cmd),
            desc=f"Processing {geo_accession}"
        ):
            cmd, out_dir = future_to_cmd[future]
            try:
                result = future.result()
                if result["status"] == "success":
                    print(f"Successfully processed {result['fastq']}")
                else:
                    print(f"Error processing {result['fastq']}: {result['error']}")
            except Exception as e:
                print(f"Error executing Kallisto command: {e}")
                raise


# %% Define function for processing
def process_benchmark_datasets(benchmark_accessions: pd.Series, base_dir: Path | str) -> None:
    ####################
    # MAIN EXECUTION
    ####################
    """
    Process all benchmark datasets through the Kallisto workflow.

    Parameters
    ----------
    benchmark_accessions : pd.Series
        Series containing GEO accessions to process
    base_dir : Path | str
        Base directory for output

    Returns
    -------
    None
    """
    # Convert base_dir to Path
    base_dir = Path(base_dir)

    # Create progress bar
    with tqdm(total=len(benchmark_accessions), desc="Processing datasets") as pbar:
        for accession in benchmark_accessions:
            try:
                print(f"\nProcessing {accession}...")

                # Get Kallisto parameters
                kallisto_params = identify_kallisto_params(accession)

                # Update output directories to use SRR IDs
                for cmd in kallisto_params.commands:
                    cmd.output = generate_output_name(cmd.fastq1)

                # Execute Kallisto commands
                execute_kallisto_commands(kallisto_params, accession, base_dir)

                print(f"Successfully processed {accession}")

            except Exception as e:
                print(f"Error processing {accession}: {str(e)}")
                continue

            finally:
                pbar.update(1)

# %% Main execution

if __name__ == "__main__":
    base_dir = Path("Benchmark_Kallisto")
    process_benchmark_datasets(benchmark, base_dir)
# %% Testing a single dataset

def test_single_dataset(geo_accession: str, base_dir: Path | str = "Kallisto_Test") -> None:
    """
    Test Kallisto quantification workflow on a single dataset.

    Parameters
    ----------
    geo_accession : str
        GEO accession number to process (e.g., 'GSE126096')
    base_dir : Path | str
        Base directory for output (default: "Kallisto_Test")

    Returns
    -------
    None
    """
    try:
        print(f"Processing {geo_accession}...")

        # Get Kallisto parameters
        kallisto_params = identify_kallisto_params(geo_accession)

        # Update output directories to use SRR IDs
        for cmd in kallisto_params.commands:
            cmd.output = generate_output_name(cmd.fastq1)

        # Execute Kallisto commands
        execute_kallisto_commands(kallisto_params, geo_accession, base_dir)

        print(f"Successfully processed {geo_accession}")

    except Exception as e:
        error_msg = f"""
Error processing {geo_accession}:
{str(e)}

Check the following locations for more details:
- {base_dir}/{geo_accession}_Kallisto/*/error.log (for Kallisto execution errors)
- {base_dir}/{geo_accession}_Kallisto/*/command_output.log (for successful runs)
- {base_dir}/{geo_accession}_Kallisto/*/justification.txt (for parameter choices)
"""
        print(error_msg)
        raise Exception(error_msg) from e

# Now to test with all small datasets
# List of GEO accession numbers with 20 or fewer samples
small_datasets = [
    "GSE119027",
    "GSE126096",
    "GSE133702",
    "GSE236761",
    "GSE242084",
    "GSE58702",
    "GSE63392",
    "GSE64712"
]

def process_small_datasets(datasets, base_dir="Kallisto_Test"):
    """
    Process multiple GEO datasets using the test_single_dataset function.

    Parameters
    ----------
    datasets : List[str]
        List of GEO accession numbers to process (e.g., 'GSE64712')
    base_dir : Path | str
        Base directory for output (default: "Kallisto_Test")

    Returns
    -------
    None
    """
    for geo_accession in datasets:
        print(f"\nStarting processing for {geo_accession}...")
        try:
            test_single_dataset(geo_accession, base_dir)
            print(f"Successfully processed {geo_accession}")
        except Exception as e:
            print(f"Failed to process {geo_accession}: {e}")
            continue
# Execute the processing for small datasets
process_small_datasets(small_datasets)

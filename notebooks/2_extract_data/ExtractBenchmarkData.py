
# %% Loading modules
import pandas as pd
from pathlib import Path
import re
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm
import time
from ratelimit import limits, sleep_and_retry

# Helper function to find scripts directory
def find_scripts_dir(start_path: Path) -> Path:
    """Find the directory containing the required scripts by searching up from start_path."""
    current = start_path
    required_scripts = ['process_geo.sh', 'download_metadata.R', 'download_FASTQs.sh']

    while current != current.parent:  # Stop at root directory
        if all((current / script).exists() for script in required_scripts):
            return current
        current = current.parent

    raise FileNotFoundError(f"Could not find directory containing required scripts: {required_scripts}")
# %% Define processing function

def process_geo_accession(accession: str,
                          n_spots: int = 20000,
                          output_dir: Path | None = None,
                          force: bool = True) -> subprocess.CompletedProcess:
    """
    Process a GEO accession using the process_geo.sh script.

    Parameters
    ----------
    accession : str
        The GEO accession number to process
    n_spots : int, optional
        Number of spots to download, by default 20000
    output_dir : Path, optional
        Output directory path. If None, uses accession_data
    force : bool, optional
        Whether to force overwrite existing files, by default True

    Returns
    -------
    subprocess.CompletedProcess
        The completed process object from running the script

    Raises
    ------
    subprocess.CalledProcessError
        If the script execution fails
    """
    # Validate accession format
    if not re.match(r'^GSE\d+$', accession):
        raise ValueError(f"Invalid GEO accession format: {accession}")

    # Set default output directory if none provided
    if output_dir is None:
        output_dir = Path(accession + "_data")

    # Construct script path - looking in the same directory as this notebook
        notebook_dir = Path().absolute()
        script_path = notebook_dir / "process_geo.sh"

        print(f"Notebook directory: {notebook_dir}")
        print(f"Looking for script at: {script_path}")
        print(f"Script exists: {script_path.exists()}")
        print(f"Script is executable: {os.access(str(script_path), os.X_OK) if script_path.exists() else False}")

        if not script_path.exists():
            # Try alternate locations relative to notebook directory
            alt_paths = [
                notebook_dir / "process_geo.sh",
                notebook_dir / ".." / "process_geo.sh",
                notebook_dir / "scripts" / "process_geo.sh"
            ]

            for path in alt_paths:
                if path.exists():
                    script_path = path
                    break
            else:
                searched_paths = [script_path] + alt_paths
                raise FileNotFoundError(
                    f"process_geo.sh script not found. Searched in:\n" +
                    "\n".join(f"- {p}" for p in searched_paths)
                )

    # Ensure script is executable
    if not os.access(str(script_path), os.X_OK):
        os.chmod(str(script_path), 0o755)
        print(f"Made script executable: {script_path}")

    # Build command with only short options
    cmd = [
        str(script_path),
        "-g", accession,
        "-o", str(output_dir),
        "-n", str(n_spots)
    ]

    if force:
        cmd.append("-f")  # Use short option for force

    print(f"Executing command: {' '.join(cmd)}")

    # Run the subprocess with explicit working directory
    try:
        result = subprocess.run(
            cmd,
            check=True,
            text=True,
            capture_output=True,
            cwd=script_path.parent  # Ensure correct working directory
        )
        print(f"Successfully processed {accession}")
        return result
    except subprocess.CalledProcessError as e:
        print(f"Error processing {accession}:")

        # Print error details from both log files
        master_log = output_dir / "master.log"
        processing_log = output_dir / "processing.log"

        if master_log.exists():
            print("\nContents of master.log:")
            print(master_log.read_text())

        if processing_log.exists():
            print("\nContents of processing.log:")
            print(processing_log.read_text())

        raise

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

# %% Example usage
# Process a single accession
result = process_geo_accession("GSE213001",
    n_spots=2000)

# %% All benchmark datasets
# Process all benchmark datasets
# for acc in benchmark:
#     try:
#         result = process_geo_accession(acc)
#     except Exception as e:
#         print(f"Failed to process {acc}: {str(e)}")

# Rate limit to 2 requests per second
@sleep_and_retry
@limits(calls=2, period=1)
def rate_limited_process(acc, n_spots=80000):
    try:
        return process_geo_accession(acc, n_spots=n_spots)
    except Exception as e:
        return f"Failed to process {acc}: {str(e)}"

# Process in parallel with progress bar
def process_all_accessions(accessions, max_workers=4):
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(rate_limited_process, acc) for acc in accessions]
        results = []
        for f in tqdm(futures, desc="Processing accessions"):
            results.append(f.result())
        return results

# Process all accessions
results = process_all_accessions(benchmark)

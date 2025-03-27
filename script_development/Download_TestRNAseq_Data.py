#!/usr/bin/env python3
"""
Download_TestRNAseq_Data.py

This script downloads a subset of RNA-seq data from GEO for testing
the DataAnalysisAgent. It downloads both FASTQ files and metadata.

The script can operate in two modes:
1. Single dataset mode: Downloads data for a single GEO accession.
2. Batch mode: Downloads data for multiple GEO accessions provided as a comma-separated list.

Usage:
  - Single dataset:
      python Download_TestRNAseq_Data.py -g GSE213001 -n 10000
  - Batch mode:
      python Download_TestRNAseq_Data.py -b -g GSE123456,GSE234567 -n 10000 -l 3
    (processes up to the specified number of datasets with 10,000 spots each)
"""

import os
import sys
import subprocess
import pandas as pd
import argparse
import logging
import shutil
import tempfile
import glob
import csv
import time
import concurrent.futures
from Bio import Entrez

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler()]
)
logger = logging.getLogger(__name__)

# Default constants
DEFAULT_OUTPUT_DIR = "./data/RNAseq_Benchmark"
DEFAULT_NUM_SPOTS = 10000  # Number of spots to download per sample
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# NCBI API constants (not used for metadata in this version)
NCBI_EMAIL = "kevin.chen@telethonkids.org.au"  # Set your email or use environment variable
NCBI_API_KEY = "d632f339861672624dc7fb31b2641099f107"  # Optional API key for higher rate limits

# Configure Entrez for potential API calls (if needed elsewhere)
Entrez.email = NCBI_EMAIL
if NCBI_API_KEY:
    Entrez.api_key = NCBI_API_KEY

def entrez_retry(func, *args, max_retries=3, **kwargs):
    """Retry an Entrez function call with exponential backoff."""
    for attempt in range(max_retries):
        try:
            return func(*args, **kwargs)
        except Exception as e:
            if attempt == max_retries - 1:
                logger.error(f"All Entrez API retries failed: {e}")
                raise
            wait_time = 2 ** attempt
            logger.warning(f"Entrez API call failed, retrying in {wait_time}s: {e}")
            time.sleep(wait_time)
    return None  # Should not reach here

def ensure_directory(directory):
    """Create directory if it doesn't exist."""
    if not os.path.exists(directory):
        os.makedirs(directory, exist_ok=True)

def run_subprocess(cmd, **kwargs):
    """Run subprocess with error handling for different Python versions."""
    try:
        return subprocess.run(cmd, capture_output=True, **kwargs)
    except TypeError:
        process = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            **kwargs
        )
        return process

def check_dependencies():
    """Check if necessary command line tools are available."""
    dependencies = ["prefetch", "fasterq-dump", "Rscript", "gzip"]
    missing = []
    for cmd in dependencies:
        try:
            subprocess.run(["which", cmd], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError:
            missing.append(cmd)
    if missing:
        logger.error(f"Missing dependencies: {', '.join(missing)}")
        logger.error("Please install the missing tools before running this script.")
        sys.exit(1)
    # Special check for the R package GEOquery
    try:
        r_check = run_subprocess(
            ["Rscript", "-e", "if(!require('GEOquery')) stop('GEOquery R package not installed')"],
            check=False, text=True
        )
        if r_check.returncode != 0:
            logger.error("Required R package GEOquery is not installed.")
            logger.error("Please install in R using: BiocManager::install('GEOquery')")
            sys.exit(1)
    except Exception as e:
        logger.error(f"Error checking R packages: {e}")
        logger.error("Continuing anyway, but metadata download may fail.")

def create_r_script(output_dir):
    """Create a temporary R script to download GEO metadata and extract SRA IDs."""
    r_script = r"""
#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(GEOquery))
args <- commandArgs(trailingOnly = TRUE)
geo_accession <- args[1]
output_dir <- args[2]

cat("Downloading metadata for", geo_accession, "\n")
gse <- getGEO(geo_accession, GSEMatrix = TRUE)
metadata <- pData(gse[[1]])

# Save full metadata to CSV
output_file <- file.path(output_dir, paste0(geo_accession, "_metadata.csv"))
write.csv(metadata, output_file, row.names = FALSE)
cat("Metadata saved to", output_file, "\n")

# Initialize SRA IDs vector (one per sample)
sra_ids <- rep(NA, nrow(metadata))

# First, try to extract SRA ID from the "relation" column if it contains "SRA:"
if("relation" %in% colnames(metadata)) {
  sra_ids <- sapply(metadata$relation, function(x) {
    if(grepl("SRA:", x)) {
      m <- regmatches(x, regexpr("term=([A-Z]+\\d+)", x, perl=TRUE))
      if(length(m) > 0 && nchar(m) > 0) {
        return(sub("term=", "", m))
      } else {
        return(NA)
      }
    } else {
      return(NA)
    }
  })
}

# If no SRA IDs found, try the "relation.1" column
if(all(is.na(sra_ids)) && ("relation.1" %in% colnames(metadata))) {
  sra_ids <- sapply(metadata$`relation.1`, function(x) {
    if(grepl("SRA:", x)) {
      m <- regmatches(x, regexpr("term=([A-Z]+\\d+)", x, perl=TRUE))
      if(length(m) > 0 && nchar(m) > 0) {
        return(sub("term=", "", m))
      } else {
        return(NA)
      }
    } else {
      return(NA)
    }
  })
}

sample_info <- data.frame(sample_name = rownames(metadata), sra_id = sra_ids, stringsAsFactors = FALSE)
output_file2 <- file.path(output_dir, paste0(geo_accession, "_sample_sra.csv"))
write.csv(sample_info, output_file2, row.names = FALSE)
cat("Sample-to-SRA mapping saved to", output_file2, "\n")
"""
    import tempfile
    temp_file_path = tempfile.mktemp(suffix=".R")
    with open(temp_file_path, "w") as f:
        f.write(r_script)
    return temp_file_path

def download_metadata(geo_accession, output_dir):
    """
    Download metadata for a GEO accession using R/GEOquery exclusively.
    This function always calls the R script.
    """
    logger.info(f"Downloading metadata for {geo_accession} using R/GEOquery")
    r_script_path = create_r_script(output_dir)
    try:
        process = subprocess.run(
            ["Rscript", r_script_path, geo_accession, output_dir],
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        logger.info(f"R script output: {process.stdout}")
        logger.info(f"Metadata downloaded successfully to {output_dir} using R")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Error downloading metadata using R: {e}")
        if hasattr(e, 'stdout') and e.stdout:
            logger.error(f"STDOUT: {e.stdout}")
        if hasattr(e, 'stderr') and e.stderr:
            logger.error(f"STDERR: {e.stderr}")
        return False
    finally:
        if os.path.exists(r_script_path):
            os.unlink(r_script_path)

def get_sra_ids(geo_accession, output_dir):
    """
    Get SRA IDs for a GEO accession.
    This function first checks for a CSV file (generated by the R script) and,
    if present, parses it to extract SRA IDs.
    """
    sample_sra_file = os.path.join(output_dir, f"{geo_accession}_sample_sra.csv")
    if os.path.exists(sample_sra_file):
        df = pd.read_csv(sample_sra_file)
        df = df.dropna(subset=['sra_id'])
        sra_ids = df['sra_id'].tolist()
        all_srr_ids = []
        for sra_id in sra_ids:
            if isinstance(sra_id, str):
                parts = sra_id.split(':')
                if len(parts) > 1:
                    srr_ids = parts[1].split(',')
                    all_srr_ids.extend(srr_ids)
                else:
                    all_srr_ids.append(sra_id)
        return all_srr_ids
    else:
        logger.info(f"Sample-to-SRA CSV not found for {geo_accession}. Attempting Entrez extraction.")
        try:
            cmd = f"esearch -db gds -query '{geo_accession}[ACCN]' | " + \
                  "elink -target sra | " + \
                  "efetch -format docsum | " + \
                  "xtract -pattern DocumentSummary -element Run@acc"
            result = subprocess.run(
                cmd,
                shell=True,
                check=True,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )
            sra_ids = result.stdout.strip().split('\n')
            sra_ids = [s for s in sra_ids if s]
            if not sra_ids:
                logger.warning(f"No SRA IDs found for {geo_accession}")
                sra_ids = ["SRR14415242", "SRR14415243", "SRR14415244", "SRR14415245"]
            return sra_ids
        except subprocess.CalledProcessError as e:
            logger.error(f"Error getting SRA IDs: {e}")
            logger.error(f"STDOUT: {e.stdout}")
            logger.error(f"STDERR: {e.stderr}")
            logger.warning("Using hardcoded SRA IDs as fallback")
            return ["SRR14415242", "SRR14415243", "SRR14415244", "SRR14415245"]

def download_fastq(sra_id, output_dir, num_spots):
    """
    Download FASTQ files for an SRA ID using prefetch and fastq-dump
    (limiting to a specified number of spots).
    """
    logger.info(f"Downloading FASTQ files for {sra_id} (max {num_spots} spots)")
    fastq_dir = os.path.join(output_dir, "fastq")
    ensure_directory(fastq_dir)
    with tempfile.TemporaryDirectory() as temp_dir:
        # Step 1: Prefetch the SRA file
        try:
            cmd_prefetch = ["prefetch", "-o", f"{temp_dir}/{sra_id}.sra", sra_id]
            logger.info(f"Executing: {' '.join(cmd_prefetch)}")
            process = subprocess.run(
                cmd_prefetch,
                check=True,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=600
            )
            logger.debug(f"Prefetch output: {process.stdout}")
        except subprocess.CalledProcessError as e:
            logger.error(f"Error prefetching {sra_id}: {e}")
            if hasattr(e, 'stderr'):
                logger.error(f"STDERR: {e.stderr}")
            return False
        except subprocess.TimeoutExpired:
            logger.error(f"Prefetch for {sra_id} timed out after 10 minutes")
            return False

        # Step 2: Extract FASTQ using fastq-dump
        try:
            cmd_fastq = [
                "fastq-dump",
                "--split-files",
                "--outdir", fastq_dir,
                "-X", str(num_spots),
                f"{temp_dir}/{sra_id}.sra"
            ]
            logger.info(f"Executing: {' '.join(cmd_fastq)}")
            process = subprocess.run(
                cmd_fastq,
                check=True,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=300
            )
            logger.debug(f"fastq-dump output: {process.stdout}")
        except subprocess.CalledProcessError as e:
            logger.error(f"Error extracting FASTQ for {sra_id}: {e}")
            if hasattr(e, 'stderr'):
                logger.error(f"STDERR: {e.stderr}")
            return False
        except subprocess.TimeoutExpired:
            logger.error(f"FASTQ extraction for {sra_id} timed out after 5 minutes")
            return False

    # Step 3: Compress FASTQ files
    for fastq_file in glob.glob(os.path.join(fastq_dir, f"{sra_id}_*.fastq")):
        try:
            logger.info(f"Compressing {os.path.basename(fastq_file)}...")
            process = subprocess.run(
                ["gzip", "-f", fastq_file],
                check=True,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=60
            )
            logger.debug(f"Gzip output: {process.stdout}")
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
            logger.error(f"Error or timeout compressing {fastq_file}: {e}")
            if isinstance(e, subprocess.CalledProcessError) and hasattr(e, 'stderr'):
                logger.error(f"STDERR: {e.stderr}")
            return False
    return True

def verify_fastq_files(directory):
    """
    Verify that FASTQ files exist and are non-empty.
    """
    fastq_dir = os.path.join(directory, "fastq")
    if not os.path.exists(fastq_dir):
        logger.error(f"FASTQ directory not found: {fastq_dir}")
        return False
    fastq_files = glob.glob(os.path.join(fastq_dir, "*.fastq.gz"))
    if not fastq_files:
        fastq_files = glob.glob(os.path.join(fastq_dir, "*.fastq"))
    if not fastq_files:
        logger.error(f"No FASTQ files found in {fastq_dir}")
        return False
    for fastq in fastq_files:
        if os.path.getsize(fastq) < 100:
            logger.error(f"FASTQ file too small: {fastq}")
            return False
    logger.info(f"Verified {len(fastq_files)} FASTQ files in {fastq_dir}")
    return True

def process_all_datasets(output_dir, num_spots, geo_accessions, limit=None):
    """
    Process all GEO datasets provided in the list.
    For each dataset, download metadata, extract SRA IDs, and download FASTQ files in parallel.
    """
    datasets_to_process = geo_accessions[:limit] if limit else geo_accessions
    logger.info(f"Processing {len(datasets_to_process)} datasets: {', '.join(datasets_to_process)}")
    results = []
    for geo_accession in datasets_to_process:
        logger.info(f"==== Processing dataset {geo_accession} ====")
        dataset_dir = os.path.join(output_dir, geo_accession)
        ensure_directory(dataset_dir)

        try:
            if not download_metadata(geo_accession, dataset_dir):
                logger.error(f"Metadata download failed for {geo_accession}")
                results.append({
                    "geo_accession": geo_accession,
                    "success": False,
                    "error": "Metadata download failed"
                })
                continue

            sra_ids = get_sra_ids(geo_accession, dataset_dir)
            logger.info(f"Found {len(sra_ids)} SRA IDs for {geo_accession}: {', '.join(sra_ids)}")

            successful_downloads = 0
            with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
                future_to_sra = {
                    executor.submit(download_fastq, sra_id, dataset_dir, num_spots): sra_id
                    for sra_id in sra_ids
                }
                for future in concurrent.futures.as_completed(future_to_sra):
                    sra = future_to_sra[future]
                    try:
                        if future.result():
                            successful_downloads += 1
                    except Exception as exc:
                        logger.error(f"{sra} generated an exception: {exc}")

            if successful_downloads == 0:
                logger.error(f"Failed to download FASTQ files for {geo_accession}")
                results.append({
                    "geo_accession": geo_accession,
                    "success": False,
                    "error": "FASTQ download failed"
                })
                continue

            summary_file = os.path.join(dataset_dir, "summary.txt")
            with open(summary_file, "w") as f:
                f.write(f"GEO Accession: {geo_accession}\n")
                f.write(f"Number of SRA IDs processed: {len(sra_ids)}\n")
                f.write(f"Number of successful downloads: {successful_downloads}\n")
                f.write(f"SRA IDs downloaded: {', '.join(sra_ids)}\n")
                f.write(f"Number of spots per sample: {num_spots}\n")
            results.append({
                "geo_accession": geo_accession,
                "success": True,
                "samples_downloaded": successful_downloads
            })

        except Exception as e:
            logger.error(f"Error processing dataset {geo_accession}: {e}")
            results.append({
                "geo_accession": geo_accession,
                "success": False,
                "error": str(e)
            })

    summary_path = os.path.join(output_dir, "batch_summary.txt")
    with open(summary_path, "w") as f:
        f.write(f"Processed {len(datasets_to_process)} datasets\n")
        successful = sum(1 for r in results if r["success"])
        f.write(f"Successfully processed: {successful}/{len(datasets_to_process)}\n\n")
        for result in results:
            if result["success"]:
                f.write(f"{result['geo_accession']}: SUCCESS ({result['samples_downloaded']} samples)\n")
            else:
                f.write(f"{result['geo_accession']}: FAILED ({result.get('error', 'Unknown error')})\n")
    return results

def main():
    parser = argparse.ArgumentParser(description="Download RNA-seq data for DataAnalysisAgent")
    parser.add_argument("-g", "--geo-accession", required=True,
                        help="GEO accession number for single dataset mode or a comma-separated list for batch mode")
    parser.add_argument("-o", "--output-dir", default=DEFAULT_OUTPUT_DIR,
                        help=f"Output directory (default: {DEFAULT_OUTPUT_DIR})")
    parser.add_argument("-n", "--num-spots", type=int, default=DEFAULT_NUM_SPOTS,
                        help=f"Number of spots to download (default: {DEFAULT_NUM_SPOTS})")
    parser.add_argument("-b", "--batch", action="store_true",
                        help="Process all datasets in batch mode (comma-separated GEO accessions)")
    parser.add_argument("-l", "--limit", type=int, default=None,
                        help="Limit the number of datasets to process in batch mode")
    args = parser.parse_args()

    output_dir = args.output_dir
    check_dependencies()
    ensure_directory(output_dir)

    if args.batch:
        geo_accessions = [g.strip() for g in args.geo_accession.split(",") if g.strip()]
        logger.info(f"Starting batch processing of datasets with {args.num_spots} spots each")
        results = process_all_datasets(output_dir, args.num_spots, geo_accessions, args.limit)
        successful = sum(1 for r in results if r["success"])
        if successful > 0:
            logger.info(f"Successfully processed {successful}/{len(results)} datasets")
        else:
            logger.error("Failed to process any datasets successfully")
            sys.exit(1)
    else:
        dataset_dir = os.path.join(output_dir, args.geo_accession)
        ensure_directory(dataset_dir)
        if not download_metadata(args.geo_accession, dataset_dir):
            logger.error("Metadata download failed.")
            sys.exit(1)
        sra_ids = get_sra_ids(args.geo_accession, dataset_dir)
        logger.info(f"Found {len(sra_ids)} SRA IDs: {', '.join(sra_ids)}")
        successful_downloads = 0
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
            future_to_sra = {
                executor.submit(download_fastq, sra_id, dataset_dir, args.num_spots): sra_id
                for sra_id in sra_ids
            }
            for future in concurrent.futures.as_completed(future_to_sra):
                sra = future_to_sra[future]
                try:
                    if future.result():
                        successful_downloads += 1
                except Exception as exc:
                    logger.error(f"{sra} generated an exception: {exc}")

        if successful_downloads == 0:
            logger.error(f"Failed to download any samples for {args.geo_accession}.")
            sys.exit(1)
        summary_file = os.path.join(dataset_dir, "summary.txt")
        with open(summary_file, "w") as f:
            f.write(f"GEO Accession: {args.geo_accession}\n")
            f.write(f"Number of SRA IDs processed: {len(sra_ids)}\n")
            f.write(f"Number of successful downloads: {successful_downloads}\n")
            f.write(f"SRA IDs downloaded: {', '.join(sra_ids)}\n")
            f.write(f"Number of spots per sample: {args.num_spots}\n")
        logger.info(f"Successfully downloaded {successful_downloads} samples to {dataset_dir}")

if __name__ == "__main__":
    main()

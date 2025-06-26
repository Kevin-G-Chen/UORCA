#!/usr/bin/env python3
"""
Optimized Dataset Identification Script
======================================

This script performs dataset identification with optimal workflow and parallelization.
Key features:
- API key-based dynamic rate limiting (9 req/sec with key, 3 req/sec without)
- Parallel processing for GEO summary retrieval, SRA fetching, and dataset validation
- Optimal workflow: get all info → validate → cluster/score only valid datasets
- Enhanced reporting with validation reasons

Workflow (optimal efficiency):
1. Parallel dataset info fetch using esummary v2.0 (title, summary, BioProject, etc.)
2. Fetch SRA metadata for ALL datasets using BioProject directly
3. Validate ALL datasets based on complete SRA criteria
4. Only perform expensive operations (clustering, relevance) on valid datasets
5. Include ALL datasets in final output for transparency

Usage:
    python OptimizedDatasetIdentification.py "research query" --max-datasets 10
"""

import argparse
import asyncio
import io
import json
import logging
import os
import statistics
import sys
import threading
import time
import warnings
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
from functools import partial

import numpy as np
import pandas as pd
import requests
from Bio import Entrez
from dotenv import load_dotenv
from openai import OpenAI
from pydantic import BaseModel
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from tqdm import tqdm

# Load environment variables
load_dotenv()

# Configure Entrez
Entrez.email = os.getenv("ENTREZ_EMAIL", "user@example.com")
if os.getenv("ENTREZ_API_KEY"):
    Entrez.api_key = os.getenv("ENTREZ_API_KEY")

# OpenAI client
client = OpenAI(api_key=os.getenv("OPENAI_API_KEY"))

# Rate limiting for API calls
class APIRateLimiter:
    """Thread-safe rate limiter for API calls with dynamic rates based on API key."""
    def __init__(self, base_delay: float = None):
        # Auto-detect optimal delay based on API key presence
        if base_delay is None:
            if os.getenv("ENTREZ_API_KEY"):
                # With API key: up to 10 requests/second = 0.1 second delay
                self.min_delay = 0.11  # Slightly conservative
                logging.info("🚀 API key detected - using fast rate limit (9 requests/sec)")
            else:
                # Without API key: up to 3 requests/second = 0.33 second delay
                self.min_delay = 0.34  # Slightly conservative
                logging.info("⚠️ No API key - using standard rate limit (3 requests/sec)")
        else:
            self.min_delay = base_delay

        self.last_call_time = 0
        self.lock = threading.Lock()

    def wait(self):
        """Wait if necessary to respect rate limits."""
        with self.lock:
            current_time = time.time()
            time_since_last_call = current_time - self.last_call_time
            if time_since_last_call < self.min_delay:
                sleep_time = self.min_delay - time_since_last_call
                time.sleep(sleep_time)
            self.last_call_time = time.time()

# Global rate limiter instance with auto-detection
api_rate_limiter = APIRateLimiter()

# Constants
MAX_RETRIES = 3
EMBEDDING_MODEL = "text-embedding-3-small"

# Configure logging
def setup_logging(verbose: bool = False) -> None:
    """Setup logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    format_str = "%(asctime)s - %(levelname)s - %(message)s"

    # Create logs directory if it doesn't exist
    logs_dir = Path("logs")
    logs_dir.mkdir(exist_ok=True)

    # Clear any existing handlers to avoid conflicts
    root_logger = logging.getLogger()
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)

    # Configure logging to both file and console with forced configuration
    file_handler = logging.FileHandler(logs_dir / "dataset_identification.log")
    console_handler = logging.StreamHandler(sys.stdout)

    logging.basicConfig(
        level=level,
        format=format_str,
        handlers=[file_handler, console_handler],
        force=True  # Force reconfiguration even if logging was already configured
    )

    # Also configure warnings to use the same format
    logging.captureWarnings(True)
    warnings_logger = logging.getLogger('py.warnings')
    warnings_logger.setLevel(level)

    # Debug: Print configuration info to stderr to ensure it's visible
    print(f"🔧 LOGGING SETUP DEBUG:", file=sys.stderr)
    print(f"   - Level: {logging.getLevelName(level)}", file=sys.stderr)
    print(f"   - Format: {format_str}", file=sys.stderr)
    print(f"   - Handlers count: {len(logging.getLogger().handlers)}", file=sys.stderr)
    print(f"   - File handler: {logs_dir / 'dataset_identification.log'}", file=sys.stderr)
    print(f"   - Console handler: stdout", file=sys.stderr)
    print("🔧 LOGGING SETUP COMPLETE", file=sys.stderr)

def load_prompt(file_path: str) -> str:
    """Load prompt from file."""
    return Path(file_path).read_text().strip()

# Pydantic models for structured output
class ExtractedTerms(BaseModel):
    extracted_terms: List[str]
    expanded_terms: List[str]

class Assessment(BaseModel):
    ID: str
    RelevanceScore: int
    Justification: str

class Assessments(BaseModel):
    assessments: List[Assessment]

def call_openai_json(prompt: str, user_input: str, response_format: BaseModel) -> Any:
    """Make OpenAI API call with JSON response format."""
    response = client.beta.chat.completions.parse(
        model="gpt-4o-mini",
        messages=[
            {"role": "system", "content": prompt},
            {"role": "user", "content": user_input}
        ],
        response_format=response_format
    )
    return response.choices[0].message.parsed

def extract_terms(research_query: str) -> ExtractedTerms:
    """Extract search terms from research query."""
    prompt = load_prompt("./main_workflow/prompts/dataset_identification/extract_terms.txt")
    return call_openai_json(prompt, research_query, ExtractedTerms)

def perform_search(term: str, max_results: int = 2000) -> List[str]:
    """Search GEO database for a single term with RNA-seq filter."""
    search = f"{term} AND (\"Expression profiling by high throughput sequencing\"[Filter])"
    handle = Entrez.esearch(
        db="gds",
        term=search,
        retmode="xml",
        retmax=max_results
    )
    search_results = Entrez.read(handle)
    handle.close()

    geo_ids = search_results.get("IdList", [])
    logging.info(f"Found {len(geo_ids)} results for term: {term}")

    return geo_ids

def get_basic_dataset_info(
    geo_ids: List[str],
    api_delay: float = 0.4,
    workers: Optional[int] = None,
) -> pd.DataFrame:
    """Get complete dataset information in parallel using esummary v2.0."""

    unique_ids = list(dict.fromkeys(geo_ids))
    if workers is None:
        workers = 6 if os.getenv("ENTREZ_API_KEY") else 2

    logging.info(
        f"Fetching complete dataset information for {len(unique_ids)} unique datasets using {workers} workers..."
    )

    datasets = []

    def fetch_info(geo_id: str) -> Optional[Dict[str, Any]]:
        """Fetch esummary information for a single GEO ID."""
        try:
            # Get UID for this GEO ID
            api_rate_limiter.wait()
            handle = Entrez.esearch(db="gds", term=geo_id)
            search_result = Entrez.read(handle)
            handle.close()
            uids = search_result.get("IdList", [])

            if not uids:
                logging.warning(f"No UID found for {geo_id}")
                return None

            uid = uids[0]

            api_rate_limiter.wait()
            handle = Entrez.esummary(db="gds", id=uid, version="2.0", retmode="xml")
            doc = Entrez.read(handle)["DocumentSummarySet"]["DocumentSummary"][0]
            handle.close()

            accession = doc["Accession"]
            if accession and accession.startswith("GSE"):
                logging.info(f"Fetched dataset {accession} with UID {uid}")
                return {
                    "ID": accession,
                    "Title": doc["title"],
                    "Summary": doc["summary"],
                    "Accession": accession,
                    "BioProject": doc["BioProject"],
                    "Species": doc["taxon"],
                    "Date": doc["PDAT"],
                    "NumSamples": doc["n_samples"],
                    "PrimaryPubMedID": (
                        doc["PubMedIds"][0] if doc.get("PubMedIds") else None
                    ),
                }
            else:
                logging.warning(f"Invalid accession for {geo_id}")
                return None

        except Exception as e:
            logging.warning(f"Error processing GEO ID {geo_id}: {e}")
            return None

    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(fetch_info, gid): gid for gid in unique_ids}
        for future in as_completed(futures):
            result = future.result()
            if result:
                datasets.append(result)

    df = pd.DataFrame(datasets)
    if not df.empty:
        df = df.drop_duplicates(subset=["ID"]).reset_index(drop=True)

    logging.info(f"Found {len(df)} unique GSE datasets with complete information")
    return df

# Streamlined SRA processing functions (old XML parsing functions removed)
def fetch_runinfo_from_bioproject(bioproject: str, api_delay: float = 0.4) -> pd.DataFrame:
    """Fetch RunInfo data using streamlined BioProject approach."""
    try:
        # 1. Search SRA and keep the server-side history
        srch = Entrez.read(Entrez.esearch(db="sra",
                                        term=f"{bioproject}[BioProject]",
                                        usehistory="y"))
        count, webenv, qk = int(srch["Count"]), srch["WebEnv"], srch["QueryKey"]

        if count == 0:
            logging.warning(f"No SRA records found for BioProject {bioproject}")
            return pd.DataFrame()

        time.sleep(api_delay)

        # 2. Fetch run-level metadata
        handle = Entrez.efetch(db="sra",
                               rettype="runinfo",
                               retmode="text",
                               WebEnv=webenv,
                               query_key=qk,
                               retmax=count)

        # 3. Convert bytes → str → DataFrame
        csv_bytes = handle.read()                # bytes
        csv_text = csv_bytes.decode('utf-8')    # str
        runs = pd.read_csv(io.StringIO(csv_text))
        handle.close()

        logging.info(f"Fetched {len(runs)} runs for BioProject {bioproject}")
        return runs

    except Exception as e:
        logging.warning(f"Error fetching RunInfo for BioProject {bioproject}: {e}")
        return pd.DataFrame()


# ==================== OPTIONAL: DATASET SIZE CALCULATION FUNCTION ====================
# The following function calculates dataset sizes from SRA runinfo data
# Remove this function if you don't need dataset size information
def calculate_dataset_sizes_from_runinfo(sra_df: pd.DataFrame) -> Dict[str, int]:
    """
    Calculate dataset sizes using the size_MB column already present in runinfo data.

    This is much faster than calling vdb-dump for each SRR individually.

    Returns:
        Dictionary mapping GEO_Accession to total size in bytes
    """
    if sra_df.empty or 'size_MB' not in sra_df.columns or 'GEO_Accession' not in sra_df.columns:
        logging.warning("Cannot calculate dataset sizes: missing required columns in runinfo data")
        return {}

    # Convert size_MB to bytes and group by GEO accession
    sra_df_copy = sra_df.copy()
    sra_df_copy['size_MB'] = pd.to_numeric(sra_df_copy['size_MB'], errors='coerce').fillna(0)
    sra_df_copy['size_bytes'] = (sra_df_copy['size_MB'] * 1024 * 1024).astype(int)

    # Group by GEO accession and sum sizes
    dataset_sizes = sra_df_copy.groupby('GEO_Accession')['size_bytes'].sum().to_dict()

    # Log results
    total_datasets = len(dataset_sizes)
    total_size_gb = sum(dataset_sizes.values()) / (1024**3)
    logging.info(f"Dataset size calculation complete: {total_datasets} datasets, total: {total_size_gb:.2f} GB")

    for geo_acc, size_bytes in dataset_sizes.items():
        logging.info(f"Dataset {geo_acc}: {size_bytes:,} bytes ({size_bytes/(1024**3):.2f} GB)")

    return dataset_sizes
# =====================================================================================


# Removed functions for streamlined workflow:
# - validate_and_size_datasets: validation now happens after SRA merge
# - fetch_geo_summaries_for_valid: all info now fetched via esummary v2.0
# - get_geo_uid, get_sra_uid, fetch_sra_xml, parse_bioprojects: replaced with direct BioProject approach

def get_embedding(text: str) -> List[float]:
    """Get embedding for text using OpenAI API."""
    try:
        response = client.embeddings.create(
            model=EMBEDDING_MODEL,
            input=text
        )
        return response.data[0].embedding
    except Exception as e:
        logging.warning(f"Error getting embedding: {e}")
        return [0.0] * 1536  # Default embedding size

def embed_datasets(datasets_df: pd.DataFrame) -> pd.DataFrame:
    """Generate embeddings for datasets."""
    logging.info("Generating embeddings for valid datasets...")

    embeddings = []
    for _, row in datasets_df.iterrows():
        # Combine title and summary for embedding
        text = f"{row.get('Title', '')} {row.get('Summary', '')}"
        embedding = get_embedding(text)
        embeddings.append(embedding)

        time.sleep(0.1)  # Small delay for API rate limiting

    # Add embeddings to dataframe
    datasets_df = datasets_df.copy()
    datasets_df['embedding'] = embeddings

    return datasets_df

def cluster_datasets(datasets_df: pd.DataFrame, n_clusters: int = None) -> pd.DataFrame:
    """Cluster datasets based on embeddings."""
    if len(datasets_df) <= 1:
        datasets_df = datasets_df.copy()
        datasets_df['Cluster'] = 0
        return datasets_df

    if n_clusters is None:
        n_clusters = min(len(datasets_df) // 3, 10)  # Heuristic for cluster count
        n_clusters = max(n_clusters, 1)

    logging.info(f"Clustering {len(datasets_df)} datasets into {n_clusters} clusters...")

    # Prepare embeddings
    embeddings = np.array(datasets_df['embedding'].tolist())

    # Standardize embeddings
    scaler = StandardScaler()
    embeddings_scaled = scaler.fit_transform(embeddings)

    # Perform clustering
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    clusters = kmeans.fit_predict(embeddings_scaled)

    # Add cluster labels
    datasets_df = datasets_df.copy()
    datasets_df['Cluster'] = clusters

    return datasets_df

def select_representative_datasets(datasets_df: pd.DataFrame, max_datasets: int) -> pd.DataFrame:
    """Select representative datasets from clusters."""
    if len(datasets_df) <= max_datasets:
        return datasets_df

    logging.info(f"Selecting {max_datasets} representative datasets from clusters...")

    representatives = []

    # Get representatives from each cluster
    for cluster_id in datasets_df['Cluster'].unique():
        cluster_datasets = datasets_df[datasets_df['Cluster'] == cluster_id]

        # For now, just take the first dataset from each cluster
        # Could implement more sophisticated selection logic here
        representatives.append(cluster_datasets.iloc[0])

    representatives_df = pd.DataFrame(representatives)

    # If we still have too many, take the top ones by some criteria
    if len(representatives_df) > max_datasets:
        # For now, just take the first max_datasets
        representatives_df = representatives_df.head(max_datasets)

    return representatives_df.reset_index(drop=True)

# ==================== OPTIONAL: REPEATED ASSESSMENT FUNCTIONS ====================
# The following functions implement repeated assessment with averaging
# Remove both assess_subbatch and repeated_relevance if you want single assessment only
async def assess_subbatch(df: pd.DataFrame, query: str, schema, key: str, rep: int, idx: int, total_batches: int, sem: asyncio.Semaphore) -> List[Assessment]:
    """Assess a sub-batch of datasets."""
    async with sem:
        try:
            # Prepare data for assessment
            assessment_data = []
            for _, row in df.iterrows():
                assessment_data.append({
                    "ID": row['ID'],
                    "Species": row.get('Species', 'Unknown'),
                    "Tissue": "Unknown",
                    "Technique": row.get('Technique', 'RNA-seq'),
                    "Summary": row.get('Summary', '')
                })

            prompt = load_prompt("./main_workflow/prompts/dataset_identification/assess_relevance.txt")
            user_input = f"Research Query: {query}\n\nDatasets:\n{json.dumps(assessment_data, indent=2)}"

            assessments = call_openai_json(prompt, user_input, Assessments)
            logging.info(f"Completed assessment batch {idx+1}/{total_batches} (repetition {rep+1})")
            return assessments.assessments
        except Exception as e:
            logging.warning(f"Error in assessment batch {idx+1}/{total_batches} (rep {rep+1}): {e}")
            # Return default assessments
            return [Assessment(ID=row['ID'], RelevanceScore=5, Justification="Assessment failed") for _, row in df.iterrows()]

async def repeated_relevance(df: pd.DataFrame, query: str, repeats: int = 3, batch_size: int = 10, openai_api_jobs: int = 3) -> pd.DataFrame:
    """Perform repeated relevance assessment with averaging."""
    print(f"📊 Starting relevance scoring: {repeats} repetitions, batch size {batch_size}, parallel API jobs: {openai_api_jobs}")
    sem = asyncio.Semaphore(openai_api_jobs)
    tasks = []

    # create sub-batches
    batches = [df.iloc[i:i+batch_size] for i in range(0, len(df), batch_size)]
    total_batches = len(batches)

    for rep in range(repeats):
        for idx, sub in enumerate(batches):
            tasks.append(assess_subbatch(sub, query, None, "assessments", rep, idx, total_batches, sem))

    all_results = await asyncio.gather(*tasks)

    # flatten and aggregate
    coll: Dict[str, Dict[str, Any]] = {}
    for result in all_results:
        for a in result:
            entry = coll.setdefault(a.ID, {'scores': [], 'justifications': []})
            entry['scores'].append(a.RelevanceScore)
            entry['justifications'].append(a.Justification)

    records: List[Dict[str, Any]] = []
    for id_, v in coll.items():
        rec = {'ID': id_, 'RelevanceScore': round(statistics.mean(v['scores']), 2)}
        for i, (score, just) in enumerate(zip(v['scores'], v['justifications']), 1):
            rec[f'Run{i}Score'] = score
            rec[f'Run{i}Justification'] = just
        records.append(rec)

    return pd.DataFrame(records)
# ===============================================================================

def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description='Identify and evaluate relevant RNA-seq datasets from GEO for a biological research query',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # === Core Parameters ===
    core_group = parser.add_argument_group('Core Options', 'Essential parameters for dataset identification')
    core_group.add_argument('-q', '--query', required=True,
                           help='Biological research query (e.g., "neuroblastoma tumor vs normal")')
    core_group.add_argument('-o', '--output', default='./dataset_identification_results',
                           help='Output directory for results')
    core_group.add_argument('-t', '--threshold', type=float, default=7.0,
                           help='Relevance score threshold (0-10) for including datasets')

    # === Search Parameters ===
    search_group = parser.add_argument_group('Search Options', 'Control dataset search and evaluation')
    search_group.add_argument('--max-datasets', type=int, default=2000,
                             help='Maximum datasets to retrieve per search term')
    search_group.add_argument('--max-evaluate', type=int, default=None,
                             help='Limit number of datasets to evaluate (useful for testing)')
    search_group.add_argument('--use-clustering', action='store_true',
                             help='Use AI clustering to select diverse, representative datasets')

    # === Output Options ===
    output_group = parser.add_argument_group('Output Options', 'Control output format and content')
    output_group.add_argument('--multi-dataset-csv', action='store_true',
                             help='Generate CSV formatted for multi-dataset batch analysis')
    output_group.add_argument('--save-clustering-details', action='store_true',
                             help='Save detailed clustering information to separate files')

    # === Advanced Parameters ===
    advanced_group = parser.add_argument_group('Advanced Options', 'Fine-tune algorithm behavior (expert users)')
    advanced_group.add_argument('--scoring-rounds', type=int, default=3,
                               help='Number of independent relevance scoring rounds for reliability')
    advanced_group.add_argument('--batch-size', type=int, default=20,
                               help='Datasets per AI evaluation batch (affects memory usage)')
    advanced_group.add_argument('--parallel-ai-jobs', type=int, default=4,
                               help='Concurrent AI evaluation jobs (limited by API rate limits)')
    advanced_group.add_argument('--parallel-embedding-jobs', type=int, default=4,
                               help='Concurrent embedding generation jobs for clustering')
    advanced_group.add_argument('--min-cluster-size', type=int, default=2,
                               help='Minimum cluster size for HDBSCAN clustering algorithm')
    advanced_group.add_argument('--max-per-cluster', type=int, default=3,
                               help='Maximum datasets to select from each cluster')
    advanced_group.add_argument('--api-delay', type=float, default=0.4,
                               help='Delay between API calls (seconds) to respect rate limits')
    advanced_group.add_argument('--info-workers', type=int, default=None,
                               help='Parallel workers for GEO summary fetching (auto-detects based on API key)')
    advanced_group.add_argument('--validation-workers', type=int, default=None,
                               help='Number of parallel workers for dataset validation (auto-detects based on API key)')
    advanced_group.add_argument('--sra-workers', type=int, default=None,
                               help='Number of parallel workers for SRA data fetching (auto-detects based on API key)')
    advanced_group.add_argument('--min-samples', type=int, default=None,
                                help='Minimum number of samples in a neighbourhood for HDBSCAN (None lets HDBSCAN default to min_cluster_size)')
    advanced_group.add_argument('--verbose', action='store_true',
                                help='Enable verbose logging (DEBUG level)')

    args = parser.parse_args()

    # Set up logging in the output directory
    from pathlib import Path
    output_path = Path(args.output)
    # Ensure output directory exists
    output_path.mkdir(parents=True, exist_ok=True)
    setup_logging(verbose=args.verbose)  # Use consistent logging setup

    # Test logging configuration
    logging.info("=" * 60)
    logging.info("🔧 LOGGING CONFIGURATION TEST")
    logging.info(f"📝 Log level: {'DEBUG' if args.verbose else 'INFO'}")
    logging.info(f"⏰ Timestamp test: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    logging.info("=" * 60)

    research_query = args.query
    logging.info(f"Starting optimized dataset identification for query: {research_query}")

    try:
        # Step 1: Extract terms from research query
        logging.info("Extracting search terms...")
        terms = extract_terms(research_query)
        logging.info(f"Extracted terms: {terms.extracted_terms}")
        logging.info(f"Expanded terms: {terms.expanded_terms}")

        # Step 2: Search GEO database
        logging.info("Searching GEO database...")
        search_terms = set(terms.extracted_terms + terms.expanded_terms)
        logging.info(f"Using {len(search_terms)} unique search terms: {', '.join(search_terms)}")

        geo_ids = []
        for term in search_terms:
            logging.info(f"Searching for: {term}")
            term_ids = perform_search(term, args.max_datasets)
            geo_ids.extend(term_ids)
            # Apply rate limiting between searches to follow NCBI guidelines
            time.sleep(args.api_delay)

        unique_geo_ids = list(dict.fromkeys(geo_ids))
        logging.info(f"Found {len(geo_ids)} total GEO IDs, {len(unique_geo_ids)} unique IDs")

        if not unique_geo_ids:
            logging.error("No datasets found in search")
            return

        # Step 3: Get basic dataset information
        basic_datasets_df = get_basic_dataset_info(
            unique_geo_ids,
            args.api_delay,
            workers=args.info_workers,
        )

        if basic_datasets_df.empty:
            logging.error("No valid GSE datasets found")
            return

        # Step 4: Complete dataset info already fetched via streamlined esummary v2.0
        logging.info("Using complete dataset information from streamlined esummary v2.0...")
        enriched_datasets_df = basic_datasets_df  # Already has title, summary, BioProject, species, etc.

        # Step 5: Fetch SRA data for ALL datasets to determine validity (optimal workflow)
        logging.info(f"🔗 Fetching SRA metadata for ALL {len(enriched_datasets_df)} datasets to determine validity...")
        to_fetch = enriched_datasets_df  # Fetch SRA data for ALL datasets

        def fetch_sra_for_dataset(row, api_delay):
            """Fetch SRA data for a single dataset with rate limiting using streamlined approach."""
            acc = row['ID']
            bioproject = row.get('BioProject', '')

            try:
                # Apply rate limiting before making API calls
                api_rate_limiter.wait()

                if bioproject:
                    # Use BioProject directly from esummary (streamlined approach)
                    df_run = fetch_runinfo_from_bioproject(bioproject, api_delay)
                    if not df_run.empty:
                        df_run.insert(0, 'GEO_Accession', acc)
                        logging.info(f"Fetched {len(df_run)} runs for {acc} using BioProject {bioproject}")
                        return df_run, True
                    else:
                        logging.warning(f"No runs found for {acc} with BioProject {bioproject}")
                        return pd.DataFrame(), False
                else:
                    logging.warning(f"No BioProject found for {acc}")
                    return pd.DataFrame(), False

            except Exception as e:
                logging.info(f'Warning, skipping {acc}: {e}')
                return pd.DataFrame(), False

        runs = []
        successful_datasets = 0
        total_runs_fetched = 0

        # Auto-determine optimal SRA worker count if not specified
        sra_workers = args.sra_workers
        if sra_workers is None:
            if os.getenv("ENTREZ_API_KEY"):
                sra_workers = 6  # More aggressive with API key
                logging.info("🚀 API key detected - using 6 SRA workers")
            else:
                sra_workers = 2  # Conservative without API key
                logging.info("⚠️ No API key - using 2 SRA workers")

        # Process SRA fetching in parallel with rate limiting
        logging.info(f"Fetching SRA data using {sra_workers} parallel workers (rate limited)...")
        # Update rate limiter for SRA operations
        api_rate_limiter = APIRateLimiter()
        with ThreadPoolExecutor(max_workers=sra_workers) as executor:
            # Create partial function with fixed api_delay
            fetch_func = partial(fetch_sra_for_dataset, api_delay=args.api_delay)

            # Submit all tasks
            future_to_row = {executor.submit(fetch_func, row): row for _, row in to_fetch.iterrows()}

            # Collect results as they complete
            completed = 0
            total = len(future_to_row)
            for future in as_completed(future_to_row):
                df_result, success = future.result()
                completed += 1

                if success and not df_result.empty:
                    runs.append(df_result)
                    successful_datasets += 1
                    total_runs_fetched += len(df_result)

                # Progress logging
                if completed % 5 == 0 or completed == total:
                    logging.info(f"SRA fetching progress: {completed}/{total} datasets processed")

        # Combine SRA data
        if runs:
            sra_df = pd.concat(runs, ignore_index=True)
            logging.info(f"Successfully processed {successful_datasets}/{len(to_fetch)} GEO datasets")
            logging.info(f"Retrieved metadata for {total_runs_fetched} total SRA runs")
        else:
            sra_df = pd.DataFrame(columns=['GEO_Accession'])
            logging.warning("No SRA runs were successfully fetched")

        # Calculate dataset sizes using runinfo data
        dataset_sizes = {}
        if not sra_df.empty:
            logging.info("🔍 Calculating dataset sizes from runinfo data...")
            dataset_sizes = calculate_dataset_sizes_from_runinfo(sra_df)

        # Add dataset sizes to enriched datasets before merging
        if dataset_sizes:
            enriched_datasets_df['DatasetSizeBytes'] = enriched_datasets_df['ID'].map(dataset_sizes)
            enriched_datasets_df['DatasetSizeGB'] = enriched_datasets_df['DatasetSizeBytes'] / (1024**3)
            # Fill NaN values with 0 for datasets where size couldn't be calculated
            enriched_datasets_df['DatasetSizeBytes'] = enriched_datasets_df['DatasetSizeBytes'].fillna(0)
            enriched_datasets_df['DatasetSizeGB'] = enriched_datasets_df['DatasetSizeGB'].fillna(0.0)
        else:
            enriched_datasets_df['DatasetSizeBytes'] = 0
            enriched_datasets_df['DatasetSizeGB'] = 0.0

        # Merge all GEO datasets with fetched SRA info
        final_results = enriched_datasets_df.merge(
            sra_df,
            left_on='ID',
            right_on='GEO_Accession',
            how='left'
        )
        logging.info(f"Merged GEO and SRA data for {len(final_results)} total records")

        # Step 6: Validate ALL datasets based on complete SRA criteria
        logging.info("Performing dataset validation based on SRA criteria...")
        validation_data = []
        for _, row in final_results.iterrows():
            gse_id = row['ID']
            lib_source = row.get('LibrarySource')
            lib_layout = row.get('LibraryLayout')
            lib_strategy = row.get('LibraryStrategy')

            # Check RNA-seq criteria
            if pd.isna(lib_source) or pd.isna(lib_layout) or pd.isna(lib_strategy):
                valid = False
                reason = "No SRA metadata available"
            elif lib_source != 'TRANSCRIPTOMIC':
                valid = False
                reason = f"Not transcriptomic (LibrarySource: {lib_source})"
            elif lib_strategy != 'RNA-Seq':
                valid = False
                reason = f"Not RNA-seq (LibraryStrategy: {lib_strategy})"
            elif lib_layout != 'PAIRED':
                valid = False
                reason = f"Not paired-end (LibraryLayout: {lib_layout})"
            else:
                valid = True
                reason = "Meets RNA-seq criteria"

            validation_data.append({
                'ID': gse_id,
                'valid_dataset': valid,
                'validation_reason': reason
            })

        validation_df = pd.DataFrame(validation_data).drop_duplicates(subset='ID')
        print(validation_df)
        final_results = final_results.merge(validation_df, on='ID', how='left', validate = 'many_to_one')

        # Step 7: Filter to ONLY valid datasets for expensive operations (clustering & relevance)
        valid_datasets_df = final_results[final_results['valid_dataset'] == True].copy()
        invalid_datasets_df = final_results[final_results['valid_dataset'] == False].copy()

        logging.info(f"Found {len(valid_datasets_df)} valid datasets and {len(invalid_datasets_df)} invalid datasets")
        logging.info(f"🚀 EFFICIENCY GAIN: Will only cluster/score {len(valid_datasets_df)} datasets instead of {len(final_results)} total")
        logging.info(f"⚡ Avoiding expensive operations on {len(invalid_datasets_df)} invalid datasets")

        if len(valid_datasets_df) == 0:
            logging.error("No valid datasets found - cannot proceed with clustering/relevance assessment")
            # Still output all datasets for transparency
            final_results['RelevanceScore'] = None
            final_results['Run1Score'] = None
            final_results['Run1Justification'] = None
            final_results['Run2Score'] = None
            final_results['Run2Justification'] = None
            final_results['Run3Score'] = None
            final_results['Run3Justification'] = None
        else:
            # Step 8: Embed and cluster ONLY valid datasets (efficient approach)
            start_time = time.time()
            logging.info(f"🧮 Embedding {len(valid_datasets_df)} valid datasets (skipping {len(invalid_datasets_df)} invalid)...")
            embedded_df = embed_datasets(valid_datasets_df)

            logging.info("🔗 Clustering valid datasets...")
            clustered_df = cluster_datasets(embedded_df)

            # Step 9: Select representatives from valid datasets only
            logging.info("🎯 Selecting representative datasets for assessment...")
            representatives_df = select_representative_datasets(clustered_df, args.max_evaluate or 10)

            # Step 10: Assess relevance of representatives only
            logging.info(f"📊 Assessing relevance of {len(representatives_df)} representative valid datasets...")
            assessed_df = asyncio.run(repeated_relevance(representatives_df, research_query, repeats=3, batch_size=10, openai_api_jobs=4))

            end_time = time.time()
            processing_time = end_time - start_time
            logging.info(f"⚡ Expensive operations completed in {processing_time:.1f}s for {len(valid_datasets_df)} valid datasets")
            estimated_saved_time = processing_time * (len(invalid_datasets_df) / len(valid_datasets_df)) if len(valid_datasets_df) > 0 else 0
            logging.info(f"💰 Estimated time saved by skipping invalid datasets: {estimated_saved_time:.1f}s")

            # Step 11: Merge assessment results back to valid datasets
            valid_with_scores = valid_datasets_df.merge(assessed_df, on='ID', how='left')

            # Combine valid (with scores) and invalid datasets for final output
            final_results = pd.concat([
                valid_with_scores,
                invalid_datasets_df
            ], ignore_index=True)

        # Remove embedding column to prevent CSV malformation
        if 'embedding' in final_results.columns:
            final_results = final_results.drop('embedding', axis=1)

        # Create final dataframe with all columns to match original implementation
        final = pd.DataFrame()
        final['ID'] = final_results['ID']
        final['Title'] = final_results.get('Title', '')
        final['Summary'] = final_results.get('Summary', '')
        final['Accession'] = final_results['ID']
        final['Species'] = final_results.get('Species', 'Unknown')
        final['Date'] = final_results.get('Date', '')
        final['NumSamples'] = final_results.get('NumSamples', 0)
        final['PrimaryPubMedID'] = final_results.get('PrimaryPubMedID', None)
        final['AllPubMedIDs'] = final_results.get('AllPubMedIDs', None)
        final['RelevanceScore'] = final_results.get('RelevanceScore', None)

        # ==================== OPTIONAL: REPEATED ASSESSMENT COLUMNS ====================
        # The following columns show individual assessment runs (can be removed for simplified output)
        final['Run1Score'] = final_results.get('Run1Score', None)
        final['Run1Justification'] = final_results.get('Run1Justification', None)
        final['Run2Score'] = final_results.get('Run2Score', None)
        final['Run2Justification'] = final_results.get('Run2Justification', None)
        final['Run3Score'] = final_results.get('Run3Score', None)
        final['Run3Justification'] = final_results.get('Run3Justification', None)
        # ===============================================================================

        # ==================== OPTIONAL: DATASET SIZE COLUMNS ====================
        # The following columns show dataset storage size (can be removed for simplified output)
        final['DatasetSizeBytes'] = final_results.get('DatasetSizeBytes', 0)
        final['DatasetSizeGB'] = final_results.get('DatasetSizeGB', 0.0)
        # =========================================================================

        # ==================== OPTIONAL: SRA METADATA COLUMNS ====================
        # The following columns show SRA technical details (can be removed for simplified output)
        final['LibrarySource'] = final_results.get('LibrarySource', None)
        final['LibraryLayout'] = final_results.get('LibraryLayout', None)
        final['LibraryStrategy'] = final_results.get('LibraryStrategy', None)
        # =========================================================================

        # Select and flag validity based on SRA criteria (matching original implementation)
        final['Valid'] = final_results.get('valid_dataset', False).apply(lambda x: 'Yes' if x else 'No')

        # Add validation reason column for debugging
        final['ValidationReason'] = final_results.get('validation_reason', 'Unknown')

        # ==================== OPTIONAL: CONVENIENCE URL COLUMNS ====================
        # The following columns create clickable URLs (can be removed for simplified output)
        def create_pubmed_url(pmid):
            if pd.notna(pmid):
                return f"https://pubmed.ncbi.nlm.nih.gov/{int(pmid)}/"
            return None

        def create_geo_url(accession):
            if pd.notna(accession):
                return f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}"
            return None

        final['PubMedURL'] = final['PrimaryPubMedID'].apply(create_pubmed_url)
        final['GEOURL'] = final['Accession'].apply(create_geo_url)
        # ============================================================================

        final = final.drop_duplicates(subset=['Accession'])

        # Save to output directory
        main_output_file = output_path / "Dataset_identification_result.csv"

        final.to_csv(main_output_file, index=False)
        message = f'Saved combined table to {main_output_file}'
        logging.info(message)
        print(message)

        # Save a separate file with the selected datasets and their clusters if clustering was used
        if args.use_clustering and 'Cluster' in final_results.columns:
            selected_path = output_path / 'selected_datasets.csv'

            # Create clustering output with simplified columns
            clustering_df = pd.DataFrame()
            clustering_df['Accession'] = final['Accession']
            clustering_df['Title'] = final['Title']
            clustering_df['Cluster'] = final_results.get('Cluster', -1)
            clustering_df['RelevanceScore'] = final['RelevanceScore']
            clustering_df = clustering_df.sort_values(['Cluster', 'RelevanceScore'], ascending=[True, False])
            clustering_df.to_csv(selected_path, index=False)
            message = f'Saved selected datasets with cluster information to {selected_path}'
            logging.info(message)
            print(message)

        if args.multi_dataset_csv:
            multi_df = final[['Accession', 'Species', 'PrimaryPubMedID', 'RelevanceScore', 'Valid', 'DatasetSizeBytes', 'DatasetSizeGB']].copy()
            multi_df = multi_df[multi_df['Valid'] == 'Yes']
            multi_df = multi_df.dropna(subset=['RelevanceScore'])  # Only include assessed datasets
            multi_df = multi_df.sort_values('RelevanceScore', ascending=False)
            multi_df = multi_df.rename(columns={'Species': 'organism'})

            # Set path for the multi-dataset CSV in the output directory
            multi_csv_path = output_path / 'batch_analysis_input.csv'

            multi_df.to_csv(multi_csv_path, index=False)
            message = f"Generated batch analysis input CSV at {multi_csv_path}"
            logging.info(message)
            print(message)

        # Print summary results
        print("\n" + "="*80)
        print("OPTIMIZED DATASET IDENTIFICATION RESULTS")
        print("="*80)
        print(f"Research Query: {research_query}")
        print(f"\n📊 WORKFLOW SUMMARY:")
        print(f"   Total datasets found: {len(final)}")
        print(f"   Valid datasets (RNA-seq, paired-end, transcriptomic): {len(final[final['Valid'] == 'Yes'])}")
        print(f"   Invalid datasets: {len(final[final['Valid'] == 'No'])}")
        print(f"   Clustered & assessed (from valid only): {len(final.dropna(subset=['RelevanceScore']))}")

        print(f"\n🔍 EFFICIENCY GAINS:")
        print(f"   ✅ Fetched complete info for ALL {len(final)} datasets")
        print(f"   ✅ Only clustered/scored {len(final[final['Valid'] == 'Yes'])} valid datasets")
        print(f"   ✅ Saved compute by skipping {len(final[final['Valid'] == 'No'])} invalid datasets")

        # Show some invalid datasets with reasons
        invalid_df = final[final['Valid'] == 'No']
        if len(invalid_df) > 0:
            print(f"\n❌ Sample invalid datasets (showing first 5):")
            for _, row in invalid_df.head(5).iterrows():
                print(f"   {row['Accession']}: {row['ValidationReason']}")

        print(f"\n🏆 Top {min(len(final.dropna(subset=['RelevanceScore'])), args.max_evaluate or 10)} Assessed Results:")
        print("-" * 80)

        assessed_final = final.dropna(subset=['RelevanceScore']).sort_values('RelevanceScore', ascending=False)
        for i, (_, row) in enumerate(assessed_final.head(args.max_evaluate or 10).iterrows()):
            print(f"\n{i+1}. Dataset: {row['Accession']}")
            print(f"   Title: {row['Title']}")
            print(f"   Samples: {row['NumSamples']}")
            print(f"   Species: {row['Species']}")
            print(f"   Valid for analysis: {'✅ Yes' if row['Valid'] == 'Yes' else '❌ No'}")
            if row['Valid'] == 'No':
                print(f"   Validation reason: {row['ValidationReason']}")
            print(f"   Relevance: {row['RelevanceScore']}/10")
            print(f"   Dataset Size: {row['DatasetSizeGB']:.2f} GB")
            print(f"   GEO URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={row['Accession']}")

    except Exception as e:
        logging.error(f"Error in main execution: {e}")
        raise

if __name__ == "__main__":
    main()

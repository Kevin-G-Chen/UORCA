#!/usr/bin/env python3
"""
Dataset Identification Script
============================

This script identifies relevant RNA-seq datasets from GEO based on research queries.
Key features:
- Automatic API rate limiting based on API key presence
- Parallel processing for efficient data retrieval
- Dataset-level validation with RNA-seq criteria
- AI-powered clustering and relevance scoring (enabled by default)
- Multi-dataset CSV output for batch analysis (enabled by default)

Requirements:
- ENTREZ_EMAIL environment variable (required by NCBI guidelines)
- OPENAI_API_KEY environment variable for AI features
- Optional: ENTREZ_API_KEY for faster processing

Workflow:
1. Extract search terms and query GEO database
2. Fetch complete dataset information using esummary v2.0
3. Retrieve SRA metadata for validation
4. Validate datasets based on RNA-seq criteria (>3 paired-end transcriptomic samples)
5. Cluster valid datasets and assess relevance using AI
6. Generate results CSV and batch analysis input

Usage:
    python DatasetIdentification.py --query "research query"
"""

import argparse
import asyncio
import datetime
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
from tqdm import tqdm

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

# Configure Entrez with explicit key setting
def _configure_entrez():
    """Configure Entrez with current environment variables."""
    if os.getenv("ENTREZ_EMAIL"):
        Entrez.email = os.getenv("ENTREZ_EMAIL")
    if os.getenv("ENTREZ_API_KEY"):
        Entrez.api_key = os.getenv("ENTREZ_API_KEY")
        # Verify it was set
        if not hasattr(Entrez, 'api_key') or not Entrez.api_key:
            logging.warning("Failed to set Entrez.api_key - NCBI rate limiting may not work as expected")
    else:
        # Explicitly set to None if not available
        Entrez.api_key = None

# Initial configuration at module load
_configure_entrez()

# OpenAI client initialization with validation
openai_api_key = os.getenv("OPENAI_API_KEY")
if not openai_api_key:
    logging.error("❌ OPENAI_API_KEY not found in environment variables")
    logging.error("🔧 Dataset identification requires OpenAI API key for AI-powered relevance assessment")
    logging.error("📝 Please add OPENAI_API_KEY=sk-proj-your-key-here to your .env file")
    logging.error("🔗 Get your API key at: https://platform.openai.com/api-keys")
    raise EnvironmentError(
        "OPENAI_API_KEY is required for dataset identification. "
        "Please set this environment variable in your .env file."
    )

try:
    client = OpenAI(api_key=openai_api_key)
    # Test the client with a minimal request to validate the key
    logging.info("✅ OpenAI API key validated successfully")
except Exception as e:
    logging.error(f"❌ Failed to initialize OpenAI client: {e}")
    logging.error("🔧 Please check that your OPENAI_API_KEY is valid")
    logging.error("🔗 Verify your API key at: https://platform.openai.com/api-keys")
    raise EnvironmentError(f"Failed to initialize OpenAI client: {e}")

# Rate limiting for API calls
class APIRateLimiter:
    """Thread-safe rate limiter for API calls with dynamic rates based on API key."""
    def __init__(self, base_delay: float = None):
        # Auto-detect optimal delay based on API key presence
        if base_delay is None:
            if os.getenv("ENTREZ_API_KEY"):
                # With API key: target 7 requests/second = 0.143 second delay
                self.min_delay = 1.0 / 7.0  # Conservative rate limiting
            else:
                # Without API key: target 2 requests/second = 0.5 second delay
                self.min_delay = 1.0 / 2.0  # Conservative rate limiting
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
def setup_logging(verbose: bool = False, suppress_sra_warnings: bool = True) -> None:
    """Setup logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    format_str = "%(asctime)s - %(levelname)s - %(message)s"

    # Create logs directory if it doesn't exist
    logs_dir = Path("logs/identification_logs")
    logs_dir.mkdir(parents=True, exist_ok=True)

    # Clear any existing handlers to avoid conflicts
    root_logger = logging.getLogger()
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)

    # Create timestamped log file
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = logs_dir / f"dataset_identification_{timestamp}.log"

    # Configure logging to both file and console with forced configuration
    file_handler = logging.FileHandler(log_file)
    console_handler = logging.StreamHandler(sys.stdout)

    logging.basicConfig(
        level=level,
        format=format_str,
        handlers=[file_handler, console_handler],
        force=True  # Force reconfiguration even if logging was already configured
    )

    # Enforce retention: keep only the 5 most recent identification log files
    try:
        existing = sorted(
            [p for p in logs_dir.glob("dataset_identification_*.log") if p.is_file()],
            key=lambda p: p.stat().st_mtime,
        )
        # If more than 5, delete the oldest ones
        if len(existing) > 5:
            to_delete = existing[: len(existing) - 5]
            for old in to_delete:
                try:
                    old.unlink()
                except Exception:
                    pass
    except Exception:
        pass

    # Also configure warnings to use the same format
    logging.captureWarnings(True)
    warnings_logger = logging.getLogger('py.warnings')
    warnings_logger.setLevel(level)

    # Optionally suppress SRA-related warnings by setting a higher threshold
    if suppress_sra_warnings and not verbose:
        # Create a custom filter to suppress specific SRA warnings
        class SRAWarningFilter(logging.Filter):
            def filter(self, record):
                # Suppress "No SRA records found" warnings unless in verbose mode
                if "No SRA records found" in record.getMessage():
                    return False
                return True

        # Apply filter to both file and console handlers
        for handler in logging.getLogger().handlers:
            handler.addFilter(SRAWarningFilter())

for name in ("openai", "openai._base_client", "httpx"):
    logging.getLogger(name).setLevel(logging.WARNING)

def load_prompt(file_path: str) -> str:
    """Load prompt from file."""
    return Path(file_path).read_text().strip()

# Query config management for Streamlit integration
def save_query_config(query: str) -> None:
    """Save the dataset identification query to a config file for Streamlit app."""
    config_dir = Path("main_workflow/reporting/.config")
    config_dir.mkdir(parents=True, exist_ok=True)

    config_file = config_dir / "dataset_query.json"
    config_data = {
        "query": query,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S")
    }

    with open(config_file, 'w') as f:
        json.dump(config_data, f, indent=2)

def load_query_config() -> Optional[str]:
    """Load the dataset identification query from config file."""
    config_file = Path("main_workflow/reporting/.config/dataset_query.json")

    if config_file.exists():
        try:
            with open(config_file, 'r') as f:
                config_data = json.load(f)
                return config_data.get("query")
        except (json.JSONDecodeError, KeyError):
            return None
    return None

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

def call_openai_json(prompt: str, user_input: str, response_format: BaseModel, model: str = "gpt-5-mini") -> Any:
    """Make OpenAI API call with JSON response format."""
    response = client.beta.chat.completions.parse(
        model=model,
        messages=[
            {"role": "system", "content": prompt},
            {"role": "user", "content": user_input}
        ],
        response_format=response_format
    )
    return response.choices[0].message.parsed

def extract_terms(research_query: str, model: str = "gpt-5-mini") -> ExtractedTerms:
    """Extract search terms from research query."""
    prompt = load_prompt("./main_workflow/prompts/dataset_identification/extract_terms.txt")
    return call_openai_json(prompt, research_query, ExtractedTerms, model)

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

    return geo_ids

def get_basic_dataset_info(
    geo_ids: List[str],
    api_delay: float = 0.4,
) -> pd.DataFrame:
    """Get complete dataset information in parallel using esummary v2.0."""

    unique_ids = list(dict.fromkeys(geo_ids))
    workers = 6 if os.getenv("ENTREZ_API_KEY") else 1

    logging.info(f"Fetching complete dataset information for {len(unique_ids)} unique datasets...")

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

    is_main_thread = threading.current_thread() is threading.main_thread()
    with ThreadPoolExecutor(max_workers=workers) as executor, \
         tqdm(total=len(unique_ids), desc="Fetching GEO summaries", unit="dataset", disable=not is_main_thread) as pbar:

        futures = {executor.submit(fetch_info, gid): gid for gid in unique_ids}
        for future in as_completed(futures):
            result = future.result()
            if result:
                datasets.append(result)
            pbar.update(1)

    df = pd.DataFrame(datasets)
    if not df.empty:
        df = df.drop_duplicates(subset=["ID"]).reset_index(drop=True)

    logging.info(f"Found {len(df)} unique GSE datasets with complete information")
    return df

# Streamlined SRA processing functions (old XML parsing functions removed)
def fetch_runinfo_from_bioproject(bioproject: str, api_delay: float = 0.4, suppress_missing_warnings: bool = True) -> pd.DataFrame:
    """Fetch RunInfo data using streamlined BioProject approach."""
    try:
        # 1. Search SRA and keep the server-side history
        srch = Entrez.read(Entrez.esearch(db="sra",
                                        term=f"{bioproject}[BioProject]",
                                        usehistory="y"))
        count, webenv, qk = int(srch["Count"]), srch["WebEnv"], srch["QueryKey"]

        if count == 0:
            if not suppress_missing_warnings:
                logging.warning(f"No SRA records found for BioProject {bioproject}")
            else:
                logging.debug(f"No SRA records found for BioProject {bioproject}")
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
        return runs

    except Exception as e:
        logging.warning(f"Error fetching RunInfo for BioProject {bioproject}: {e}")
        return pd.DataFrame()



def calculate_dataset_sizes_from_runinfo(sra_df: pd.DataFrame) -> Dict[str, int]:
    """
    Calculate dataset sizes using only valid RNA-seq samples (TRANSCRIPTOMIC, RNA-Seq, PAIRED).

    Returns:
        Dictionary mapping GEO_Accession to total size in bytes (valid samples only)
    """
    if sra_df.empty or 'size_MB' not in sra_df.columns or 'GEO_Accession' not in sra_df.columns:
        logging.warning("Cannot calculate dataset sizes: missing required columns in runinfo data")
        return {}

    # Filter to only valid RNA-seq samples
    valid_samples = sra_df[
        (sra_df['LibrarySource'] == 'TRANSCRIPTOMIC') &
        (sra_df['LibraryStrategy'] == 'RNA-Seq') &
        (sra_df['LibraryLayout'] == 'PAIRED')
    ].copy()

    if valid_samples.empty:
        logging.info("No valid RNA-seq samples found for size calculation")
        return {}

    # Convert size_MB to bytes and group by GEO accession
    valid_samples['size_MB'] = pd.to_numeric(valid_samples['size_MB'], errors='coerce').fillna(0)
    valid_samples['size_bytes'] = (valid_samples['size_MB'] * 1024 * 1024).astype(int)

    dataset_sizes = valid_samples.groupby('GEO_Accession')['size_bytes'].sum().to_dict()

    # Log results
    total_datasets = len(dataset_sizes)
    total_size_gb = sum(dataset_sizes.values()) / (1024**3)

    return dataset_sizes




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
    embeddings = []

    # Use progress bar for embedding generation (disable if not main thread)
    is_main_thread = threading.current_thread() is threading.main_thread()
    for _, row in tqdm(datasets_df.iterrows(), total=len(datasets_df), desc="Generating embeddings", unit="dataset", disable=not is_main_thread):
        # Combine title and summary for embedding
        text = f"{row.get('Title', '')} {row.get('Summary', '')}"
        embedding = get_embedding(text)
        embeddings.append(embedding)

        time.sleep(0.1)  # Small delay for API rate limiting

    # Add embeddings to dataframe
    datasets_df = datasets_df.copy()
    datasets_df['embedding'] = embeddings

    return datasets_df

def cluster_datasets(datasets_df: pd.DataFrame, n_clusters: int = None, cluster_divisor: int = 10) -> pd.DataFrame:
    """Cluster datasets based on embeddings with simple divisor-based cluster count."""
    if len(datasets_df) <= 1:
        datasets_df = datasets_df.copy()
        datasets_df['Cluster'] = 0
        return datasets_df

    if n_clusters is None:
        n_clusters = max(1, len(datasets_df) // cluster_divisor)  # Simple: total/divisor

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

def select_representative_datasets(datasets_df: pd.DataFrame, max_assess: int = 300) -> pd.DataFrame:
    """Select datasets for assessment with budget constraint, distributed proportionally across clusters."""
    n_clusters = datasets_df['Cluster'].nunique()

    if len(datasets_df) <= max_assess:
        logging.info(f"Assessing all {len(datasets_df)} datasets (within budget of {max_assess})")
        return datasets_df  # Assess all if within budget

    # Distribute assessment budget proportionally across clusters
    datasets_per_cluster = max(1, max_assess // n_clusters)

    logging.info(f"Selecting up to {datasets_per_cluster} datasets from each of {n_clusters} clusters (budget: {max_assess})")

    representatives = []
    total_selected = 0

    for cluster_id in sorted(datasets_df['Cluster'].unique()):
        cluster_datasets = datasets_df[datasets_df['Cluster'] == cluster_id]
        n_to_take = min(datasets_per_cluster, len(cluster_datasets))
        cluster_representatives = cluster_datasets.head(n_to_take)
        representatives.extend(cluster_representatives.to_dict('records'))
        total_selected += n_to_take

    logging.info(f"Selected {total_selected} total datasets for relevance assessment")
    return pd.DataFrame(representatives).reset_index(drop=True)


async def assess_subbatch(df: pd.DataFrame, query: str, schema, key: str, rep: int, idx: int, total_batches: int, sem: asyncio.Semaphore, model: str) -> List[Assessment]:
    """Assess a subbatch of datasets."""
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

            assessments = call_openai_json(prompt, user_input, Assessments, model)
            return assessments.assessments
        except Exception as e:
            logging.warning(f"Error in assessment batch {idx+1}/{total_batches} (rep {rep+1}): {e}")
            # Return default assessments
            return [Assessment(ID=row['ID'], RelevanceScore=5, Justification="Assessment failed") for _, row in df.iterrows()]

async def repeated_relevance(df: pd.DataFrame, query: str, repeats: int = 3, batch_size: int = 10, openai_api_jobs: int = 3, model: str = "gpt-5-mini") -> pd.DataFrame:
    """Perform repeated relevance assessment with averaging."""
    logging.info(f"Starting relevance scoring: {repeats} repetitions, batch size {batch_size}, parallel API jobs: {openai_api_jobs}")
    sem = asyncio.Semaphore(openai_api_jobs)

    # create sub-batches
    batches = [df.iloc[i:i+batch_size] for i in range(0, len(df), batch_size)]
    total_batches = len(batches)
    total_tasks = repeats * total_batches

    # Create tasks with progress tracking
    logging.info(f"Processing {total_tasks} assessment batches...")

    async def assess_with_progress(rep, idx, sub, pbar):
        """Wrapper to update progress bar when task completes."""
        result = await assess_subbatch(sub, query, None, "assessments", rep, idx, total_batches, sem, model)
        pbar.update(1)
        pbar.set_postfix({"Rep": f"{rep+1}/{repeats}", "Batch": f"{(idx+1)+rep*total_batches}/{total_tasks}"})
        return result

    is_main_thread = threading.current_thread() is threading.main_thread()
    with tqdm(total=total_tasks, desc="Assessing dataset relevance", unit="batch", disable=not is_main_thread) as pbar:
        tasks = []
        for rep in range(repeats):
            for idx, sub in enumerate(batches):
                tasks.append(assess_with_progress(rep, idx, sub, pbar))

        # Execute all tasks and gather results
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


def main():
    """
    Main function for dataset identification.

    Process:
    1. Extract search terms from research query
    2. Search GEO database for relevant datasets
    3. Validate datasets for RNA-seq compatibility
    4. Cluster valid datasets (count = total_datasets / cluster_divisor)
    5. Select representative datasets for assessment (budget distributed across clusters)
    6. Assess relevance of representatives using AI
    7. Generate comprehensive results CSV (all datasets) and batch analysis CSV (filtered)
    """
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
    search_group.add_argument('-m', '--max-per-term', type=int, default=500,
                             help='Maximum datasets to retrieve per search term')
    search_group.add_argument('-a', '--max-assess', type=int, default=300,
                             help='Maximum number of datasets to assess for relevance (distributed across clusters)')

    # === Advanced Parameters ===
    advanced_group = parser.add_argument_group('Advanced Options', 'Fine-tune algorithm behavior (expert users)')
    advanced_group.add_argument('--cluster-divisor', type=int, default=10,
                               help='Divisor for cluster count (total_datasets / divisor). Smaller values = more clusters.')
    advanced_group.add_argument('-r', '--rounds', type=int, default=3,
                               help='Number of independent relevance scoring rounds for reliability')
    advanced_group.add_argument('-b', '--batch-size', type=int, default=20,
                               help='Datasets per AI evaluation batch (affects memory usage)')
    advanced_group.add_argument('--model', type=str, default='gpt-5-mini',
                                   help='OpenAI model to use for relevance assessment')
    advanced_group.add_argument('-v', '--verbose', action='store_true',
                               help='Enable verbose logging (DEBUG level)')

    args = parser.parse_args()

    # Check for required email (Entrez guidelines)
    if not os.getenv("ENTREZ_EMAIL"):
        logging.error("Email is required for NCBI Entrez API access. Please set ENTREZ_EMAIL environment variable.")
        logging.error("Please set the ENTREZ_EMAIL environment variable with your email address.")
        logging.error("This is required by NCBI guidelines for API usage.")
        return

    # Auto-determine API delay based on API key presence
    api_delay = 0.25 if os.getenv("ENTREZ_API_KEY") else 0.6

    # Set up logging in the output directory
    from pathlib import Path
    output_path = Path(args.output)
    # Ensure output directory exists
    output_path.mkdir(parents=True, exist_ok=True)
    setup_logging(verbose=args.verbose, suppress_sra_warnings=not args.verbose)

    # Reconfigure Entrez to ensure API key is set in this process/thread
    _configure_entrez()

    # Log API key detection status and rate limiting configuration
    if os.getenv("ENTREZ_API_KEY"):
        logging.info("NCBI API key detected - using faster rate limits")
        logging.info("Rate Limiting Configuration:")
        logging.info(f"API delay: {1.0/7.0:.3f}s (~7 req/sec)")
        logging.info(f"GEO workers: 6, SRA workers: 6")
    else:
        logging.info("No NCBI API key found - using standard rate limits")
        logging.info("Rate Limiting Configuration:")
        logging.info(f"API delay: {1.0/2.0:.3f}s (~2 req/sec)")
        logging.info(f"GEO workers: 1, SRA workers: 1")

    research_query = args.query
    logging.info(f"Starting dataset identification for query: {research_query}")

    # Save query to config file for Streamlit app
    save_query_config(research_query)

    # Track timing for metadata
    start_time = datetime.datetime.now()

    try:
        # Step 1: Extract terms from research query
        logging.info("Determining search terms...")
        terms = extract_terms(research_query, args.model)

        # Step 2: Search GEO database
        search_terms = set(terms.extracted_terms + terms.expanded_terms)
        logging.info(f"Searching GEO using {len(search_terms)} unique search terms: {', '.join(search_terms)}")

        geo_ids = []
        # Detect if running in a background thread (e.g., from GUI) - disable tqdm to avoid BrokenPipeError
        is_main_thread = threading.current_thread() is threading.main_thread()
        search_iterator = tqdm(search_terms, desc="Searching GEO database", unit="term", disable=not is_main_thread)

        for term in search_iterator:
            term_ids = perform_search(term, args.max_per_term)
            geo_ids.extend(term_ids)
            # Apply rate limiting between searches to follow NCBI guidelines
            time.sleep(api_delay)

        unique_geo_ids = list(dict.fromkeys(geo_ids))
        logging.info(f"Found {len(geo_ids)} total GEO IDs, {len(unique_geo_ids)} unique IDs")

        if not unique_geo_ids:
            logging.error("No datasets found in search")
            return

        # Step 3: Get basic dataset information
        basic_datasets_df = get_basic_dataset_info(
            unique_geo_ids,
            api_delay,
        )

        if basic_datasets_df.empty:
            logging.error("No valid GSE datasets found")
            return

        # Step 4: Use complete dataset information from esummary v2.0
        enriched_datasets_df = basic_datasets_df

        # Step 5: Fetch SRA data for ALL datasets to determine validity (optimal workflow)
        logging.info(f"Fetching SRA metadata for {len(enriched_datasets_df)} datasets to determine validity...")
        to_fetch = enriched_datasets_df  # Fetch SRA data for ALL datasets

        def fetch_sra_for_dataset(row, api_delay):
            """Fetch SRA data for a single dataset with rate limiting using streamlined approach."""
            acc = row['ID']
            bioproject = row.get('BioProject', '')

            try:
                # Apply rate limiting before making API calls
                api_rate_limiter.wait()

                if bioproject:
                    df_run = fetch_runinfo_from_bioproject(bioproject, api_delay, suppress_missing_warnings=True)
                    if not df_run.empty:
                        df_run.insert(0, 'GEO_Accession', acc)
                        return df_run, True
                    else:
                        return pd.DataFrame(), False
                else:
                    logging.warning(f"No BioProject found for {acc}")
                    return pd.DataFrame(), False

            except Exception as e:
                logging.warning(f'Error processing {acc}: {e}')
                return pd.DataFrame(), False

        runs = []
        successful_datasets = 0
        total_runs_fetched = 0

        # Auto-determine optimal SRA worker count based on API key
        if os.getenv("ENTREZ_API_KEY"):
            sra_workers = 6  # Increased workers with conservative rate limiting
        else:
            sra_workers = 1  # Very conservative without API key
        api_rate_limiter = APIRateLimiter()

        # Track timing for performance summary
        sra_start_time = time.time()

        with ThreadPoolExecutor(max_workers=sra_workers) as executor:
            fetch_func = partial(fetch_sra_for_dataset, api_delay=api_delay)
            future_to_row = {executor.submit(fetch_func, row): row for _, row in to_fetch.iterrows()}

            # Use tqdm progress bar for SRA fetching (disable if not main thread)
            is_main_thread = threading.current_thread() is threading.main_thread()
            with tqdm(total=len(future_to_row), desc="Fetching SRA data", unit="dataset", disable=not is_main_thread) as pbar:
                for future in as_completed(future_to_row):
                    df_result, success = future.result()
                    pbar.update(1)

                    if success and not df_result.empty:
                        runs.append(df_result)
                        successful_datasets += 1
                        total_runs_fetched += len(df_result)
                        pbar.set_postfix(successful=successful_datasets)

        # Performance summary
        sra_end_time = time.time()
        sra_elapsed_time = sra_end_time - sra_start_time
        datasets_per_second = len(to_fetch) / sra_elapsed_time if sra_elapsed_time > 0 else 0
        success_rate = (successful_datasets / len(to_fetch)) * 100 if len(to_fetch) > 0 else 0

        if runs:
            sra_df = pd.concat(runs, ignore_index=True)
        else:
            sra_df = pd.DataFrame(columns=['GEO_Accession'])
            logging.warning("No SRA runs were successfully fetched")

        # Calculate dataset sizes using runinfo data
        dataset_sizes = {}
        if not sra_df.empty:
            logging.info("Calculating dataset sizes from runinfo data...")
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

        # Step 6: Validate ALL datasets based on complete SRA criteria (dataset-level validation)
        validation_data = []

        # Group by dataset (GSE ID) for proper dataset-level validation
        grouped_datasets = list(final_results.groupby('ID'))
        is_main_thread = threading.current_thread() is threading.main_thread()
        for gse_id, group in tqdm(grouped_datasets, desc="Validating datasets", unit="dataset", disable=not is_main_thread):
            # Count samples that meet RNA-seq criteria
            rnaseq_samples = 0
            total_samples = len(group)

            # Track all library types present in this dataset
            lib_sources = set()
            lib_layouts = set()
            lib_strategies = set()

            # Check each sample in the dataset
            for _, row in group.iterrows():
                lib_source = row.get('LibrarySource')
                lib_layout = row.get('LibraryLayout')
                lib_strategy = row.get('LibraryStrategy')

                # Track all library types (for reporting)
                if pd.notna(lib_source):
                    lib_sources.add(lib_source)
                if pd.notna(lib_layout):
                    lib_layouts.add(lib_layout)
                if pd.notna(lib_strategy):
                    lib_strategies.add(lib_strategy)

                # Check if this sample meets RNA-seq criteria
                if (pd.notna(lib_source) and lib_source == 'TRANSCRIPTOMIC' and
                    pd.notna(lib_strategy) and lib_strategy == 'RNA-Seq' and
                    pd.notna(lib_layout) and lib_layout == 'PAIRED'):
                    rnaseq_samples += 1

            # Dataset-level validation: need >3 samples meeting RNA-seq criteria
            if total_samples == 0:
                valid = False
                reason = "No SRA metadata available"
            elif rnaseq_samples == 0:
                reason_parts = []
                if 'TRANSCRIPTOMIC' not in lib_sources:
                    reason_parts.append(f"LibrarySource: {'/'.join(lib_sources) if lib_sources else 'None'}")
                if 'RNA-Seq' not in lib_strategies:
                    reason_parts.append(f"LibraryStrategy: {'/'.join(lib_strategies) if lib_strategies else 'None'}")
                if 'PAIRED' not in lib_layouts:
                    reason_parts.append(f"LibraryLayout: {'/'.join(lib_layouts) if lib_layouts else 'None'}")
                valid = False
                reason = f"No RNA-seq samples found ({'; '.join(reason_parts)})"
            elif rnaseq_samples <= 3:
                valid = False
                reason = f"Insufficient RNA-seq samples: {rnaseq_samples} (need >3)"
            else:
                valid = True
                # Include info about mixed library types if present
                mixed_info = []
                if len(lib_sources) > 1:
                    mixed_info.append(f"Sources: {'/'.join(lib_sources)}")
                if len(lib_strategies) > 1:
                    mixed_info.append(f"Strategies: {'/'.join(lib_strategies)}")
                if len(lib_layouts) > 1:
                    mixed_info.append(f"Layouts: {'/'.join(lib_layouts)}")

                if mixed_info:
                    reason = f"Valid: {rnaseq_samples}/{total_samples} RNA-seq samples (mixed types: {'; '.join(mixed_info)})"
                else:
                    reason = f"Valid: {rnaseq_samples}/{total_samples} RNA-seq samples"

            validation_data.append({
                'ID': gse_id,
                'valid_dataset': valid,
                'validation_reason': reason,
                'rnaseq_samples': rnaseq_samples,
                'total_samples': total_samples
            })

        validation_df = pd.DataFrame(validation_data)

        # Merge validation results back to final_results
        final_results = final_results.merge(validation_df, on='ID', how='left', validate='many_to_one')

        # Step 7: Create dataset-level summary
        # Create one row per dataset (remove duplicates from multiple samples)
        dataset_summary = final_results.groupby('ID').agg({
            'Title': 'first',
            'Summary': 'first',
            'Accession': 'first',
            'Species': 'first',
            'Date': 'first',
            'NumSamples': 'first',
            'PrimaryPubMedID': 'first',
            'DatasetSizeBytes': 'first',
            'DatasetSizeGB': 'first',
            'valid_dataset': 'first',
            'validation_reason': 'first',
            'rnaseq_samples': 'first',
            'total_samples': 'first'
        }).reset_index()

        # Filter to valid and invalid datasets for processing
        valid_datasets_df = dataset_summary[dataset_summary['valid_dataset'] == True].copy()
        invalid_datasets_df = dataset_summary[dataset_summary['valid_dataset'] == False].copy()

        logging.info(f"Found {len(valid_datasets_df)} valid datasets and {len(invalid_datasets_df)} invalid datasets")

        if len(valid_datasets_df) == 0:
            logging.error("No valid datasets found - cannot proceed with clustering/relevance assessment")
            # Still output all datasets for transparency
            final_results = pd.concat([valid_datasets_df, invalid_datasets_df], ignore_index=True)
        else:
            # Step 8: Embed and cluster valid datasets
            logging.info(f"Embedding {len(valid_datasets_df)} valid datasets...")
            embedded_df = embed_datasets(valid_datasets_df)

            logging.info("Clustering valid datasets...")
            clustered_df = cluster_datasets(embedded_df, cluster_divisor=args.cluster_divisor)

            # Step 9: Select representative datasets
            logging.info("Selecting representative datasets for assessment...")
            representatives_df = select_representative_datasets(clustered_df, args.max_assess)

            # Step 10: Assess relevance of representatives
            logging.info(f"Assessing relevance of {len(representatives_df)} representative datasets...")
            assessed_df = asyncio.run(repeated_relevance(representatives_df, research_query, repeats=args.rounds, batch_size=args.batch_size, openai_api_jobs=4, model=args.model))

            # Step 11: Merge assessment results back to valid datasets
            valid_with_scores = valid_datasets_df.merge(assessed_df, on='ID', how='left')

            # Combine valid (with scores) and invalid datasets for final output
            final_results = pd.concat([
                valid_with_scores,
                invalid_datasets_df
            ], ignore_index=True)

        # Note: We include ALL datasets in the main results file regardless of threshold
        # The threshold will be applied only to the batch_analysis_input.csv file
        logging.info(f"Main results will include all {len(final_results)} datasets (threshold applied only to batch analysis file)")

        # Remove embedding column to prevent CSV malformation
        if 'embedding' in final_results.columns:
            final_results = final_results.drop('embedding', axis=1)

        # Create simplified final dataframe with only requested columns
        final = pd.DataFrame()

        # Helper functions for URLs
        def create_pubmed_url(pmid):
            if pd.notna(pmid):
                return f"https://pubmed.ncbi.nlm.nih.gov/{int(pmid)}/"
            return None

        def create_geo_url(accession):
            if pd.notna(accession):
                return f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}"
            return None

        # Create the simplified output columns
        final['GEO_Accession'] = final_results['ID']
        final['Title'] = final_results.get('Title', '')
        final['Summary'] = final_results.get('Summary', '')
        final['Species'] = final_results.get('Species', 'Unknown')
        final['Date'] = final_results.get('Date', '')
        final['Samples'] = final_results.get('rnaseq_samples', 0)  # Number of valid RNA-seq samples
        final['Valid'] = final_results.get('valid_dataset', False).apply(lambda x: 'Yes' if x else 'No')
        final['Validation_Result'] = final_results.get('validation_reason', 'Unknown')
        final['DatasetSizeGB'] = final_results.get('DatasetSizeGB', 0.0)
        final['RelevanceScore'] = final_results.get('RelevanceScore', None)
        final['Justification'] = final_results.get('Run1Justification', None)
        final['PubMed_URL'] = final_results.get('PrimaryPubMedID').apply(create_pubmed_url)
        final['GEO_URL'] = final_results['ID'].apply(create_geo_url)

        final = final.drop_duplicates(subset=['GEO_Accession'])

        # Save to output directory
        main_output_file = output_path / "Dataset_identification_result.csv"

        final.to_csv(main_output_file, index=False)
        message = f'Saved combined table to {main_output_file}'
        logging.info(message)

        # Generate multi-dataset CSV for batch analysis (default behavior)
        # This CSV only includes datasets that are valid, assessed, and above threshold
        multi_df = final[['GEO_Accession', 'Species', 'RelevanceScore', 'Valid', 'DatasetSizeGB']].copy()

        # Filter for batch analysis CSV (more restrictive than main results)
        valid_datasets = multi_df[multi_df['Valid'] == 'Yes']
        assessed_datasets = valid_datasets.dropna(subset=['RelevanceScore'])  # Only assessed datasets
        above_threshold = assessed_datasets[assessed_datasets['RelevanceScore'] >= args.threshold]  # Apply threshold

        multi_df = above_threshold.sort_values('RelevanceScore', ascending=False)
        multi_df = multi_df.rename(columns={'Species': 'organism', 'GEO_Accession': 'Accession'})

        # Log filtering steps for transparency
        logging.info(f"Batch analysis filtering: {len(valid_datasets)} valid → {len(assessed_datasets)} assessed → {len(above_threshold)} above threshold ({args.threshold})")

        # Set path for the multi-dataset CSV in the output directory
        multi_csv_path = output_path / 'selected_datasets.csv'

        multi_df.to_csv(multi_csv_path, index=False)
        message = f"Generated batch analysis input CSV at {multi_csv_path} with {len(multi_df)} datasets (threshold {args.threshold} applied)"
        logging.info(message)

        # Log summary results
        logging.info("="*80)
        logging.info("DATASET IDENTIFICATION RESULTS")
        logging.info("="*80)
        logging.info(f"Research Query: {research_query}")
        logging.info(f"Total datasets found: {len(final)}")
        logging.info(f"Valid datasets: {len(final[final['Valid'] == 'Yes'])}")
        logging.info(f"Invalid datasets: {len(final[final['Valid'] == 'No'])}")
        assessed_count = len(final.dropna(subset=['RelevanceScore']))
        logging.info(f"Assessed datasets: {assessed_count}")
        not_assessed_count = len(final[final['Valid'] == 'Yes']) - assessed_count
        if not_assessed_count > 0:
            logging.info(f"Valid but not assessed: {not_assessed_count} (not selected as cluster representatives)")

        if args.threshold > 0:
            above_threshold_count = len(final[final['RelevanceScore'].notna() & (final['RelevanceScore'] >= args.threshold)])
            logging.info(f"Above threshold ({args.threshold}): {above_threshold_count}")
            logging.info(f"Included in batch analysis CSV: {len(multi_df)}")

        assessed_final = final.dropna(subset=['RelevanceScore']).sort_values('RelevanceScore', ascending=False)

        # Generate metadata JSON file
        end_time = datetime.datetime.now()

        # Calculate metadata
        total_datasets_assessed = len(final.dropna(subset=['RelevanceScore'])) if 'RelevanceScore' in final.columns else 0
        datasets_deemed_relevant = len(final[final['RelevanceScore'].notna() & (final['RelevanceScore'] >= args.threshold)]) if 'RelevanceScore' in final.columns and args.threshold > 0 else 0

        metadata = {
            "input_query": research_query,
            "start_time": start_time.isoformat(),
            "end_time": end_time.isoformat(),
            "total_datasets_assessed": total_datasets_assessed,
            "datasets_deemed_relevant": datasets_deemed_relevant,
            "threshold_used": args.threshold
        }

        # Save metadata JSON
        metadata_file = output_path / "identification_metadata.json"
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)

        logging.info(f"Saved identification metadata to: {metadata_file}")

    except Exception as e:
        logging.error(f"Error in main execution: {e}")
        raise

if __name__ == "__main__":
    main()

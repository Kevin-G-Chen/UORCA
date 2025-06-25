#!/usr/bin/env python3
"""
Optimized Dataset Identification Script
======================================

This script performs dataset identification with early filtering to improve efficiency.
Key optimizations:
- Early RNA-seq filtering before expensive operations
- Species consistency validation
- Only fetch summaries for valid datasets
- Enhanced reporting with filtering details

Usage:
    python OptimizedDatasetIdentification.py "research query" --max-datasets 10
"""

import argparse
import asyncio
import json
import logging
import os
import statistics
import sys
import threading
import time
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

    # Configure logging to both file and console
    logging.basicConfig(
        level=level,
        format=format_str,
        handlers=[
            logging.FileHandler(logs_dir / "dataset_identification.log"),
            logging.StreamHandler(sys.stdout)
        ]
    )

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

def get_basic_dataset_info(geo_ids: List[str], api_delay: float = 0.4) -> pd.DataFrame:
    """Get basic dataset information without full summaries."""
    datasets = []
    unique_ids = list(dict.fromkeys(geo_ids))  # Ensure unique IDs to avoid redundant requests

    logging.info(f"Fetching basic dataset information for {len(unique_ids)} unique datasets...")

    # Split into batches of 10 for efficient fetching and to follow NCBI guidelines
    batches = [unique_ids[i:i+10] for i in range(0, len(unique_ids), 10)]
    logging.info(f"Processing {len(batches)} batches of GEO summaries")

    for batch in batches:
        try:
            handle = Entrez.esummary(db="gds", id=','.join(batch), retmode="xml")
            summaries = Entrez.read(handle)
            handle.close()

            for i, rec in enumerate(summaries):
                # Build basic dataset info using same approach as original script
                accession = rec.get('Accession', '')

                # Only include datasets with valid accession
                if accession:
                    datasets.append({
                        "ID": accession,  # Use accession directly as ID like original script
                        "GEO_ID": batch[i] if i < len(batch) else batch[0],  # Map to corresponding batch ID
                        "Title": rec.get('title', ''),
                        "Accession": accession
                    })

            time.sleep(api_delay)

        except Exception as e:
            logging.warning(f"Error processing batch {batch}: {e}")
            continue

    # Remove duplicates by ID (accession)
    df = pd.DataFrame(datasets)
    if not df.empty:
        df = df.drop_duplicates(subset=['ID']).reset_index(drop=True)

    logging.info(f"Found {len(df)} unique GSE datasets")
    return df

# XML and SRA processing functions
def strip_ns(root: ET.Element) -> None:
    """Remove XML namespaces for easier XPath."""
    for el in root.iter():
        if isinstance(el.tag, str) and "}" in el.tag:
            el.tag = el.tag.split("}", 1)[1]

def get_geo_uid(gse: str) -> str:
    """Get GEO UID from GSE accession."""
    handle = Entrez.esearch(db="gds", term=gse)
    result = Entrez.read(handle)
    handle.close()
    if result["IdList"]:
        return result["IdList"][0]
    raise ValueError(f"No GEO UID found for {gse}")

def get_sra_uid(geo_uid: str) -> str:
    """Get SRA UID from GEO UID."""
    handle = Entrez.elink(dbfrom="gds", db="sra", id=geo_uid)
    result = Entrez.read(handle)
    handle.close()

    if result and result[0].get("LinkSetDb"):
        links = result[0]["LinkSetDb"]
        for link_set in links:
            if link_set.get("Link"):
                return link_set["Link"][0]["Id"]

    raise ValueError(f"No SRA UID found for GEO UID {geo_uid}")

def fetch_sra_xml(sra_uid: str) -> str:
    """Fetch SRA XML data."""
    handle = Entrez.efetch(db="sra", id=sra_uid, rettype="xml", retmode="text")
    xml_data = handle.read()
    handle.close()
    return xml_data

def parse_bioprojects(xml_str: str) -> List[str]:
    """Parse BioProject accessions from SRA XML."""
    root = ET.fromstring(xml_str)
    strip_ns(root)
    nodes = root.findall(".//STUDY") + root.findall(".//PROJECT")
    accs = {n.attrib["accession"] for n in nodes if "accession" in n.attrib}
    if not accs:
        raise ValueError("No BioProject accession found in SRA XML")
    return sorted(accs)

def fetch_runinfo(gse: str, api_delay: float = 0.4, retmax: int = 100_000) -> pd.DataFrame:
    """Fetch RunInfo data for a GSE accession with no retries for faster processing."""
    try:
        # GSE -> GEO UID
        geo_uid = get_geo_uid(gse)
        time.sleep(api_delay)

        # GEO UID -> SRA UID
        sra_uid = get_sra_uid(geo_uid)
        time.sleep(api_delay)

        # SRA UID -> XML -> BioProjects
        xml_data = fetch_sra_xml(sra_uid)
        bioprojects = parse_bioprojects(xml_data)
        time.sleep(api_delay)

        # BioProjects -> RunInfo
        frames = []
        for bp in bioprojects:
            handle = Entrez.esearch(db="sra", term=f"{bp}[BioProject]", retmax=retmax)
            uid_list = Entrez.read(handle)["IdList"]
            handle.close()

            for uid in uid_list:
                handle = Entrez.efetch(db="sra", id=uid, rettype="runinfo", retmode="text")
                try:
                    df = pd.read_csv(handle, low_memory=False)
                    frames.append(df)
                except Exception as e:
                    logging.warning(f"Error reading RunInfo for {uid}: {e}")
                handle.close()
                time.sleep(api_delay)

        return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

    except Exception as e:
        logging.warning(f"Error fetching RunInfo for {gse}: {e}")
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


def validate_and_size_datasets(datasets_df: pd.DataFrame, api_delay: float = 0.4, max_workers: int = None) -> pd.DataFrame:
    """Validate datasets and calculate sizes with RNA-seq filtering using parallel processing."""

    # Auto-determine optimal worker count based on API key if not specified
    if max_workers is None:
        if os.getenv("ENTREZ_API_KEY"):
            max_workers = 8  # More aggressive with API key
            logging.info("🚀 API key detected - using 8 validation workers")
        else:
            max_workers = 3  # Conservative without API key
            logging.info("⚠️ No API key - using 3 validation workers")

    # Update global rate limiter
    global api_rate_limiter
    api_rate_limiter = APIRateLimiter()

    def process_dataset_with_delay(gse_id, api_delay):
        """Process a single dataset with API delay and rate limiting."""
        try:
            logging.info(f"Processing {gse_id}...")
            # Apply rate limiting before making API calls
            api_rate_limiter.wait()
            runinfo_df = fetch_runinfo(gse_id, api_delay)

            if runinfo_df.empty:
                return {
                    'ID': gse_id,
                    'original_samples': 0,
                    'filtered_samples': 0,
                    'species_count': 0,
                    'species_list': [],
                    'valid_dataset': False,
                    'skip_reason': 'No runinfo data available'
                }

            original_samples = len(runinfo_df)

            # Apply RNA-seq filters
            rna_seq_filter = (
                (runinfo_df['LibrarySource'] == 'TRANSCRIPTOMIC') &
                (runinfo_df['LibraryStrategy'] == 'RNA-Seq')
            )
            filtered_df = runinfo_df[rna_seq_filter]
            filtered_samples = len(filtered_df)

            # Check if any samples remain after filtering
            if filtered_samples == 0:
                return {
                    'ID': gse_id,
                    'original_samples': original_samples,
                    'filtered_samples': 0,
                    'species_count': 0,
                    'species_list': [],
                    'valid_dataset': False,
                    'skip_reason': 'No RNA-seq samples after filtering'
                }

            # Check species consistency
            unique_species = filtered_df['TaxID'].nunique()
            species_list = filtered_df['TaxID'].unique().tolist()

            if unique_species > 1:
                return {
                    'ID': gse_id,
                    'original_samples': original_samples,
                    'filtered_samples': filtered_samples,
                    'species_count': unique_species,
                    'species_list': species_list,
                    'valid_dataset': False,
                    'skip_reason': f'Multiple species detected ({unique_species})'
                }

            # Check minimum sample count
            if filtered_samples < 3:
                return {
                    'ID': gse_id,
                    'original_samples': original_samples,
                    'filtered_samples': filtered_samples,
                    'species_count': unique_species,
                    'species_list': species_list,
                    'valid_dataset': False,
                    'skip_reason': f'Insufficient samples ({filtered_samples} < 3)'
                }

            return {
                'ID': gse_id,
                'original_samples': original_samples,
                'filtered_samples': filtered_samples,
                'species_count': unique_species,
                'species_list': species_list,
                'valid_dataset': True,
                'skip_reason': None
            }

        except Exception as e:
            logging.warning(f"Error processing {gse_id}: {e}")
            return {
                'ID': gse_id,
                'original_samples': 0,
                'filtered_samples': 0,
                'species_count': 0,
                'species_list': [],
                'valid_dataset': False,
                'skip_reason': f'Error: {str(e)}'
            }

    # Process datasets in parallel
    logging.info(f"Validating datasets with RNA-seq filtering using {max_workers} parallel workers...")
    gse_ids = datasets_df['ID'].tolist()
    validation_data = []

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Create partial function with fixed api_delay
        process_func = partial(process_dataset_with_delay, api_delay=api_delay)

        # Submit all tasks
        future_to_gse = {executor.submit(process_func, gse_id): gse_id for gse_id in gse_ids}

        # Collect results as they complete
        for future in as_completed(future_to_gse):
            validation_info = future.result()
            validation_data.append(validation_info)

            # Progress logging
            completed = len(validation_data)
            total = len(gse_ids)
            if completed % 10 == 0 or completed == total:
                logging.info(f"Validation progress: {completed}/{total} datasets processed")

    # Create results dataframe
    validation_df = pd.DataFrame(validation_data)

    # Merge with original dataset info
    result_df = datasets_df.merge(validation_df, on='ID', how='left')

    # Log summary
    total_datasets = len(result_df)
    valid_datasets = result_df['valid_dataset'].sum()
    logging.info(f"Dataset validation summary: {valid_datasets}/{total_datasets} datasets are valid")

    return result_df

def fetch_geo_summaries_for_valid(valid_datasets_df: pd.DataFrame, api_delay: float = 0.4) -> pd.DataFrame:
    """Fetch detailed GEO summaries only for valid datasets."""
    logging.info(f"Fetching GEO summaries for {len(valid_datasets_df)} valid datasets...")

    enriched_datasets = []

    for _, row in valid_datasets_df.iterrows():
        gse_id = row['ID']
        geo_id = row['GEO_ID']

        try:
            # Fetch detailed summary
            handle = Entrez.esummary(db="gds", id=geo_id)
            summary = Entrez.read(handle)[0]
            handle.close()

            # Extract relevant information
            title = summary.get("title", "")
            summary_text = summary.get("summary", "")
            organism = summary.get("taxon", "")

            enriched_datasets.append({
                **row.to_dict(),
                "Title": title,
                "Summary": summary_text,
                "Species": organism,
                "Technique": "RNA-seq"  # We know this from filtering
            })

            time.sleep(api_delay)

        except Exception as e:
            logging.warning(f"Error fetching summary for {gse_id}: {e}")
            # Keep the dataset but with limited info
            enriched_datasets.append({
                **row.to_dict(),
                "Summary": "Summary not available",
                "Species": "Unknown",
                "Technique": "RNA-seq"
            })

    return pd.DataFrame(enriched_datasets)

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
    advanced_group.add_argument('--validation-workers', type=int, default=None,
                               help='Number of parallel workers for dataset validation (auto-detects based on API key)')
    advanced_group.add_argument('--sra-workers', type=int, default=None,
                               help='Number of parallel workers for SRA data fetching (auto-detects based on API key)')
    advanced_group.add_argument('--min-samples', type=int, default=None,
                                help='Minimum number of samples in a neighbourhood for HDBSCAN (None lets HDBSCAN default to min_cluster_size)')

    args = parser.parse_args()

    # Set up logging in the output directory
    from pathlib import Path
    output_path = Path(args.output)
    # Ensure output directory exists
    output_path.mkdir(parents=True, exist_ok=True)
    setup_logging(verbose=False)  # Use consistent logging setup

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
        basic_datasets_df = get_basic_dataset_info(unique_geo_ids, args.api_delay)

        if basic_datasets_df.empty:
            logging.error("No valid GSE datasets found")
            return

        # Step 4: EARLY VALIDATION - Filter datasets with RNA-seq criteria
        logging.info("Performing early validation with RNA-seq filtering...")
        validated_datasets_df = validate_and_size_datasets(basic_datasets_df, args.api_delay, max_workers=args.validation_workers)

        # Split into valid and invalid datasets
        valid_datasets_df = validated_datasets_df[validated_datasets_df['valid_dataset'] == True].copy()
        invalid_datasets_df = validated_datasets_df[validated_datasets_df['valid_dataset'] == False].copy()

        logging.info(f"Valid datasets: {len(valid_datasets_df)}")
        logging.info(f"Invalid datasets: {len(invalid_datasets_df)}")

        if len(valid_datasets_df) == 0:
            logging.error("No valid datasets found after filtering")

            # Still show the invalid datasets for transparency
            print("\n" + "="*80)
            print("DATASET IDENTIFICATION RESULTS")
            print("="*80)
            print(f"Research Query: {research_query}")
            print(f"No valid datasets found after RNA-seq filtering")
            print(f"\nInvalid datasets ({len(invalid_datasets_df)}):")

            for _, row in invalid_datasets_df.iterrows():
                print(f"\nDataset: {row['ID']}")
                print(f"  Original samples: {row['original_samples']}")
                print(f"  RNA-seq samples: {row['filtered_samples']}")
                print(f"  Species count: {row['species_count']}")
                print(f"  Valid for analysis: ❌ No")
                print(f"  Skip reason: {row['skip_reason']}")

            return

        # Step 5: Fetch detailed summaries only for valid datasets
        enriched_valid_df = fetch_geo_summaries_for_valid(valid_datasets_df, args.api_delay)

        # Step 6: Embed and cluster only valid datasets
        logging.info("Embedding valid datasets...")
        embedded_df = embed_datasets(enriched_valid_df)

        logging.info("Clustering valid datasets...")
        clustered_df = cluster_datasets(embedded_df)

        # Step 7: Select representatives
        logging.info("Selecting representative datasets...")
        representatives_df = select_representative_datasets(clustered_df, args.max_evaluate or 10)

        # Step 8: Assess relevance with repeated scoring
        # ==================== REPLACE WITH SIMPLE ASSESSMENT IF DESIRED ====================
        # To use single assessment instead of repeated assessment, replace the line below with:
        # assessed_df = asyncio.run(assess_single_relevance(representatives_df, research_query))
        # (and implement assess_single_relevance function similar to the original assess_relevance)
        logging.info("Assessing relevance with repeated scoring...")
        assessed_df = asyncio.run(repeated_relevance(representatives_df, research_query))
        # ===================================================================================

        # ==================== OPTIONAL: SRA DATA FETCHING SECTION ====================
        # The following section fetches SRA metadata and calculates dataset sizes
        # Remove this entire section if you don't need dataset sizes or SRA metadata
        to_fetch = assessed_df[assessed_df['RelevanceScore'] >= args.threshold]
        logging.info(f"🔗 Fetching SRA metadata for {len(to_fetch)} datasets (threshold >= {args.threshold})...")

        def fetch_sra_for_dataset(row, api_delay):
            """Fetch SRA data for a single dataset with rate limiting."""
            acc = row['ID']
            try:
                # Apply rate limiting before making API calls
                api_rate_limiter.wait()
                g_uid = get_geo_uid(acc)
                api_rate_limiter.wait()
                s_uid = get_sra_uid(g_uid)
                api_rate_limiter.wait()
                xml = fetch_sra_xml(s_uid)
                bioprojects = parse_bioprojects(xml)

                # Handle multiple BioProjects for a single GEO accession
                acc_frames = []
                for bp in bioprojects:
                    df_run = fetch_runinfo(bp, api_delay)
                    if not df_run.empty:
                        df_run.insert(0, 'GEO_Accession', acc)
                        acc_frames.append(df_run)

                if acc_frames:
                    combined_df = pd.concat(acc_frames, ignore_index=True)
                    logging.info(f"Fetched {len(combined_df)} runs for {acc} from {len(bioprojects)} BioProject(s)")
                    return combined_df, True
                else:
                    logging.warning(f"No runs found for {acc}")
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
        # ===============================================================================

        # Merge assessed datasets with original data to get all required columns
        geo_full = representatives_df.merge(assessed_df, on='ID', how='left')

        # ==================== OPTIONAL: DATASET SIZE AND SRA MERGE ====================
        # The following section adds dataset sizes and merges with SRA data
        # If SRA fetching section above is removed, replace this with: merged = geo_full
        if dataset_sizes:
            geo_full['DatasetSizeBytes'] = geo_full['ID'].map(dataset_sizes)
            geo_full['DatasetSizeGB'] = geo_full['DatasetSizeBytes'] / (1024**3)
            # Fill NaN values with 0 for datasets where size couldn't be calculated
            geo_full['DatasetSizeBytes'] = geo_full['DatasetSizeBytes'].fillna(0)
            geo_full['DatasetSizeGB'] = geo_full['DatasetSizeGB'].fillna(0.0)
        else:
            geo_full['DatasetSizeBytes'] = 0
            geo_full['DatasetSizeGB'] = 0.0

        # Merge all GEO datasets with fetched SRA info
        merged = geo_full.merge(
            sra_df,
            left_on='ID',
            right_on='GEO_Accession',
            how='left'
        )
        logging.info(f"Merged GEO and SRA data for {len(merged)} total records")
        # ===============================================================================

        # Add invalid datasets with default scores
        invalid_for_output = invalid_datasets_df.copy()
        invalid_for_output['RelevanceScore'] = 0
        invalid_for_output['Justification'] = "Dataset excluded due to filtering criteria"
        invalid_for_output['Cluster'] = -1
        # ==================== OPTIONAL: Add columns for invalid datasets ====================
        # Add missing optional columns for invalid datasets (remove this section if columns above are removed)
        for col in ['Run1Score', 'Run2Score', 'Run3Score', 'Run1Justification', 'Run2Justification', 'Run3Justification']:
            invalid_for_output[col] = None
        invalid_for_output['DatasetSizeBytes'] = 0
        invalid_for_output['DatasetSizeGB'] = 0.0
        invalid_for_output['LibrarySource'] = None
        invalid_for_output['LibraryLayout'] = None
        invalid_for_output['LibraryStrategy'] = None
        # ===================================================================================

        # Combine valid (merged) and invalid datasets for final output
        final_results = pd.concat([
            merged,
            invalid_for_output
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
        final['NumSamples'] = final_results.get('filtered_samples', final_results.get('NumSamples', 0))
        final['PrimaryPubMedID'] = None  # We don't have PubMed IDs in this workflow
        final['AllPubMedIDs'] = None
        final['RelevanceScore'] = final_results['RelevanceScore']

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
        final['Valid'] = final.apply(
            lambda row: 'Yes' if (row.get('LibraryLayout') == 'PAIRED' and row.get('LibrarySource') == 'TRANSCRIPTOMIC') else 'No',
            axis=1
        )

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
        print(f"Total datasets found: {len(final)}")
        print(f"Valid datasets: {len(final[final['Valid'] == 'Yes'])}")
        print(f"Invalid datasets: {len(final[final['Valid'] == 'No'])}")

        print(f"\nTop {min(len(final), args.max_evaluate or 10)} Results:")
        print("-" * 80)

        for i, (_, row) in enumerate(final.head(args.max_evaluate or 10).iterrows()):
            print(f"\n{i+1}. Dataset: {row['Accession']}")
            print(f"   Title: {row['Title']}")
            print(f"   RNA-seq samples: {row['NumSamples']}")
            print(f"   Species: {row['Species']}")
            print(f"   Valid for analysis: {'✅ Yes' if row['Valid'] == 'Yes' else '❌ No'}")
            print(f"   Relevance: {row['RelevanceScore']}/10")
            print(f"   Dataset Size: {row['DatasetSizeGB']:.2f} GB")
            print(f"   GEO URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={row['Accession']}")

    except Exception as e:
        logging.error(f"Error in main execution: {e}")
        raise

if __name__ == "__main__":
    main()

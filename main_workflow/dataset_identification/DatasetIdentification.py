#!/usr/bin/env python3
"""
UORCA Dataset Identification - AI-powered discovery of relevant RNA-seq datasets

This script intelligently identifies RNA-seq datasets from GEO that are relevant to your biological
research question. It combines multiple AI techniques for comprehensive dataset discovery:

Key Features:
 • AI-powered relevance scoring using GPT-4 with multiple evaluation rounds
 • Optional dataset clustering for diverse, representative selection
 • Fast dataset size calculation using existing NCBI metadata
 • Comprehensive output with dataset validation and batch analysis preparation

Main Workflow:
 1. Extract biological terms from your research query and search GEO
 2. Score dataset relevance using AI with detailed biological justification
 3. Optionally cluster datasets by similarity to ensure diverse representation
 4. Fetch SRA metadata and calculate dataset sizes for selected datasets
 5. Generate comprehensive results with validity flags for downstream analysis

Usage Examples:
  # Basic usage - find neuroblastoma datasets
  python DatasetIdentification.py -q "neuroblastoma tumor vs normal tissue"

  # Use clustering for diverse selection with batch analysis output
  python DatasetIdentification.py -q "cardiac development mouse" --use-clustering --multi-dataset-csv

  # Lower threshold for broader discovery
  python DatasetIdentification.py -q "MYCN amplification" --threshold 6.0

Output Files:
 • Dataset_identification_result.csv: Complete results with relevance scores and metadata
 • batch_analysis_input.csv: Filtered datasets ready for UORCA batch analysis (optional)
 • dataset_clusters.csv: Clustering details (if --save-clustering-details used)

Technical Improvements:
 • Fast dataset size calculation using existing runinfo metadata (no subprocess calls)
 • Robust fetching with retmax=100,000 to get all UIDs in one call
 • Support for multiple BioProjects per GEO series (common in complex studies)
 • Improved retry logic and error handling for reliable NCBI API interactions
 • Optimized performance: 10-100x faster than previous vdb-dump approach
"""
from __future__ import annotations
import argparse
import asyncio
import datetime
import json
import os
import statistics
import sys
import time
import logging
from pathlib import Path
from typing import Any, Dict, List, Tuple, Optional
from collections import defaultdict

import pandas as pd
import numpy as np
from pandas._libs.algos import Infinity
from Bio import Entrez
from dotenv import load_dotenv
from openai import OpenAI
from pydantic import BaseModel, ConfigDict
from tqdm import tqdm
import xml.etree.ElementTree as ET
import hdbscan
from sklearn.preprocessing import normalize



# ----------------------
# Environment & prompts
# ----------------------
load_dotenv()
Entrez.email = os.getenv("ENTREZ_EMAIL")
Entrez.api_key = os.getenv("ENTREZ_API_KEY")
openai_api_key = os.getenv("OPENAI_API_KEY")
if not openai_api_key:
    raise RuntimeError("OPENAI_API_KEY not set")
client = OpenAI(api_key=openai_api_key)
PROMPT_DIR = os.getenv("PROMPT_DIR", "./main_workflow/prompts/dataset_identification")

# Set up logging
logger = logging.getLogger(__name__)

def setup_logging(output_dir=None):
    """Configure logging to both console and file if output_dir is provided."""
    log_format = "%(asctime)s - %(levelname)s - %(message)s"
    handlers = [logging.StreamHandler()]

    # Add file handler if output directory is provided
    if output_dir:
        log_dir = Path(output_dir)
        log_dir.mkdir(parents=True, exist_ok=True)
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        log_file = log_dir / f"dataset_identification_{timestamp}.log"
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(logging.Formatter(log_format))
        handlers.append(file_handler)
        print(f"Logging to file: {log_file}")

    # Configure root logger
    logging.basicConfig(
        level=logging.WARNING,
        format=log_format,
        handlers=handlers,
        force=True
    )

    return logger

def load_prompt(fname: str) -> str:
    path = Path(PROMPT_DIR) / fname
    return path.read_text().strip()

# ----------------------
# Pydantic models
# ----------------------
class ExtractedTerms(BaseModel):
    model_config = ConfigDict(extra="forbid")
    extracted_terms: List[str]
    expanded_terms: List[str]

class Assessment(BaseModel):
    model_config = ConfigDict(extra="forbid")
    ID: str
    RelevanceScore: float
    Justification: str

class Assessments(BaseModel):
    model_config = ConfigDict(extra="forbid")
    assessments: List[Assessment]

TERMS_SCHEMA = ExtractedTerms.model_json_schema()
ASSESS_SCHEMA = Assessments.model_json_schema()

# ----------------------
# OpenAI JSON helper
# ----------------------
def call_openai_json(prompt: str, schema: Dict[str, Any], name: str) -> dict:
    resp = client.responses.create(
        model="gpt-4.1-mini",
        input=prompt,
        text={"format": {"type": "json_schema", "name": name, "schema": schema, "strict": True}},
    )
    return json.loads(resp.output_text)

# ----------------------
# Step 1 & 2: GEO search
# ----------------------
def extract_terms(query: str) -> ExtractedTerms:
    sys_prompt = load_prompt("extract_terms.txt")
    prompt = f"{sys_prompt}\nQuery: {query}"
    data = call_openai_json(prompt, TERMS_SCHEMA, "extracted_terms")
    return ExtractedTerms.model_validate(data)

def perform_search(term: str, retmax: int = 10000) -> List[str]:
    search = f"{term} AND (\"Expression profiling by high throughput sequencing\"[Filter])"
    with Entrez.esearch(db="gds", term=search, retmode="xml", retmax=retmax) as h:
        out = Entrez.read(h)
    return out.get("IdList", [])


def fetch_geo_summaries(ids: List[str]) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    unique_ids = list(dict.fromkeys(ids))  # Ensure unique IDs to avoid redundant requests
    print(f"⬇️ Fetching GEO summaries for {len(unique_ids)} unique datasets...")

    # Split into batches of 10 for efficient fetching and to follow NCBI guidelines
    batches = [unique_ids[i:i+10] for i in range(0, len(unique_ids), 10)]
    logger.info(f"Processing {len(batches)} batches of GEO summaries")

    for batch in tqdm(batches, desc="GEO summaries"):
        with Entrez.esummary(db="gds", id=','.join(batch), retmode="xml") as h:
            summaries = Entrez.read(h)
        for rec in summaries:
            # Safely extract PubMed IDs if available
            pmids = []
            try:
                if 'PubMedIds' in rec and rec['PubMedIds']:
                    pmids = [int(pmid) for pmid in rec['PubMedIds']]
            except (ValueError, TypeError) as e:
                logger.warning(f"Error extracting PubMed IDs for {rec.get('Accession', '')}: {e}")

            # Build row with optional PubMed data
            row_data = {
                'ID': rec.get('Accession', ''),
                'Title': rec.get('title', ''),
                'Summary': rec.get('summary', ''),
                'Accession': rec.get('Accession', ''),
                'Species': rec.get('taxon', ''),
                'Date': rec.get('PDAT', '')
            }

            # Add PubMed information when available
            if pmids:
                row_data['PrimaryPubMedID'] = pmids[0]
                row_data['AllPubMedIDs'] = ','.join(str(pid) for pid in pmids)
            else:
                row_data['PrimaryPubMedID'] = None
                row_data['AllPubMedIDs'] = None

            rows.append(row_data)

        # NCBI recommends no more than 3 requests per second without an API key
        # With an API key, 10 requests per second are allowed
        # We'll use 0.4s sleep (conservative for ~2.5 requests/second)
        time.sleep(0.4)

    result_df = pd.DataFrame(rows)
    logger.info(f"Successfully fetched {len(result_df)} GEO summaries")
    return result_df

# ----------------------
# Step 2.5: Dataset Clustering
# ----------------------
def embed_datasets(df: pd.DataFrame, max_workers: int = 4) -> np.ndarray:
    """
    Create embeddings for each dataset based on its title and summary.
    Uses parallel processing to speed up embedding generation.

    Args:
        df: DataFrame containing GEO datasets with 'Title' and 'Summary' columns
        max_workers: Maximum number of parallel workers for embedding

    Returns:
        numpy array of embeddings with shape (n_datasets, embedding_dim)
    """
    from concurrent.futures import ThreadPoolExecutor, as_completed

    texts = []
    for _, row in df.iterrows():
        # Combine title and summary for a richer embedding
        text = f"{row['Title']}. {row['Summary']}"
        if pd.notna(text) and len(text.strip()) > 0:
            texts.append(text)
        else:
            texts.append("No information available")

    def get_embedding(text, idx):
        """Get embedding for a single text"""
        try:
            response = client.embeddings.create(
                model="text-embedding-3-small",
                input=text
            )
            embedding = response.data[0].embedding
            return idx, embedding
        except Exception as e:
            logger.warning(f"Error creating embedding at index {idx}: {e}")
            # Add a zero vector as fallback
            return idx, [0.0] * 1536  # Default dimension for embedding models

    # Initialize an array to store embeddings in the correct order
    embeddings = [None] * len(texts)

    # Use ThreadPoolExecutor for parallel processing
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all embedding tasks
        future_to_idx = {
            executor.submit(get_embedding, text, idx): idx
            for idx, text in enumerate(texts)
        }

        # Process results as they complete
        for future in tqdm(as_completed(future_to_idx), total=len(texts), desc=f"Creating embeddings (using {max_workers} workers)"):
            idx, embedding = future.result()
            embeddings[idx] = embedding

    # Stack embeddings into a matrix
    embedding_matrix = np.vstack(embeddings)

    # Normalize embeddings for better clustering
    normalized_embeddings = normalize(embedding_matrix)

    return normalized_embeddings

def cluster_datasets(embeddings: np.ndarray, min_cluster_size: int = 2, min_samples: int = None) -> Tuple[List[int], Dict[int, List[int]]]:
    """
    Cluster dataset embeddings using HDBSCAN.

    Args:
        embeddings: Matrix of dataset embeddings
        min_cluster_size: Minimum size of clusters (replaces eps parameter)
        min_samples: Optional minimum number of samples required for a core point

    Returns:
        Tuple of (cluster labels, dict mapping cluster labels to indices)
    """
    # Apply HDBSCAN clustering
    clusterer = hdbscan.HDBSCAN(
        min_cluster_size=min_cluster_size,
        min_samples=min_samples if min_samples is not None else min_cluster_size,
        metric="euclidean"
    )
    labels = clusterer.fit_predict(embeddings)

    # Group indices by cluster
    clusters = defaultdict(list)
    for idx, label in enumerate(labels):
        clusters[int(label)].append(idx)

    return labels.tolist(), clusters

def select_representative_datasets(df: pd.DataFrame, clusters: Dict[int, List[int]],
                                  max_per_cluster: int = 3) -> pd.DataFrame:
    """
    Select representative datasets from each cluster for relevance assessment.

    Args:
        df: DataFrame containing GEO datasets
        clusters: Dict mapping cluster labels to lists of indices
        max_per_cluster: Maximum number of datasets to select from each cluster

    Returns:
        DataFrame with selected representative datasets
    """
    selected_indices = []

    # For each cluster, select up to max_per_cluster representatives
    for cluster_label, indices in clusters.items():
        # Always include at least one dataset from each cluster
        # For larger clusters, select more representatives up to max_per_cluster
        n_select = min(max_per_cluster, len(indices))

        # For now, just take the first n_select items
        # In a more sophisticated implementation, we could select the most central items
        selected_indices.extend(indices[:n_select])

    # Get the selected datasets
    selected_df = df.iloc[selected_indices].copy()

    logger.info(f"Selected {len(selected_df)} representative datasets from {len(clusters)} clusters")

    return selected_df

# ----------------------
# Step 3: Relevance scoring
# ----------------------
async def assess_subbatch(sub_df: pd.DataFrame, query: str, schema: Dict[str, Any], name: str, rep: int, batch_idx: int, total_batches: int, sem: asyncio.Semaphore) -> List[Assessment]:
    async with sem:
        print(f"⏳ Relevance rep {rep+1}, batch {batch_idx+1}/{total_batches}")
        sys_prompt = load_prompt("assess_relevance.txt")
        prompt = f"{sys_prompt}\nResearch query: \"{query}\"\nDatasets to assess (batch {batch_idx+1}):\n" + sub_df.to_json(orient='records')
        data = await asyncio.to_thread(call_openai_json, prompt, schema, name)
        return [Assessment.model_validate(a) for a in data['assessments']]

async def repeated_relevance(df: pd.DataFrame, query: str, repeats: int, batch_size: int, openai_api_jobs: int) -> pd.DataFrame:
    print(f"📊 Starting relevance scoring: {repeats} repetitions, batch size {batch_size}, parallel API jobs: {openai_api_jobs}")
    sem = asyncio.Semaphore(openai_api_jobs)
    tasks = []
    # create sub-batches
    batches = [df.iloc[i:i+batch_size] for i in range(0, len(df), batch_size)]
    total_batches = len(batches)
    for rep in range(repeats):
        for idx, sub in enumerate(batches):
            tasks.append(assess_subbatch(sub, query, ASSESS_SCHEMA, "assessments", rep, idx, total_batches, sem))
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

# ----------------------
# SRA extraction utils
# ----------------------
def strip_ns(root: ET.Element):
    for e in root.iter():
        if isinstance(e.tag, str) and '}' in e.tag:
            e.tag = e.tag.split('}', 1)[1]


def get_geo_uid(acc: str) -> str:
    with Entrez.esearch(db="gds", term=acc) as h:
        rec = Entrez.read(h)
    # Apply rate limiting to follow NCBI guidelines
    time.sleep(args.api_delay)
    return rec['IdList'][0]


def get_sra_uid(geo_uid: str) -> str:
    max_retries = 3
    for attempt in range(max_retries):
        try:
            with Entrez.elink(dbfrom="gds", db="sra", id=geo_uid) as h:
                rec = Entrez.read(h)
            # Apply rate limiting to follow NCBI guidelines
            time.sleep(args.api_delay)
            return rec[0]['LinkSetDb'][0]['Link'][0]['Id']
        except Exception as e:
            if attempt < max_retries - 1:
                logger.warning(f"Failed to get SRA UID for {geo_uid}, attempt {attempt + 1}/{max_retries}: {e}")
                time.sleep(10)  # Wait 10 seconds before retry
            else:
                raise


def fetch_sra_xml(uid: str) -> str:
    max_retries = 3
    for attempt in range(max_retries):
        try:
            with Entrez.efetch(db="sra", id=uid, rettype="xml", retmode="text") as h:
                xml_data = h.read()
            # Apply rate limiting to follow NCBI guidelines
            time.sleep(args.api_delay)
            return xml_data
        except Exception as e:
            if attempt < max_retries - 1:
                logger.warning(f"Failed to fetch SRA XML for {uid}, attempt {attempt + 1}/{max_retries}: {e}")
                time.sleep(10)  # Wait 10 seconds before retry
            else:
                raise


def parse_bioprojects(xml_str: str) -> List[str]:
    """
    Return ALL BioProject accessions found in the SRA XML, because a single
    GEO series can map to multiple studies.
    """
    root = ET.fromstring(xml_str)
    strip_ns(root)
    nodes = root.findall('.//STUDY') + root.findall('.//PROJECT')
    bioprojects = sorted({
        n.attrib['accession'] for n in nodes if 'accession' in n.attrib
    })
    if not bioprojects:
        raise ValueError('No BioProject accession found')
    return bioprojects


def fetch_runinfo(proj_acc: str, api_delay: float = 0.4, retmax: int = 100_000) -> pd.DataFrame:
    """
    BioProject → DataFrame containing *all* runs in that project.

    Uses retmax=100_000 to grab every UID in one call (safe upper bound for
    almost all projects). Concatenates the RunInfo CSV from every UID.
    """
    max_retries = 3
    last_err: Optional[Exception] = None

    for attempt in range(max_retries):
        try:
            # 1. Find *all* SRA UIDs linked to this BioProject
            with Entrez.esearch(
                db="sra",
                term=f"{proj_acc}[BioProject]",
                retmax=retmax
            ) as h:
                rec = Entrez.read(h)
            time.sleep(api_delay)

            uid_list = rec["IdList"]
            if not uid_list:
                logger.warning(f"No UIDs found for BioProject {proj_acc}")
                return pd.DataFrame()  # empty result

            # 2. Fetch RunInfo for every UID and concatenate
            frames = []
            for uid in uid_list:
                with Entrez.efetch(db="sra", id=uid, rettype="runinfo", retmode="text") as h:
                    frames.append(pd.read_csv(h, low_memory=False))
                time.sleep(api_delay)

            df = pd.concat(frames, ignore_index=True)
            logger.info(f"Fetched {len(df)} SRA runs for {proj_acc} from {len(uid_list)} UIDs")
            return df

        except Exception as e:
            last_err = e
            if attempt < max_retries - 1:
                logger.warning(f"Failed to fetch runinfo for {proj_acc}, attempt {attempt + 1}/{max_retries}: {e}")
                time.sleep(10)  # Wait 10 seconds before retry
            else:
                raise last_err

def calculate_dataset_sizes_from_runinfo(sra_df: pd.DataFrame) -> Dict[str, int]:
    """
    Calculate dataset sizes using the size_MB column already present in runinfo data.

    This is much faster than calling vdb-dump for each SRR individually.

    Returns:
        Dictionary mapping GEO_Accession to total size in bytes
    """
    if sra_df.empty or 'size_MB' not in sra_df.columns or 'GEO_Accession' not in sra_df.columns:
        logger.warning("Cannot calculate dataset sizes: missing required columns in runinfo data")
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
    logger.info(f"Dataset size calculation complete: {total_datasets} datasets, total: {total_size_gb:.2f} GB")

    for geo_acc, size_bytes in dataset_sizes.items():
        logger.info(f"Dataset {geo_acc}: {size_bytes:,} bytes ({size_bytes/(1024**3):.2f} GB)")

    return dataset_sizes

# ----------------------
# Main
# ----------------------
def main():
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
                               help='Delay between NCBI API calls (seconds) to respect rate limits')
    advanced_group.add_argument('--min-samples', type=int, default=None,
                                help='Minimum number of samples in a neighbourhood for HDBSCAN (None lets HDBSCAN default to min_cluster_size)')


    global args
    args = parser.parse_args()

    # Set up logging in the output directory
    output_path = Path(args.output)
    # Ensure output directory exists
    output_path.mkdir(parents=True, exist_ok=True)
    setup_logging(output_path)

    message = "🧠 Extracting terms …"
    print(message)
    logger.info(message)

    terms = extract_terms(args.query)
    search_terms = set(terms.extracted_terms + terms.expanded_terms)
    log_message = f"Extracted {len(terms.extracted_terms)} direct terms and {len(terms.expanded_terms)} expanded terms"
    logger.info(log_message)
    log_message = f"Using {len(search_terms)} unique search terms: {', '.join(search_terms)}"
    logger.info(log_message)

    message = "🔍 Searching GEO …"
    print(message)
    logger.info(message)
    geo_ids: List[str] = []
    for term in search_terms:
        message = f"Searching for: {term}"
        print(message)
        logger.info(message)
        term_ids = perform_search(term, args.max_datasets)
        geo_ids.extend(term_ids)
        logger.info(f"Found {len(term_ids)} results for term: {term}")
        # Apply rate limiting between searches to follow NCBI guidelines
        time.sleep(args.api_delay)

    unique_geo_ids = list(dict.fromkeys(geo_ids))
    logger.info(f"Found {len(geo_ids)} total GEO IDs, {len(unique_geo_ids)} unique IDs")

    if not unique_geo_ids:
        message = 'No GEO IDs found'
        print(message)
        logger.info(message)
        sys.exit(0)

    geo_df = fetch_geo_summaries(unique_geo_ids)
    logger.info(f"Fetched metadata for {len(geo_df)} GEO datasets")

    # Store original dataframe before any filtering for later reference
    original_geo_df = geo_df.copy()

    # Log some basic dataset statistics
    species_counts = geo_df['Species'].value_counts()
    logger.info(f"Species distribution: {', '.join([f'{s}({c})' for s, c in species_counts.head(5).items()])}")
    if len(species_counts) > 5:
        logger.info(f"... and {len(species_counts)-5} more species")

    if args.use_clustering:
        message = "🔬 Clustering datasets by similarity..."
        print(message)
        logger.info(message)
        logger.info(f"Starting embedding and clustering for {len(geo_df)} datasets")

        # Create embeddings
        embeddings = embed_datasets(geo_df, max_workers=args.parallel_embedding_jobs)
        logger.info(f"Created embeddings with shape {embeddings.shape}")

        # Cluster the datasets
        labels, clusters = cluster_datasets(embeddings, min_cluster_size=args.min_cluster_size, min_samples=args.min_samples)

        # Count datasets per cluster
        cluster_counts = {}
        for label, indices in clusters.items():
            cluster_counts[label] = len(indices)

        # Log cluster distribution
        logger.info(f"HDBSCAN clustering results: {len(clusters)} clusters formed")
        for label, count in sorted(cluster_counts.items()):
            logger.info(f"Cluster {label}: {count} datasets")

        # Add cluster labels to the dataframe
        geo_df['ClusterLabel'] = labels

        # Select representative datasets from each cluster
        representative_df = select_representative_datasets(
            geo_df, clusters, max_per_cluster=args.max_per_cluster
        )
        logger.info(f"Selected {len(representative_df)} representative datasets from {len(clusters)} clusters")

        # Save clustering information if requested
        if args.save_clustering_details:
            cluster_info = geo_df[['ID', 'Accession', 'Title', 'ClusterLabel']].copy()
            cluster_path = Path(args.output) / 'dataset_clusters.csv'
            cluster_info.to_csv(cluster_path, index=False)
            print(f"💾 Saved clustering information to {cluster_path}")

        # Use the representative datasets for relevance assessment
        assess_df = representative_df
        message = f"📊 Selected {len(representative_df)} representative datasets for relevance assessment"
        print(message)
        logger.info(message)
    else:
        # Use all datasets (or max_evaluate limit) for relevance assessment
        if args.max_evaluate is not None:
            assess_df = geo_df.head(args.max_evaluate)
        else:
            assess_df = geo_df

    rel_df = asyncio.run(repeated_relevance(assess_df, args.query, args.scoring_rounds, args.batch_size, args.parallel_ai_jobs))

    # Merge relevance scores back to the original dataset
    if args.use_clustering:
        # Merge with the full dataset to include all datasets
        geo_full = original_geo_df.merge(rel_df, on='ID', how='left')
    else:
        geo_full = geo_df.merge(rel_df, on='ID', how='left')

    # Filter for SRA fetch
    to_fetch = geo_full[geo_full['RelevanceScore'] >= args.threshold]
    message = f"🔗 Fetching SRA metadata for {len(to_fetch)} datasets (threshold >= {args.threshold}) …"
    print(message)
    logger.info(message)
    logger.info(f"Selected {len(to_fetch)}/{len(geo_full)} datasets above relevance threshold {args.threshold}")

    runs: List[pd.DataFrame] = []
    successful_datasets = 0
    total_runs_fetched = 0

    for acc in tqdm(to_fetch['Accession'].unique(), desc='Processing SRA'):
        try:
            g_uid = get_geo_uid(acc)
            s_uid = get_sra_uid(g_uid)
            xml = fetch_sra_xml(s_uid)
            bioprojects = parse_bioprojects(xml)

            # Handle multiple BioProjects for a single GEO accession
            acc_frames = []
            for bp in bioprojects:
                df_run = fetch_runinfo(bp, args.api_delay)
                if not df_run.empty:
                    df_run.insert(0, 'GEO_Accession', acc)
                    acc_frames.append(df_run)

            if acc_frames:
                combined_df = pd.concat(acc_frames, ignore_index=True)
                runs.append(combined_df)
                successful_datasets += 1
                total_runs_fetched += len(combined_df)
                logger.info(f"Fetched {len(combined_df)} runs for {acc} from {len(bioprojects)} BioProject(s)")
            else:
                logger.warning(f"No runs found for {acc}")

        except Exception as e:
            logger.info(f'Warning, skipping {acc}: {e}')

        if runs:
            sra_df = pd.concat(runs, ignore_index=True)
            logger.info(f"Successfully processed {successful_datasets}/{len(to_fetch['Accession'].unique())} GEO datasets")
            logger.info(f"Retrieved metadata for {total_runs_fetched} total SRA runs")
        else:
            sra_df = pd.DataFrame(columns=['GEO_Accession'])
            logger.warning("No SRA runs were successfully fetched")

        # Calculate dataset sizes using runinfo data (much faster than vdb-dump)
        dataset_sizes = {}
        if not sra_df.empty:
            logger.info("🔍 Calculating dataset sizes from runinfo data...")
            dataset_sizes = calculate_dataset_sizes_from_runinfo(sra_df)

    # Add dataset sizes to geo_full before merging
    if dataset_sizes:
        geo_full['DatasetSizeBytes'] = geo_full['Accession'].map(dataset_sizes)
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
        left_on='Accession',
        right_on='GEO_Accession',
        how='left'
    )
    logger.info(f"Merged GEO and SRA data for {len(merged)} total records")

    # Select and flag validity
    cols_to_keep = list(geo_full.columns) + ['LibrarySource', 'LibraryLayout', 'LibraryStrategy']
    final = merged[cols_to_keep].copy()

    # Count non-null values for SRA fields
    sra_fields_present = final[['LibrarySource', 'LibraryLayout', 'LibraryStrategy']].notna().sum()
    logger.info(f"Records with SRA fields: LibrarySource={sra_fields_present['LibrarySource']}, "
                f"LibraryLayout={sra_fields_present['LibraryLayout']}, "
                f"LibraryStrategy={sra_fields_present['LibraryStrategy']}")

    final['Valid'] = final.apply(
        lambda row: 'Yes' if (row.get('LibraryLayout') == 'PAIRED' \
                              and row.get('LibrarySource') == 'TRANSCRIPTOMIC') else 'No',
        axis=1
    )

    # Log detailed distribution of LibraryLayout and LibrarySource values
    summary_message = "\n📊 Distribution of Library properties:"
    print(summary_message)
    logger.info(summary_message)

    # Calculate the number of samples with SRA metadata
    samples_with_sra = final['LibraryLayout'].notna().sum()

    # Count LibraryLayout values
    layout_counts = final['LibraryLayout'].value_counts().to_dict()
    summary_message = "\nLibraryLayout distribution:"
    logger.info(summary_message)
    print(summary_message)
    for layout, count in layout_counts.items():
        percentage = count/samples_with_sra*100 if samples_with_sra > 0 else 0
        message = f"  - {layout if pd.notna(layout) else 'NA'}: {count} samples ({percentage:.1f}% of SRA samples)"
        logger.info(message)
        print(message)

    # Count LibrarySource values
    source_counts = final['LibrarySource'].value_counts().to_dict()
    summary_message = "\nLibrarySource distribution:"
    logger.info(summary_message)
    print(summary_message)
    for source, count in source_counts.items():
        percentage = count/samples_with_sra*100 if samples_with_sra > 0 else 0
        message = f"  - {source if pd.notna(source) else 'NA'}: {count} samples ({percentage:.1f}% of SRA samples)"
        logger.info(message)
        print(message)

    # Count LibraryStrategy values
    strategy_counts = final['LibraryStrategy'].value_counts().to_dict()
    summary_message = "\nLibraryStrategy distribution:"
    logger.info(summary_message)
    print(summary_message)
    for strategy, count in strategy_counts.items():
        percentage = count/samples_with_sra*100 if samples_with_sra > 0 else 0
        message = f"  - {strategy if pd.notna(strategy) else 'NA'}: {count} samples ({percentage:.1f}% of SRA samples)"
        logger.info(message)
        print(message)

    # Print validity summary
    valid_counts = final['Valid'].value_counts().to_dict()
    summary_message = "\nValidity summary:"
    logger.info(summary_message)
    print(summary_message)
    yes_count = valid_counts.get('Yes', 0)
    message = f"  - Valid samples (PAIRED + TRANSCRIPTOMIC): {yes_count} ({yes_count/len(final)*100:.1f}% of total)"
    logger.info(message)
    print(message)
    no_count = valid_counts.get('No', 0)
    message = f"  - Invalid samples: {no_count} ({no_count/len(final)*100:.1f}% of total)"
    logger.info(message)
    print(message)
    print()

    # For convenience, generate URLs to PubMed and GEO
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

    final = final.drop_duplicates(subset=['Accession'])


    # Save to output directory
    main_output_file = output_path / "Dataset_identification_result.csv"

    final.to_csv(main_output_file, index=False)
    message = f'Saved combined table to {main_output_file}'
    logger.info(message)
    print(message)

    # Save a separate file with the selected datasets and their clusters if clustering was used
    if args.use_clustering and 'ClusterLabel' in geo_full.columns:
        selected_path = output_path / 'selected_datasets.csv'

        selected_cols = ['ID', 'Accession', 'Title', 'ClusterLabel', 'RelevanceScore']
        selected_df = geo_full[selected_cols].sort_values(['ClusterLabel', 'RelevanceScore'], ascending=[True, False])
        selected_df.to_csv(selected_path, index=False)
        message = f'Saved selected datasets with cluster information to {selected_path}'
        logger.info(message)
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
        logger.info(message)
        print(message)

if __name__ == '__main__':
    main()

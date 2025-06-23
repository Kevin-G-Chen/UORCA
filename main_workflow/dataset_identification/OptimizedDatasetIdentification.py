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
import json
import logging
import os
import sys
import time
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

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

# Load environment variables
load_dotenv()

# Configure Entrez
Entrez.email = os.getenv("ENTREZ_EMAIL", "user@example.com")
if os.getenv("ENTREZ_API_KEY"):
    Entrez.api_key = os.getenv("ENTREZ_API_KEY")

# OpenAI client
client = OpenAI(api_key=os.getenv("OPENAI_API_KEY"))

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
    """Fetch RunInfo data for a GSE accession."""
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

def validate_and_size_datasets(datasets_df: pd.DataFrame, api_delay: float = 0.4) -> pd.DataFrame:
    """Validate datasets and calculate sizes with RNA-seq filtering."""
    def process_dataset(gse_id):
        try:
            logging.info(f"Processing {gse_id}...")
            runinfo_df = fetch_runinfo(gse_id, api_delay)

            if runinfo_df.empty:
                return {
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
                    'original_samples': original_samples,
                    'filtered_samples': filtered_samples,
                    'species_count': unique_species,
                    'species_list': species_list,
                    'valid_dataset': False,
                    'skip_reason': f'Insufficient samples ({filtered_samples} < 3)'
                }

            return {
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
                'original_samples': 0,
                'filtered_samples': 0,
                'species_count': 0,
                'species_list': [],
                'valid_dataset': False,
                'skip_reason': f'Error: {str(e)}'
            }

    # Process all datasets
    logging.info("Validating datasets with RNA-seq filtering...")
    validation_data = []

    for _, row in datasets_df.iterrows():
        gse_id = row['ID']
        validation_info = process_dataset(gse_id)
        validation_data.append({
            'ID': gse_id,
            **validation_info
        })

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

async def assess_relevance(datasets_df: pd.DataFrame, research_query: str) -> pd.DataFrame:
    """Assess relevance of datasets to research query."""
    logging.info(f"Assessing relevance of {len(datasets_df)} datasets...")

    # Prepare data for assessment
    assessment_data = []
    for _, row in datasets_df.iterrows():
        assessment_data.append({
            "ID": row['ID'],
            "Species": row.get('Species', 'Unknown'),
            "Tissue": "Unknown",  # Would need more parsing to extract
            "Technique": row.get('Technique', 'RNA-seq'),
            "Summary": row.get('Summary', '')
        })

    # Call OpenAI for assessment
    prompt = load_prompt("./main_workflow/prompts/dataset_identification/assess_relevance.txt")
    user_input = f"Research Query: {research_query}\n\nDatasets:\n{json.dumps(assessment_data, indent=2)}"

    try:
        assessments = call_openai_json(prompt, user_input, Assessments)

        # Create assessment lookup
        assessment_lookup = {a.ID: (a.RelevanceScore, a.Justification) for a in assessments.assessments}

        # Add assessments to dataframe
        datasets_df = datasets_df.copy()
        datasets_df['RelevanceScore'] = datasets_df['ID'].map(lambda x: assessment_lookup.get(x, (0, "No assessment"))[0])
        datasets_df['Justification'] = datasets_df['ID'].map(lambda x: assessment_lookup.get(x, (0, "No assessment"))[1])

        return datasets_df

    except Exception as e:
        logging.warning(f"Error in relevance assessment: {e}")
        # Fallback: assign default scores
        datasets_df = datasets_df.copy()
        datasets_df['RelevanceScore'] = 5
        datasets_df['Justification'] = "Assessment failed"
        return datasets_df

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
                               help='Delay between NCBI API calls (seconds) to respect rate limits')
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
        validated_datasets_df = validate_and_size_datasets(basic_datasets_df, args.api_delay)

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

        # Step 8: Assess relevance
        logging.info("Assessing relevance...")
        import asyncio
        assessed_df = asyncio.run(assess_relevance(representatives_df, research_query))

        # Step 9: Combine results and sort
        # Add invalid datasets with default scores
        invalid_for_output = invalid_datasets_df.copy()
        invalid_for_output['RelevanceScore'] = 0
        invalid_for_output['Justification'] = "Dataset excluded due to filtering criteria"
        invalid_for_output['Cluster'] = -1

        # Combine valid and invalid datasets for final output
        final_results = pd.concat([
            assessed_df,
            invalid_for_output
        ], ignore_index=True)

        # Remove embedding column to prevent CSV malformation
        if 'embedding' in final_results.columns:
            final_results = final_results.drop('embedding', axis=1)

        # Create simplified final dataframe with only essential columns
        final = pd.DataFrame()
        final['Title'] = final_results.get('Title', '')
        final['Accession'] = final_results['ID']
        final['Samples'] = final_results['filtered_samples']
        final['Species'] = final_results.get('Species', 'Unknown')
        final['Summary'] = final_results.get('Summary', '')
        final['valid_dataset'] = final_results['valid_dataset']
        final['skip_reason'] = final_results['skip_reason']
        final['RelevanceScore'] = final_results['RelevanceScore']
        final['Justification'] = final_results['Justification']
        final['PrimaryPubMed'] = None  # We don't have PubMed IDs in this workflow
        final['DatasetSizeGB'] = 0.0  # Default to 0 since we don't calculate this

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
            multi_df = final[['Accession', 'Species', 'PrimaryPubMed', 'RelevanceScore', 'valid_dataset', 'DatasetSizeGB']].copy()
            multi_df = multi_df[multi_df['valid_dataset'] == True]
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
        print(f"Valid datasets: {len(final[final['valid_dataset'] == True])}")
        print(f"Invalid datasets: {len(final[final['valid_dataset'] == False])}")

        print(f"\nTop {min(len(final), args.max_evaluate or 10)} Results:")
        print("-" * 80)

        for i, (_, row) in enumerate(final.head(args.max_evaluate or 10).iterrows()):
            print(f"\n{i+1}. Dataset: {row['Accession']}")
            print(f"   Title: {row['Title']}")
            print(f"   RNA-seq samples: {row['Samples']}")
            print(f"   Species: {row['Species']}")
            print(f"   Valid for analysis: {'✅ Yes' if row['valid_dataset'] else '❌ No'}")
            if not row['valid_dataset']:
                print(f"   Skip reason: {row['skip_reason']}")
            print(f"   Relevance: {row['RelevanceScore']}/10")
            print(f"   Justification: {row['Justification']}")
            print(f"   GEO URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={row['Accession']}")

    except Exception as e:
        logging.error(f"Error in main execution: {e}")
        raise

if __name__ == "__main__":
    main()

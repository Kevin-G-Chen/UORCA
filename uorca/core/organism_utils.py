"""Organism grouping and filtering utilities."""
from typing import Dict, List, Set
import pandas as pd


def group_datasets_by_organism(
    analysis_info: Dict[str, Dict],
    selected_datasets: List[str]
) -> Dict[str, List[str]]:
    """
    Group datasets by their organism.

    Args:
        analysis_info: Analysis metadata dict (from ResultsIntegrator)
        selected_datasets: List of dataset IDs

    Returns:
        Dict mapping organism names to lists of dataset IDs
    """
    organism_groups = {}
    for dataset_id in selected_datasets:
        if dataset_id in analysis_info:
            organism = analysis_info[dataset_id].get('organism', 'Unknown')
            if organism not in organism_groups:
                organism_groups[organism] = []
            organism_groups[organism].append(dataset_id)
    return organism_groups


def group_contrasts_by_organism(
    analysis_info: Dict[str, Dict],
    selected_contrasts: List[tuple]
) -> Dict[str, List[tuple]]:
    """
    Group contrasts by their dataset's organism.

    Args:
        analysis_info: Analysis metadata dict
        selected_contrasts: List of (analysis_id, contrast_id) tuples

    Returns:
        Dict mapping organism names to lists of contrast tuples
    """
    organism_groups = {}
    for analysis_id, contrast_id in selected_contrasts:
        if analysis_id in analysis_info:
            organism = analysis_info[analysis_id].get('organism', 'Unknown')
            if organism not in organism_groups:
                organism_groups[organism] = []
            organism_groups[organism].append((analysis_id, contrast_id))
    return organism_groups


def filter_genes_by_organism(
    deg_data: Dict[str, Dict[str, pd.DataFrame]],
    cpm_data: Dict[str, pd.DataFrame],
    analysis_info: Dict[str, Dict],
    genes: List[str],
    organism: str,
    selected_contrasts: List[tuple]
) -> List[str]:
    """
    Filter genes to only include those found in datasets of a specific organism.

    Args:
        deg_data: DEG data dict (analysis_id -> contrast_id -> DataFrame)
        cpm_data: CPM data dict (analysis_id -> DataFrame)
        analysis_info: Analysis metadata dict
        genes: List of gene names to filter
        organism: Target organism name
        selected_contrasts: List of (analysis_id, contrast_id) tuples

    Returns:
        List of genes present in the organism's datasets
    """
    organism_genes: Set[str] = set()

    # Collect all genes from datasets of this organism
    for analysis_id, contrast_id in selected_contrasts:
        if analysis_id in analysis_info and analysis_info[analysis_id].get('organism') == organism:
            # Check DEG data
            if analysis_id in deg_data and contrast_id in deg_data[analysis_id]:
                deg_df = deg_data[analysis_id][contrast_id]
                if 'Gene' in deg_df.columns:
                    organism_genes.update(deg_df['Gene'].tolist())

            # Check CPM data
            if analysis_id in cpm_data:
                cpm_df = cpm_data[analysis_id]
                if 'Gene' in cpm_df.columns:
                    organism_genes.update(cpm_df['Gene'].tolist())

    # Return only genes in input list AND found in organism's data
    return [gene for gene in genes if gene in organism_genes]


def filter_genes_by_organism_datasets(
    cpm_data: Dict[str, pd.DataFrame],
    analysis_info: Dict[str, Dict],
    genes: List[str],
    organism: str,
    selected_datasets: List[str]
) -> List[str]:
    """
    Filter genes to only include those found in datasets of a specific organism (for datasets).

    Args:
        cpm_data: CPM data dict (analysis_id -> DataFrame)
        analysis_info: Analysis metadata dict
        genes: List of gene names to filter
        organism: Target organism name
        selected_datasets: List of dataset IDs for this organism

    Returns:
        List of genes that are present in the organism's datasets
    """
    organism_genes: Set[str] = set()

    # Collect all genes from datasets of this organism
    for analysis_id in selected_datasets:
        if analysis_id in analysis_info and analysis_info[analysis_id].get('organism') == organism:
            # Check CPM data (primary source for expression plots)
            if analysis_id in cpm_data:
                cpm_df = cpm_data[analysis_id]
                if 'Gene' in cpm_df.columns:
                    organism_genes.update(cpm_df['Gene'].tolist())

    # Return only genes in input list AND found in organism's data
    return [gene for gene in genes if gene in organism_genes]


def get_organism_display_name(organism: str) -> str:
    """
    Get a user-friendly display name for an organism.

    Args:
        organism: Scientific or common name of organism

    Returns:
        Cleaned display name for UI
    """
    if not organism or organism == 'Unknown':
        return 'Unknown Species'

    organism = organism.strip()

    # Capitalize first letter of each word for scientific names
    if len(organism.split()) <= 2:
        return organism.title()
    else:
        return organism


def get_organisms_from_datasets(
    analysis_info: Dict[str, Dict],
    selected_datasets: List[str]
) -> List[str]:
    """
    Get unique organisms from selected datasets.

    Args:
        analysis_info: Analysis metadata dict
        selected_datasets: List of dataset IDs

    Returns:
        Sorted list of unique organism names
    """
    organisms = set()
    for dataset_id in selected_datasets:
        if dataset_id in analysis_info:
            organism = analysis_info[dataset_id].get('organism', 'Unknown')
            organisms.add(organism)
    return sorted(list(organisms))
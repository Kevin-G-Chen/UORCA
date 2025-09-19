"""
Enhanced Ortholog Mapper Module for UORCA Explorer using Ensembl Compara.

This module provides functionality to map genes between species using a combination
of MyGene API for symbol resolution and Ensembl Compara REST API for accurate
ortholog detection. It provides orthology types, confidence scores, and percent identity.
"""

import logging
import pandas as pd
import mygene
import requests
import time
from typing import List, Dict, Set, Tuple, Optional, Any
from functools import lru_cache
import streamlit as st

logger = logging.getLogger(__name__)

# Initialize MyGene client
MG = mygene.MyGeneInfo()

# Ensembl REST API base URL
ENSEMBL_REST = "https://rest.ensembl.org"

# Species mappings for taxonomy IDs and Ensembl names
SPECIES_MAP = {
    'Homo sapiens': {'taxid': 9606, 'ensembl': 'homo_sapiens'},
    'human': {'taxid': 9606, 'ensembl': 'homo_sapiens'},
    'Mus musculus': {'taxid': 10090, 'ensembl': 'mus_musculus'},
    'mouse': {'taxid': 10090, 'ensembl': 'mus_musculus'},
    'Rattus norvegicus': {'taxid': 10116, 'ensembl': 'rattus_norvegicus'},
    'rat': {'taxid': 10116, 'ensembl': 'rattus_norvegicus'},
    'Danio rerio': {'taxid': 7955, 'ensembl': 'danio_rerio'},
    'zebrafish': {'taxid': 7955, 'ensembl': 'danio_rerio'},
    'Drosophila melanogaster': {'taxid': 7227, 'ensembl': 'drosophila_melanogaster'},
    'fruit fly': {'taxid': 7227, 'ensembl': 'drosophila_melanogaster'},
    'Caenorhabditis elegans': {'taxid': 6239, 'ensembl': 'caenorhabditis_elegans'},
    'worm': {'taxid': 6239, 'ensembl': 'caenorhabditis_elegans'},
}


def get_species_info(organism: str) -> Dict[str, Any]:
    """
    Get taxonomy ID and Ensembl name for an organism.

    Args:
        organism: Organism name (case-insensitive)

    Returns:
        Dict with 'taxid' and 'ensembl' keys, or None if not found
    """
    organism_lower = organism.lower()

    # Try exact match first
    for key, info in SPECIES_MAP.items():
        if key.lower() == organism_lower:
            return info

    # Try partial match
    for key, info in SPECIES_MAP.items():
        if key.lower() in organism_lower or organism_lower in key.lower():
            return info

    logger.warning(f"No species info found for organism: {organism}")
    return None


def get_taxid_from_organism(organism: str) -> Optional[int]:
    """Legacy function for backward compatibility."""
    info = get_species_info(organism)
    return info['taxid'] if info else None


@st.cache_data(show_spinner=False, ttl=3600)
def get_orthologs_by_symbol(
    symbols: Tuple[str, ...],
    source_species: str,
    target_species: str,
    sleep_between: float = 0.08
) -> pd.DataFrame:
    """
    Get orthologs for gene symbols using Ensembl REST API.

    Args:
        symbols: Tuple of gene symbols
        source_species: Source species in Ensembl format (e.g., 'mus_musculus')
        target_species: Target species in Ensembl format (e.g., 'homo_sapiens')
        sleep_between: Sleep time between requests for rate limiting

    Returns:
        DataFrame with ortholog relationships
    """
    rows = []

    for symbol in symbols:
        # Use the documented endpoint format
        url = f"{ENSEMBL_REST}/homology/symbol/{source_species}/{symbol}"
        params = {
            "type": "orthologues",
            "target_species": target_species,
            "format": "condensed",
            "content-type": "application/json"
        }

        try:
            r = requests.get(
                url,
                headers={"Content-Type": "application/json"},
                params=params,
                timeout=20
            )

            if r.status_code == 429:  # Too many requests
                time.sleep(1.0)
                r = requests.get(
                    url,
                    headers={"Content-Type": "application/json"},
                    params=params,
                    timeout=20
                )

            if r.status_code == 200:
                data = r.json().get("data", [])

                for entry in data:
                    src_id = entry.get("id")
                    for h in entry.get("homologies", []):
                        if h.get("species") == target_species:
                            rows.append({
                                "source_symbol": symbol,
                                "source_ensembl": src_id,
                                "target_ensembl": h.get("id"),
                                "target_symbol": h.get("display_name"),  # May be available
                                "orthology_type": h.get("type"),
                                "perc_id": h.get("perc_id"),
                                "confidence": h.get("confidence", 0)
                            })

            time.sleep(sleep_between)  # Rate limiting

        except Exception as e:
            logger.warning(f"Error getting orthologs for {symbol}: {e}")
            continue

    return pd.DataFrame(rows)


@st.cache_data(show_spinner=False, ttl=3600)
def ensembl_to_symbols(
    ensembl_ids: Tuple[str, ...],
    species_taxid: int
) -> pd.DataFrame:
    """
    Map Ensembl Gene IDs back to gene symbols using MyGene.

    Args:
        ensembl_ids: Tuple of Ensembl Gene IDs
        species_taxid: NCBI taxonomy ID

    Returns:
        DataFrame with symbol mapping
    """
    if not ensembl_ids:
        return pd.DataFrame()

    try:
        res = MG.querymany(
            list(ensembl_ids),
            scopes="ensembl.gene",
            species=species_taxid,
            fields="symbol,entrezgene,name,type_of_gene",
            as_dataframe=True
        )

        if res.empty:
            return pd.DataFrame()

        res = res.reset_index().rename(columns={"query": "ensembl_gene_id"})

        # Filter for protein-coding genes
        res = res[res["type_of_gene"] == "protein-coding"].copy()

        return res[["ensembl_gene_id", "symbol", "entrezgene", "name"]]

    except Exception as e:
        logger.error(f"Error mapping Ensembl IDs to symbols: {e}")
        return pd.DataFrame()


def expand_genes_with_orthologs_ensembl(
    input_genes: List[str],
    source_organism: str,
    target_organism: str,
    orthology_filter: str = "all"
) -> Tuple[List[str], pd.DataFrame]:
    """
    Expand gene list with orthologs using Ensembl Compara.

    Args:
        input_genes: List of input gene symbols
        source_organism: Source organism name
        target_organism: Target organism name
        orthology_filter: Filter for orthology type ('all', 'one2one', 'one2many', 'many2many')

    Returns:
        Tuple of (expanded_genes_list, detailed_mapping_dataframe)
    """
    # Get species information
    source_info = get_species_info(source_organism)
    target_info = get_species_info(target_organism)

    if not source_info or not target_info:
        logger.error(f"Species not supported: {source_organism} or {target_organism}")
        return input_genes, pd.DataFrame()

    # Get orthologs directly by symbol
    orthologs = get_orthologs_by_symbol(
        tuple(input_genes),
        source_info['ensembl'],
        target_info['ensembl']
    )

    if orthologs.empty:
        logger.info("No orthologs found")
        # Create empty result with input genes
        result_df = pd.DataFrame({
            "query_symbol": input_genes,
            "source_symbol": input_genes,
            "target_symbol": pd.NA,
            "orthology_type": pd.NA,
            "perc_id": pd.NA,
            "confidence": pd.NA
        })
        return input_genes, result_df

    # Apply orthology type filter if specified
    if orthology_filter != "all":
        if orthology_filter == "one2one":
            orthologs = orthologs[orthologs["orthology_type"] == "ortholog_one2one"]
        elif orthology_filter == "one2many":
            orthologs = orthologs[orthologs["orthology_type"] == "ortholog_one2many"]
        elif orthology_filter == "many2many":
            orthologs = orthologs[orthologs["orthology_type"] == "ortholog_many2many"]

    # If target_symbol is not available, resolve via MyGene
    if "target_symbol" not in orthologs.columns or orthologs["target_symbol"].isna().all():
        target_ensembl_ids = orthologs["target_ensembl"].dropna().unique().tolist()
        if target_ensembl_ids:
            target_symbols = ensembl_to_symbols(tuple(target_ensembl_ids), target_info['taxid'])
            if not target_symbols.empty:
                orthologs = orthologs.merge(
                    target_symbols[["ensembl_gene_id", "symbol"]],
                    left_on="target_ensembl",
                    right_on="ensembl_gene_id",
                    how="left"
                )
                orthologs["target_symbol"] = orthologs["symbol"]
                orthologs = orthologs.drop(columns=["symbol", "ensembl_gene_id"])

    # Create expanded gene list
    expanded_genes = set(input_genes)

    if "target_symbol" in orthologs.columns:
        target_syms = orthologs["target_symbol"].dropna().unique().tolist()
        expanded_genes.update(target_syms)

    # Clean up the output dataframe
    orthologs["query_symbol"] = orthologs["source_symbol"]

    output_cols = [
        "query_symbol", "source_symbol", "source_ensembl",
        "target_symbol", "target_ensembl",
        "orthology_type", "perc_id", "confidence"
    ]

    for col in output_cols:
        if col not in orthologs.columns:
            orthologs[col] = pd.NA

    final_df = orthologs[output_cols].sort_values(
        ["query_symbol", "target_symbol"],
        na_position="last"
    ).reset_index(drop=True)

    # Add missing input genes to the result
    found_genes = set(final_df["query_symbol"].dropna().unique())
    missing_genes = set(input_genes) - found_genes

    if missing_genes:
        missing_df = pd.DataFrame({
            "query_symbol": list(missing_genes),
            "source_symbol": list(missing_genes),
            "source_ensembl": pd.NA,
            "target_symbol": pd.NA,
            "target_ensembl": pd.NA,
            "orthology_type": pd.NA,
            "perc_id": pd.NA,
            "confidence": pd.NA
        })
        final_df = pd.concat([final_df, missing_df], ignore_index=True)

    return list(expanded_genes), final_df


def expand_genes_all_vs_all(
    input_genes: List[str],
    organisms: List[str],
    return_mapping: bool = False
) -> Tuple[List[str], Optional[Dict[str, Dict[str, List[str]]]]]:
    """
    Expand genes using all-vs-all approach with Ensembl Compara.

    Args:
        input_genes: List of input gene symbols
        organisms: List of all organism names
        return_mapping: If True, return detailed mapping

    Returns:
        Tuple of (expanded_gene_list, mapping_dict)
    """
    if not input_genes or len(organisms) < 2:
        return input_genes, None

    expanded_genes = set(input_genes)
    combined_mapping = {} if return_mapping else None

    # Try each pair of organisms
    for source_org in organisms:
        for target_org in organisms:
            if source_org == target_org:
                continue

            # Get orthologs for this pair
            expanded, df = expand_genes_with_orthologs_ensembl(
                input_genes,
                source_org,
                target_org,
                orthology_filter="all"
            )

            # Add new genes to the set
            expanded_genes.update(expanded)

            # Build mapping if requested
            if return_mapping and not df.empty:
                for _, row in df.iterrows():
                    query = row["query_symbol"]
                    target_sym = row["target_symbol"]

                    if pd.notna(query) and pd.notna(target_sym):
                        # Don't include exact self-matches (same case and name)
                        # But DO include case-different orthologs like Cd74 -> CD74
                        if query != target_sym:
                            if query not in combined_mapping:
                                combined_mapping[query] = {}
                            if target_org not in combined_mapping[query]:
                                combined_mapping[query][target_org] = []
                            if target_sym not in combined_mapping[query][target_org]:
                                combined_mapping[query][target_org].append(target_sym)

    logger.info(f"All-vs-all expansion: {len(input_genes)} -> {len(expanded_genes)} genes")

    return list(expanded_genes), combined_mapping


def create_ortholog_mapping_table(
    input_genes: List[str],
    mapping: Dict[str, Dict[str, List[str]]],
    organisms: List[str]
) -> pd.DataFrame:
    """
    Create a DataFrame table showing ortholog mappings.

    Args:
        input_genes: Original input gene list
        mapping: Ortholog mapping dictionary
        organisms: List of all organisms (not used in new format)

    Returns:
        DataFrame with input genes and their orthologs
    """
    data = []

    for gene in input_genes:
        row = {'Input Gene': gene}

        # Collect all unique orthologs (excluding exact self-matches)
        all_orthologs = []
        if gene in mapping:
            for org in mapping[gene]:
                for ortholog in mapping[gene][org]:
                    # Don't filter case-different orthologs (Cd74 -> CD74 is valid)
                    if ortholog != gene and ortholog not in all_orthologs:
                        all_orthologs.append(ortholog)

        # Add orthologs as numbered columns
        for i, ortholog in enumerate(all_orthologs, 1):
            row[f'Ortholog {i}'] = ortholog

        # If no orthologs found, add empty column
        if not all_orthologs:
            row['Ortholog 1'] = ''

        data.append(row)

    df = pd.DataFrame(data)
    return df.fillna('')


def create_hierarchical_gene_list(
    input_genes: List[str],
    mapping: Dict[str, Dict[str, List[str]]],
    expanded_genes: List[str]
) -> str:
    """
    Create a simple text list of all genes (one per line, no indentation).

    Args:
        input_genes: Original input gene list
        mapping: Ortholog mapping dictionary
        expanded_genes: All genes including orthologs

    Returns:
        Formatted string with one gene per line
    """
    # Collect all unique genes
    all_genes = set(input_genes)

    # Add all orthologs (excluding exact self-matches)
    for gene in input_genes:
        if gene in mapping:
            for org, orthologs in mapping[gene].items():
                for ortholog in orthologs:
                    # Don't filter case-different orthologs
                    if ortholog != gene:
                        all_genes.add(ortholog)

    # Sort alphabetically
    sorted_genes = sorted(all_genes)

    return '\n'.join(sorted_genes)


def get_ortholog_summary(
    mapping: Dict[str, Dict[str, List[str]]]
) -> str:
    """
    Generate a summary of ortholog mappings.

    Args:
        mapping: Ortholog mapping dictionary

    Returns:
        Summary string
    """
    if not mapping:
        return "No ortholog mappings found."

    total_orthologs = 0
    genes_with_orthologs = 0

    for gene, orgs in mapping.items():
        if orgs:
            has_orthologs = False
            for org, orthologs in orgs.items():
                if orthologs:
                    has_orthologs = True
                    total_orthologs += len(orthologs)
            if has_orthologs:
                genes_with_orthologs += 1

    lines = [
        f"Found orthologs for {genes_with_orthologs}/{len(mapping)} input genes",
        f"Total orthologs added: {total_orthologs}"
    ]

    return "\n".join(lines)
"""
Ortholog Mapper Module for UORCA Explorer.

This module provides functionality to map genes between species using the MyGene API.
It's designed to expand gene lists by finding orthologs across different species,
specifically for use in the heatmap and expression plots tabs.
"""

import logging
import pandas as pd
import mygene
from typing import List, Dict, Set, Tuple, Optional
from functools import lru_cache
import streamlit as st

logger = logging.getLogger(__name__)

# Species mappings: organism name -> NCBI taxonomy ID
SPECIES_TAXONOMY_MAP = {
    'Homo sapiens': 9606,
    'human': 9606,
    'Mus musculus': 10090,
    'mouse': 10090,
    'Rattus norvegicus': 10116,
    'rat': 10116,
    'Danio rerio': 7955,
    'zebrafish': 7955,
    'Drosophila melanogaster': 7227,
    'fruit fly': 7227,
    'fly': 7227,
    'Caenorhabditis elegans': 6239,
    'worm': 6239,
    'Macaca mulatta': 9544,
    'rhesus macaque': 9544,
    'macaque': 9544,
    'Sus scrofa': 9823,
    'pig': 9823,
    'Gallus gallus': 9031,
    'chicken': 9031,
    'Bos taurus': 9913,
    'cow': 9913,
    'cattle': 9913,
    'Oryctolagus cuniculus': 9986,
    'rabbit': 9986,
    'Canis lupus familiaris': 9615,
    'dog': 9615,
    'Felis catus': 9685,
    'cat': 9685,
    'Equus caballus': 9796,
    'horse': 9796,
    'Ovis aries': 9940,
    'sheep': 9940,
    'Pan troglodytes': 9598,
    'chimpanzee': 9598,
    'chimp': 9598,
    'Xenopus tropicalis': 8364,
    'xenopus': 8364,
    'frog': 8364
}


def get_taxid_from_organism(organism: str) -> Optional[int]:
    """
    Get NCBI taxonomy ID from organism name.

    Args:
        organism: Organism name (case-insensitive)

    Returns:
        Taxonomy ID or None if not found
    """
    # Try exact match first (case-insensitive)
    organism_lower = organism.lower()
    for key, taxid in SPECIES_TAXONOMY_MAP.items():
        if key.lower() == organism_lower:
            return taxid

    # Try partial match
    for key, taxid in SPECIES_TAXONOMY_MAP.items():
        if key.lower() in organism_lower or organism_lower in key.lower():
            return taxid

    logger.warning(f"No taxonomy ID found for organism: {organism}")
    return None


@st.cache_data(show_spinner=False, ttl=3600)
def cached_ortholog_query(
    genes: Tuple[str, ...],
    input_species: int,
    target_species: Tuple[int, ...]
) -> pd.DataFrame:
    """
    Cached function to query MyGene for ortholog mappings.

    Args:
        genes: Tuple of gene symbols to query
        input_species: Input species taxonomy ID
        target_species: Tuple of target species taxonomy IDs

    Returns:
        DataFrame with ortholog mappings
    """
    mg = mygene.MyGeneInfo()

    # Query for input genes
    fields = "symbol,name,entrezgene,ensembl.gene,homologene.genes,homologene.id"

    try:
        result_df = mg.querymany(
            list(genes),
            scopes="symbol",
            species=input_species,
            fields=fields,
            as_dataframe=True,
            returnall=False,
        )

        if result_df.empty:
            return pd.DataFrame()

        # Process results
        result_df = result_df.copy()
        result_df.index.name = "query_symbol"
        result_df.reset_index(inplace=True)

        return result_df
    except Exception as e:
        logger.error(f"Error querying MyGene: {e}")
        return pd.DataFrame()


def extract_orthologs_for_species(
    homologene_genes: List[List],
    target_taxids: Set[int]
) -> Dict[int, List[int]]:
    """
    Extract orthologs for specific species from homologene data.

    Args:
        homologene_genes: List of [taxid, entrez_id] pairs
        target_taxids: Set of target taxonomy IDs

    Returns:
        Dict mapping taxid to list of entrez IDs
    """
    orthologs = {taxid: [] for taxid in target_taxids}

    if not isinstance(homologene_genes, (list, tuple)):
        return orthologs

    for pair in homologene_genes:
        try:
            taxid, eid = pair
            taxid = int(taxid)
            if taxid in target_taxids and pd.notnull(eid):
                orthologs[taxid].append(int(eid))
        except Exception:
            continue

    return orthologs


@st.cache_data(show_spinner=False, ttl=3600)
def cached_entrez_to_symbol(
    entrez_ids: Tuple[int, ...],
    species: int
) -> Dict[int, str]:
    """
    Convert Entrez IDs to gene symbols.

    Args:
        entrez_ids: Tuple of Entrez IDs
        species: Species taxonomy ID

    Returns:
        Dict mapping Entrez ID to gene symbol
    """
    if not entrez_ids:
        return {}

    mg = mygene.MyGeneInfo()

    try:
        result_df = mg.querymany(
            list(entrez_ids),
            scopes="entrezgene",
            species=species,
            fields="symbol,entrezgene",
            as_dataframe=True
        )

        if result_df.empty:
            return {}

        # Create mapping
        mapping = {}
        for _, row in result_df.iterrows():
            if pd.notnull(row.get('entrezgene')) and pd.notnull(row.get('symbol')):
                try:
                    entrez = int(row['entrezgene'])
                    mapping[entrez] = row['symbol']
                except (ValueError, TypeError):
                    continue

        return mapping
    except Exception as e:
        logger.error(f"Error converting Entrez IDs to symbols: {e}")
        return {}


def expand_genes_with_orthologs(
    input_genes: List[str],
    input_organism: str,
    target_organisms: List[str],
    return_mapping: bool = False
) -> Tuple[List[str], Optional[Dict[str, Dict[str, List[str]]]]]:
    """
    Expand a gene list by finding orthologs in target species.

    Args:
        input_genes: List of input gene symbols
        input_organism: Name of input organism
        target_organisms: List of target organism names
        return_mapping: If True, return detailed mapping information

    Returns:
        Tuple of (expanded_gene_list, mapping_dict)
        mapping_dict is None if return_mapping is False
    """
    if not input_genes or not target_organisms:
        return input_genes, None

    # Get taxonomy IDs
    input_taxid = get_taxid_from_organism(input_organism)
    if not input_taxid:
        logger.warning(f"Could not determine taxonomy ID for input organism: {input_organism}")
        return input_genes, None

    target_taxids = {}
    for org in target_organisms:
        taxid = get_taxid_from_organism(org)
        if taxid and taxid != input_taxid:  # Don't map to same species
            target_taxids[org] = taxid

    if not target_taxids:
        logger.info("No valid target species for ortholog mapping")
        return input_genes, None

    # Query MyGene for input genes
    genes_tuple = tuple(input_genes)
    target_taxids_tuple = tuple(target_taxids.values())

    df = cached_ortholog_query(genes_tuple, input_taxid, target_taxids_tuple)

    if df.empty:
        return input_genes, None

    # Process results
    expanded_genes = set(input_genes)  # Start with original genes
    mapping = {} if return_mapping else None

    # Extract orthologs for each gene
    for _, row in df.iterrows():
        query_gene = row.get('query_symbol')
        homologene_genes = row.get('homologene.genes')

        if pd.isnull(homologene_genes).all() if hasattr(pd.isnull(homologene_genes), 'all') else pd.isnull(homologene_genes):
            continue

        # Extract orthologs for target species
        orthologs = extract_orthologs_for_species(
            homologene_genes,
            set(target_taxids.values())
        )

        # Convert Entrez IDs to symbols and add to expanded list
        if return_mapping and query_gene:
            mapping[query_gene] = {}

        for org, taxid in target_taxids.items():
            if taxid in orthologs and orthologs[taxid]:
                # Convert Entrez IDs to symbols
                entrez_tuple = tuple(orthologs[taxid])
                symbol_map = cached_entrez_to_symbol(entrez_tuple, taxid)

                symbols = list(symbol_map.values())
                expanded_genes.update(symbols)

                if return_mapping and query_gene:
                    mapping[query_gene][org] = symbols

    # Log statistics
    n_original = len(input_genes)
    n_expanded = len(expanded_genes)
    n_added = n_expanded - n_original
    logger.info(f"Ortholog expansion: {n_original} input genes -> {n_expanded} total genes (+{n_added} orthologs)")

    return list(expanded_genes), mapping


def expand_genes_all_vs_all(
    input_genes: List[str],
    organisms: List[str],
    return_mapping: bool = False
) -> Tuple[List[str], Optional[Dict[str, Dict[str, List[str]]]]]:
    """
    Expand a gene list by treating each gene as potentially from any species
    and finding orthologs in all other species (all-vs-all approach).

    Args:
        input_genes: List of input gene symbols
        organisms: List of all organism names in the datasets
        return_mapping: If True, return detailed mapping information

    Returns:
        Tuple of (expanded_gene_list, mapping_dict)
        mapping_dict is None if return_mapping is False
    """
    if not input_genes or len(organisms) < 2:
        return input_genes, None

    expanded_genes = set(input_genes)  # Start with original genes
    combined_mapping = {} if return_mapping else None

    # Try each organism as the potential source
    for source_org in organisms:
        # Get target organisms (all except source)
        target_orgs = [org for org in organisms if org != source_org]

        # Try to expand genes assuming they're from source_org
        expanded, mapping = expand_genes_with_orthologs(
            input_genes,
            source_org,
            target_orgs,
            return_mapping=True
        )

        # Add any new genes found
        if expanded:
            expanded_genes.update(expanded)

        # Merge mappings if requested
        if return_mapping and mapping:
            for gene, org_map in mapping.items():
                if gene not in combined_mapping:
                    combined_mapping[gene] = {}
                for target_org, orthologs in org_map.items():
                    if target_org not in combined_mapping[gene]:
                        combined_mapping[gene][target_org] = []
                    # Add unique orthologs
                    for ortholog in orthologs:
                        if ortholog not in combined_mapping[gene][target_org]:
                            combined_mapping[gene][target_org].append(ortholog)

    # Log statistics
    n_original = len(input_genes)
    n_expanded = len(expanded_genes)
    n_added = n_expanded - n_original
    logger.info(f"All-vs-all ortholog expansion: {n_original} input genes -> {n_expanded} total genes (+{n_added} orthologs)")

    return list(expanded_genes), combined_mapping


def get_ortholog_summary(
    mapping: Dict[str, Dict[str, List[str]]]
) -> str:
    """
    Generate a human-readable summary of ortholog mappings.

    Args:
        mapping: Ortholog mapping dictionary

    Returns:
        Summary string
    """
    if not mapping:
        return "No ortholog mappings found."

    lines = []
    total_orthologs = 0
    genes_with_orthologs = 0

    for gene, orgs in mapping.items():
        if orgs:
            genes_with_orthologs += 1
            for org, orthologs in orgs.items():
                total_orthologs += len(orthologs)

    lines.append(f"Found orthologs for {genes_with_orthologs}/{len(mapping)} input genes")
    lines.append(f"Total orthologs added: {total_orthologs}")

    # Show some examples
    examples = []
    for gene, orgs in list(mapping.items())[:3]:
        for org, orthologs in orgs.items():
            if orthologs:
                examples.append(f"  {gene} -> {', '.join(orthologs[:2])} ({org})")
                if len(examples) >= 5:
                    break
        if len(examples) >= 5:
            break

    if examples:
        lines.append("\nExamples:")
        lines.extend(examples)

    return "\n".join(lines)


def filter_genes_by_available_species(
    genes: List[str],
    available_organisms: List[str],
    primary_organism: Optional[str] = None
) -> Tuple[List[str], List[str]]:
    """
    Determine which genes should be looked up for orthologs.

    Args:
        genes: Input gene list
        available_organisms: List of organisms in the selected datasets
        primary_organism: Optional primary organism to infer from genes

    Returns:
        Tuple of (genes_to_expand, target_organisms)
    """
    # If only one organism, no ortholog expansion needed
    if len(set(available_organisms)) <= 1:
        return [], []

    # If we can infer the primary organism from the data, use it
    # Otherwise, we'll need to try multiple species as input
    if primary_organism:
        target_orgs = [org for org in available_organisms if org != primary_organism]
        return genes, target_orgs
    else:
        # For now, assume human as default input species if multiple organisms present
        # This could be enhanced with gene name pattern detection
        if 'Homo sapiens' in available_organisms or 'human' in available_organisms:
            target_orgs = [org for org in available_organisms if org not in ['Homo sapiens', 'human']]
            return genes, target_orgs
        elif 'Mus musculus' in available_organisms or 'mouse' in available_organisms:
            target_orgs = [org for org in available_organisms if org not in ['Mus musculus', 'mouse']]
            return genes, target_orgs
        else:
            # Use first organism as input
            primary = available_organisms[0]
            target_orgs = [org for org in available_organisms if org != primary]
            return genes, target_orgs
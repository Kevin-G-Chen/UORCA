"""
Service for interacting with NCBI GEO (Gene Expression Omnibus) API.
"""
import os
import time
import logging
import requests
import pandas as pd
import xml.etree.ElementTree as ET
from typing import Dict, List, Optional, Any, Union, Tuple
from datetime import datetime
import json
import io
import gzip
from Bio import Entrez

from ..models.geo_models import (
    GEOSearchQuery, 
    GEOSearchResult, 
    GEODataset, 
    GEOSeries,
    GEOPlatform,
    GEOSample,
    ExpressionMatrix,
    GEODatasetWithExpression,
    OmicsType
)
from ..config import get_settings

logger = logging.getLogger(__name__)

class GEOService:
    """Service for interacting with NCBI GEO API."""
    
    def __init__(self):
        """Initialize the GEO service."""
        self.settings = get_settings()
        self.base_url = self.settings.geo_api_base_url
        self.last_request_time = 0
        
        # Set up Entrez email if provided
        if self.settings.geo_api_email:
            Entrez.email = self.settings.geo_api_email
    
    def _respect_rate_limit(self):
        """Respect the NCBI rate limit by adding delay if needed."""
        current_time = time.time()
        time_since_last_request = current_time - self.last_request_time
        
        # Calculate required delay based on rate limit
        required_delay = 1.0 / self.settings.geo_api_rate_limit
        
        if time_since_last_request < required_delay:
            sleep_time = required_delay - time_since_last_request
            logger.debug(f"Rate limiting: sleeping for {sleep_time:.2f} seconds")
            time.sleep(sleep_time)
        
        self.last_request_time = time.time()
    
    def search_datasets(self, query: GEOSearchQuery) -> GEOSearchResult:
        """
        Search for GEO datasets using the provided query.
        
        Args:
            query: The search query parameters
            
        Returns:
            The search results containing dataset and series IDs
        """
        self._respect_rate_limit()
        
        # Prepare the search URL and parameters
        url = f"{self.base_url}/esearch.fcgi"
        params = {
            "db": self.settings.geo_api_db,
            "term": query.term,
            "retmax": str(query.max_results),
            "retmode": "json"
        }
        
        try:
            logger.info(f"Searching GEO with query: {query.term}")
            response = requests.get(url, params=params)
            response.raise_for_status()
            
            data = response.json()
            
            # Extract the results
            result = GEOSearchResult(
                query=query,
                total_count=int(data['esearchresult']['count']),
                dataset_ids=[],
                series_ids=[]
            )
            
            # Parse the IDs and categorize them
            for id_str in data['esearchresult'].get('idlist', []):
                # We need to fetch the summary to determine if it's a dataset or series
                summary = self.get_summary(id_str)
                if summary.get('entrytype', '').lower() == 'gds':
                    result.dataset_ids.append(summary.get('accession', ''))
                elif summary.get('entrytype', '').lower() == 'gse':
                    result.series_ids.append(summary.get('accession', ''))
            
            logger.info(f"Found {len(result.dataset_ids)} datasets and {len(result.series_ids)} series")
            return result
            
        except requests.RequestException as e:
            logger.error(f"Error searching GEO: {str(e)}")
            raise
    
    def get_summary(self, id_str: str) -> Dict[str, Any]:
        """
        Get the summary information for a GEO entity by ID.
        
        Args:
            id_str: The GEO entity ID
            
        Returns:
            Dictionary containing the summary information
        """
        self._respect_rate_limit()
        
        # Prepare the summary URL and parameters
        url = f"{self.base_url}/esummary.fcgi"
        params = {
            "db": self.settings.geo_api_db,
            "id": id_str,
            "retmode": "json"
        }
        
        try:
            logger.debug(f"Fetching summary for ID: {id_str}")
            response = requests.get(url, params=params)
            response.raise_for_status()
            
            data = response.json()
            
            # Extract the summary for the first (and only) result
            summary = data['result'][id_str]
            return summary
            
        except requests.RequestException as e:
            logger.error(f"Error fetching summary for {id_str}: {str(e)}")
            raise
    
    def get_dataset(self, dataset_id: str) -> GEODataset:
        """
        Get detailed information about a GEO dataset.
        
        Args:
            dataset_id: The GEO dataset ID (GDS*)
            
        Returns:
            GEODataset object with detailed information
        """
        self._respect_rate_limit()
        
        if not dataset_id.startswith("GDS"):
            raise ValueError(f"Invalid dataset ID format: {dataset_id}")
        
        # Prepare the efetch URL and parameters
        url = f"{self.base_url}/efetch.fcgi"
        params = {
            "db": self.settings.geo_api_db,
            "id": dataset_id.replace("GDS", ""),
            "retmode": "xml"
        }
        
        try:
            logger.info(f"Fetching dataset: {dataset_id}")
            response = requests.get(url, params=params)
            response.raise_for_status()
            
            # Parse the XML response
            root = ET.fromstring(response.content)
            
            # Extract the dataset information
            title = root.findtext(".//Item[@Name='title']", "")
            summary = root.findtext(".//Item[@Name='summary']", "")
            organism = root.findtext(".//Item[@Name='taxon']", "")
            
            # Extract platform information
            platform_elem = root.find(".//Item[@Name='platform']")
            platform = GEOPlatform(
                id=platform_elem.findtext(".//Item[@Name='accession']", ""),
                title=platform_elem.findtext(".//Item[@Name='title']", ""),
                technology=platform_elem.findtext(".//Item[@Name='technology']", ""),
                organism=organism
            )
            
            # Extract sample IDs
            sample_ids = []
            for sample_elem in root.findall(".//Item[@Name='sample']/Item"):
                if sample_elem.get("Name") == "accession":
                    sample_ids.append(sample_elem.text)
            
            # Extract dates
            pub_date_str = root.findtext(".//Item[@Name='pubmed_date']", "")
            submission_date_str = root.findtext(".//Item[@Name='submission_date']", "")
            update_date_str = root.findtext(".//Item[@Name='update_date']", "")
            
            # Convert date strings to datetime objects if available
            pub_date = datetime.strptime(pub_date_str, "%Y/%m/%d") if pub_date_str else None
            submission_date = datetime.strptime(submission_date_str, "%Y/%m/%d") if submission_date_str else None
            update_date = datetime.strptime(update_date_str, "%Y/%m/%d") if update_date_str else None
            
            # Determine omics type based on title and summary
            omics_type = self._determine_omics_type(title, summary)
            
            # Create the GEODataset object
            dataset = GEODataset(
                id=dataset_id,
                title=title,
                summary=summary,
                organism=organism,
                platform=platform,
                samples=sample_ids,
                publication_date=pub_date,
                submission_date=submission_date,
                last_update_date=update_date,
                omics_type=omics_type
            )
            
            return dataset
            
        except requests.RequestException as e:
            logger.error(f"Error fetching dataset {dataset_id}: {str(e)}")
            raise
    
    def get_series(self, series_id: str) -> GEOSeries:
        """
        Get detailed information about a GEO series.
        
        Args:
            series_id: The GEO series ID (GSE*)
            
        Returns:
            GEOSeries object with detailed information
        """
        self._respect_rate_limit()
        
        if not series_id.startswith("GSE"):
            raise ValueError(f"Invalid series ID format: {series_id}")
        
        # Series information is fetched differently using the GEO web API
        url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
        params = {
            "acc": series_id,
            "targ": "self",
            "form": "xml",
            "view": "full"
        }
        
        try:
            logger.info(f"Fetching series: {series_id}")
            response = requests.get(url, params=params)
            response.raise_for_status()
            
            # Parse the XML response
            root = ET.fromstring(response.content)
            
            # Extract series information
            series_elem = root.find(".//Series")
            
            if series_elem is None:
                raise ValueError(f"Series element not found in response for {series_id}")
            
            title = series_elem.findtext("Title", "")
            summary = series_elem.findtext("Summary", "")
            organism = series_elem.findtext("Organism", "")
            
            # Extract platform information
            platform_elems = root.findall(".//Platform")
            platforms = []
            
            for platform_elem in platform_elems:
                platform = GEOPlatform(
                    id=platform_elem.findtext("Accession", ""),
                    title=platform_elem.findtext("Title", ""),
                    technology=platform_elem.findtext("Technology", ""),
                    organism=platform_elem.findtext("Organism", "")
                )
                platforms.append(platform)
            
            # Extract sample IDs
            sample_ids = []
            for sample_elem in root.findall(".//Sample"):
                sample_ids.append(sample_elem.findtext("Accession", ""))
            
            # Extract dates
            pub_date_str = series_elem.findtext("Pubmed-Date", "")
            submission_date_str = series_elem.findtext("Submission-Date", "")
            update_date_str = series_elem.findtext("Last-Update-Date", "")
            
            # Convert date strings to datetime objects if available
            try:
                pub_date = datetime.strptime(pub_date_str, "%Y-%m-%d") if pub_date_str else None
            except ValueError:
                pub_date = None
                
            try:
                submission_date = datetime.strptime(submission_date_str, "%Y-%m-%d") if submission_date_str else None
            except ValueError:
                submission_date = None
                
            try:
                update_date = datetime.strptime(update_date_str, "%Y-%m-%d") if update_date_str else None
            except ValueError:
                update_date = None
            
            # Determine omics type based on title and summary
            omics_type = self._determine_omics_type(title, summary)
            
            # Create the GEOSeries object
            series = GEOSeries(
                id=series_id,
                title=title,
                summary=summary,
                organism=organism,
                platforms=platforms,
                samples=sample_ids,
                publication_date=pub_date,
                submission_date=submission_date,
                last_update_date=update_date,
                omics_type=omics_type
            )
            
            return series
            
        except requests.RequestException as e:
            logger.error(f"Error fetching series {series_id}: {str(e)}")
            raise
    
    def get_sample(self, sample_id: str) -> GEOSample:
        """
        Get detailed information about a GEO sample.
        
        Args:
            sample_id: The GEO sample ID (GSM*)
            
        Returns:
            GEOSample object with detailed information
        """
        self._respect_rate_limit()
        
        if not sample_id.startswith("GSM"):
            raise ValueError(f"Invalid sample ID format: {sample_id}")
        
        # Sample information is fetched similar to series
        url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
        params = {
            "acc": sample_id,
            "targ": "self",
            "form": "xml",
            "view": "full"
        }
        
        try:
            logger.debug(f"Fetching sample: {sample_id}")
            response = requests.get(url, params=params)
            response.raise_for_status()
            
            # Parse the XML response
            root = ET.fromstring(response.content)
            
            # Extract sample information
            sample_elem = root.find(".//Sample")
            
            if sample_elem is None:
                raise ValueError(f"Sample element not found in response for {sample_id}")
            
            title = sample_elem.findtext("Title", "")
            
            # Extract sample attributes (characteristics and other metadata)
            attributes = {}
            
            # Process characteristics
            for char_elem in sample_elem.findall(".//Characteristics"):
                tag = char_elem.get("tag", "characteristic")
                value = char_elem.text.strip() if char_elem.text else ""
                attributes[tag] = value
            
            # Add other relevant metadata
            for elem_name in ["Source", "Organism", "Molecule", "Label", "Description"]:
                value = sample_elem.findtext(elem_name, "")
                if value:
                    attributes[elem_name.lower()] = value
            
            # Create the GEOSample object
            sample = GEOSample(
                id=sample_id,
                title=title,
                attributes=attributes
            )
            
            return sample
            
        except requests.RequestException as e:
            logger.error(f"Error fetching sample {sample_id}: {str(e)}")
            raise
    
    def get_expression_data(self, dataset_id: str) -> Optional[ExpressionMatrix]:
        """
        Get expression data for a GEO dataset.
        
        Args:
            dataset_id: The GEO dataset ID (GDS*)
            
        Returns:
            ExpressionMatrix object containing the expression data, or None if not available
        """
        self._respect_rate_limit()
        
        if not dataset_id.startswith("GDS"):
            raise ValueError(f"Invalid dataset ID format: {dataset_id}")
        
        # Prepare the URL for fetching expression data (SOFT format)
        url = f"https://ftp.ncbi.nlm.nih.gov/geo/datasets/GDSnnn/{dataset_id}/soft/{dataset_id}_full.soft.gz"
        
        try:
            logger.info(f"Fetching expression data for dataset: {dataset_id}")
            response = requests.get(url, stream=True)
            response.raise_for_status()
            
            # Decompress the gzipped content
            content = gzip.decompress(response.content).decode('utf-8')
            
            # Parse the SOFT file to extract expression data
            gene_ids = []
            sample_ids = []
            values = []
            
            # Flag to indicate when we're in the dataset table section
            in_table = False
            header_parsed = False
            
            for line in content.split('\n'):
                line = line.strip()
                
                # Look for the start of the dataset table
                if line.startswith('!dataset_table_begin'):
                    in_table = True
                    continue
                
                # Look for the end of the dataset table
                if line.startswith('!dataset_table_end'):
                    in_table = False
                    break
                
                # Process table data
                if in_table:
                    if not header_parsed and line.startswith('#'):
                        # Parse header line to get sample IDs
                        header_parts = line[1:].split('\t')  # Remove the # and split
                        # First column is usually ID_REF
                        sample_ids = header_parts[1:]
                        header_parsed = True
                    elif not line.startswith('#') and not line.startswith('!'):
                        # Parse data line
                        parts = line.split('\t')
                        if len(parts) > 1:
                            gene_id = parts[0]
                            gene_values = [float(v) if v and v != 'null' else float('nan') for v in parts[1:]]
                            
                            if len(gene_values) == len(sample_ids):
                                gene_ids.append(gene_id)
                                values.append(gene_values)
            
            if not gene_ids or not sample_ids:
                logger.warning(f"No expression data found for dataset: {dataset_id}")
                return None
            
            # Create and return the expression matrix
            matrix = ExpressionMatrix(
                gene_ids=gene_ids,
                sample_ids=sample_ids,
                values=values
            )
            
            logger.info(f"Extracted expression data: {len(gene_ids)} genes, {len(sample_ids)} samples")
            return matrix
            
        except requests.RequestException as e:
            logger.error(f"Error fetching expression data for {dataset_id}: {str(e)}")
            return None
    
    def get_dataset_with_expression(self, dataset_id: str) -> GEODatasetWithExpression:
        """
        Get a GEO dataset with its expression data.
        
        Args:
            dataset_id: The GEO dataset ID (GDS*)
            
        Returns:
            GEODatasetWithExpression object containing both metadata and expression data
        """
        # Get dataset metadata
        dataset = self.get_dataset(dataset_id)
        
        # Get expression data
        expression_matrix = self.get_expression_data(dataset_id)
        
        # Get sample metadata for each sample
        sample_metadata = {}
        for sample_id in dataset.samples:
            try:
                sample = self.get_sample(sample_id)
                sample_metadata[sample_id] = sample.attributes
            except Exception as e:
                logger.warning(f"Error fetching metadata for sample {sample_id}: {str(e)}")
        
        # Create and return the combined object
        return GEODatasetWithExpression(
            metadata=dataset,
            expression_matrix=expression_matrix,
            sample_metadata=sample_metadata
        )
    
    def get_series_with_expression(self, series_id: str) -> GEODatasetWithExpression:
        """
        Get a GEO series with expression data.
        
        Args:
            series_id: The GEO series ID (GSE*)
            
        Returns:
            GEODatasetWithExpression object containing both metadata and expression data
        """
        # Series expression data is more complex and may require additional processing
        # This is a simplified implementation
        
        # Get series metadata
        series = self.get_series(series_id)
        
        # For GSE, we need to manually construct the expression matrix from supplementary files
        # This is a simplified approach; in practice, you might need more sophisticated parsing
        
        # Create a placeholder for expression data
        expression_matrix = None
        
        # Get sample metadata for each sample
        sample_metadata = {}
        for sample_id in series.samples:
            try:
                sample = self.get_sample(sample_id)
                sample_metadata[sample_id] = sample.attributes
            except Exception as e:
                logger.warning(f"Error fetching metadata for sample {sample_id}: {str(e)}")
        
        # Create and return the combined object
        return GEODatasetWithExpression(
            metadata=series,
            expression_matrix=expression_matrix,
            sample_metadata=sample_metadata
        )
    
    def _determine_omics_type(self, title: str, summary: str) -> OmicsType:
        """
        Determine the omics type based on title and summary.
        
        Args:
            title: The dataset/series title
            summary: The dataset/series summary
            
        Returns:
            OmicsType enum value
        """
        # Combine text for searching
        text = (title + " " + summary).lower()
        
        # Look for keywords associated with each omics type
        if any(keyword in text for keyword in ["transcriptome", "rna-seq", "microarray", "gene expression", "mrna"]):
            return OmicsType.TRANSCRIPTOMICS
        elif any(keyword in text for keyword in ["proteome", "protein", "mass spec", "proteomics"]):
            return OmicsType.PROTEOMICS
        elif any(keyword in text for keyword in ["genome", "dna-seq", "exome", "snp", "variant", "mutation"]):
            return OmicsType.GENOMICS
        elif any(keyword in text for keyword in ["metabolome", "metabolite", "metabolomics"]):
            return OmicsType.METABOLOMICS
        elif any(keyword in text for keyword in ["methylation", "chip-seq", "histone", "chromatin", "epigenetic"]):
            return OmicsType.EPIGENOMICS
        else:
            return OmicsType.UNKNOWN

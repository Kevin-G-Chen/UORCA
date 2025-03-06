"""
Unit tests for the GEO service.
"""
import unittest
from unittest.mock import patch, MagicMock
import json
import os
import sys
from pathlib import Path

# Add the parent directory to the path so we can import the module
sys.path.append(str(Path(__file__).parent.parent))

from UORCA.services.geo_service import GEOService
from UORCA.models.geo_models import GEOSearchQuery, OmicsType

class TestGEOService(unittest.TestCase):
    """Tests for the GEO service."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.service = GEOService()
    
    @patch('UORCA.services.geo_service.requests.get')
    def test_search_datasets(self, mock_get):
        """Test searching for datasets."""
        # Mock the response from the GEO API
        mock_response = MagicMock()
        mock_response.json.return_value = {
            "esearchresult": {
                "count": "2",
                "idlist": ["123", "456"]
            }
        }
        mock_response.raise_for_status = MagicMock()
        mock_get.return_value = mock_response
        
        # Mock the get_summary method
        self.service.get_summary = MagicMock()
        self.service.get_summary.side_effect = [
            {"entrytype": "GDS", "accession": "GDS123"},
            {"entrytype": "GSE", "accession": "GSE456"}
        ]
        
        # Create a search query
        query = GEOSearchQuery(
            term="cancer",
            max_results=10,
            omics_type=OmicsType.TRANSCRIPTOMICS,
            organism="human"
        )
        
        # Call the method
        result = self.service.search_datasets(query)
        
        # Check that requests.get was called with the expected arguments
        mock_get.assert_called_once()
        args, kwargs = mock_get.call_args
        self.assertEqual(kwargs['params']['db'], self.service.settings.geo_api_db)
        self.assertEqual(kwargs['params']['term'], query.term)
        self.assertEqual(kwargs['params']['retmax'], str(query.max_results))
        
        # Check that get_summary was called for each ID
        self.assertEqual(self.service.get_summary.call_count, 2)
        
        # Check the result
        self.assertEqual(result.total_count, 2)
        self.assertEqual(result.dataset_ids, ["GDS123"])
        self.assertEqual(result.series_ids, ["GSE456"])
    
    @patch('UORCA.services.geo_service.requests.get')
    def test_get_summary(self, mock_get):
        """Test getting a summary."""
        # Mock the response
        mock_response = MagicMock()
        mock_response.json.return_value = {
            "result": {
                "123": {
                    "entrytype": "GDS",
                    "accession": "GDS123",
                    "title": "Test Dataset"
                }
            }
        }
        mock_response.raise_for_status = MagicMock()
        mock_get.return_value = mock_response
        
        # Call the method
        result = self.service.get_summary("123")
        
        # Check that requests.get was called with the expected arguments
        mock_get.assert_called_once()
        args, kwargs = mock_get.call_args
        self.assertEqual(kwargs['params']['db'], self.service.settings.geo_api_db)
        self.assertEqual(kwargs['params']['id'], "123")
        
        # Check the result
        self.assertEqual(result["entrytype"], "GDS")
        self.assertEqual(result["accession"], "GDS123")
        self.assertEqual(result["title"], "Test Dataset")
    
    @patch('UORCA.services.geo_service.requests.get')
    def test_get_dataset(self, mock_get):
        """Test getting a dataset."""
        # Mock the response
        with open(Path(__file__).parent / "test_data" / "gds_response.xml", "r") as f:
            xml_content = f.read()
        
        mock_response = MagicMock()
        mock_response.content = xml_content.encode('utf-8')
        mock_response.raise_for_status = MagicMock()
        mock_get.return_value = mock_response
        
        # Call the method
        result = self.service.get_dataset("GDS123")
        
        # Check that requests.get was called with the expected arguments
        mock_get.assert_called_once()
        args, kwargs = mock_get.call_args
        self.assertEqual(kwargs['params']['db'], self.service.settings.geo_api_db)
        self.assertEqual(kwargs['params']['id'], "123")
        
        # Check the result
        self.assertEqual(result.id, "GDS123")
        # Add more assertions based on the content of your test XML file

    def test_determine_omics_type(self):
        """Test determining omics type from title and summary."""
        # Test transcriptomics
        omics_type = self.service._determine_omics_type(
            "Gene expression profiling of cancer cells",
            "RNA-seq analysis of tumor samples"
        )
        self.assertEqual(omics_type, OmicsType.TRANSCRIPTOMICS)
        
        # Test proteomics
        omics_type = self.service._determine_omics_type(
            "Proteome analysis of cancer cells",
            "Mass spectrometry-based proteomics"
        )
        self.assertEqual(omics_type, OmicsType.PROTEOMICS)
        
        # Test genomics
        omics_type = self.service._determine_omics_type(
            "Genome sequencing of cancer cells",
            "Whole-genome sequencing of tumor samples"
        )
        self.assertEqual(omics_type, OmicsType.GENOMICS)
        
        # Test metabolomics
        omics_type = self.service._determine_omics_type(
            "Metabolome analysis of cancer cells",
            "Metabolomics study of tumor samples"
        )
        self.assertEqual(omics_type, OmicsType.METABOLOMICS)
        
        # Test epigenomics
        omics_type = self.service._determine_omics_type(
            "Methylation analysis of cancer cells",
            "ChIP-seq study of histone modifications"
        )
        self.assertEqual(omics_type, OmicsType.EPIGENOMICS)
        
        # Test unknown
        omics_type = self.service._determine_omics_type(
            "Analysis of cancer cells",
            "Study of tumor samples"
        )
        self.assertEqual(omics_type, OmicsType.UNKNOWN)

if __name__ == '__main__':
    # Create test data directory if it doesn't exist
    os.makedirs(Path(__file__).parent / "test_data", exist_ok=True)
    
    # Create a sample XML file for testing
    gds_xml = """
    <DocSum>
        <Id>123</Id>
        <Item Name="title">Test Dataset</Item>
        <Item Name="summary">Test Summary</Item>
        <Item Name="taxon">human</Item>
        <Item Name="platform">
            <Item Name="accession">GPL123</Item>
            <Item Name="title">Test Platform</Item>
            <Item Name="technology">Microarray</Item>
        </Item>
        <Item Name="sample">
            <Item Name="accession">GSM123</Item>
        </Item>
        <Item Name="sample">
            <Item Name="accession">GSM456</Item>
        </Item>
    </DocSum>
    """
    
    with open(Path(__file__).parent / "test_data" / "gds_response.xml", "w") as f:
        f.write(gds_xml)
    
    unittest.main()

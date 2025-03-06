"""
Unit tests for the Analysis service.
"""
import unittest
from unittest.mock import patch, MagicMock
import json
import os
import sys
import tempfile
import pandas as pd
from pathlib import Path

# Add the parent directory to the path so we can import the module
sys.path.append(str(Path(__file__).parent.parent))

from UORCA.services.analysis_service import AnalysisService
from UORCA.services.llm_service import LLMService
from UORCA.models.geo_models import (
    GEODataset, 
    GEOPlatform, 
    ExpressionMatrix, 
    GEODatasetWithExpression,
    OmicsType
)
from UORCA.models.analysis_models import (
    AnalysisRequest,
    AnalysisType,
    AnalysisParameter,
    AnalysisStatus,
    AnalysisResult,
    DifferentialExpressionResult,
    PathwayEnrichmentResult,
    Insight,
    InsightType
)

class TestAnalysisService(unittest.TestCase):
    """Tests for the Analysis service."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.service = AnalysisService()
        
        # Create a temp directory for testing
        self.temp_dir = tempfile.TemporaryDirectory()
        self.service.settings.analysis_cache_dir = self.temp_dir.name
        
        # Create a sample dataset for testing
        platform = GEOPlatform(
            id="GPL123",
            title="Test Platform",
            technology="Microarray"
        )
        
        metadata = GEODataset(
            id="GDS123",
            title="Test Dataset",
            summary="Gene expression data from cancer samples",
            organism="human",
            platform=platform,
            samples=["GSM123", "GSM456", "GSM789"],
            omics_type=OmicsType.TRANSCRIPTOMICS
        )
        
        # Create a simple expression matrix
        expression_matrix = ExpressionMatrix(
            gene_ids=["GENE1", "GENE2", "GENE3", "GENE4", "GENE5"],
            sample_ids=["SAMPLE1", "SAMPLE2", "SAMPLE3"],
            values=[
                [1.0, 2.0, 3.0],
                [4.0, 5.0, 6.0],
                [7.0, 8.0, 9.0],
                [10.0, 11.0, 12.0],
                [13.0, 14.0, 15.0]
            ]
        )
        
        # Sample metadata
        sample_metadata = {
            "GSM123": {"group": "control", "tissue": "normal"},
            "GSM456": {"group": "treatment", "tissue": "tumor"},
            "GSM789": {"group": "treatment", "tissue": "tumor"}
        }
        
        self.dataset = GEODatasetWithExpression(
            metadata=metadata,
            expression_matrix=expression_matrix,
            sample_metadata=sample_metadata
        )
        
        # Create a sample analysis request
        self.request = AnalysisRequest(
            research_question="What genes are differentially expressed in cancer?",
            analysis_type=AnalysisType.DIFFERENTIAL_EXPRESSION,
            dataset_ids=["GDS123"],
            omics_types=[OmicsType.TRANSCRIPTOMICS],
            parameters=[
                AnalysisParameter(
                    name="group_var",
                    value="group",
                    description="Variable to use for grouping samples"
                )
            ],
            max_datasets=1
        )
    
    def tearDown(self):
        """Clean up after tests."""
        self.temp_dir.cleanup()
    
    @patch.object(LLMService, 'design_analysis')
    @patch.object(LLMService, 'generate_analysis_code')
    def test_analyze_dataset(self, mock_generate_code, mock_design_analysis):
        """Test analyzing a dataset."""
        # Mock the LLM service methods
        mock_design_analysis.return_value = {
            "analysis_type": "differential_expression",
            "description": "Compare gene expression between treatment and control groups",
            "parameters": [
                {
                    "name": "group_var",
                    "value": "group",
                    "description": "Variable to use for grouping samples"
                }
            ],
            "code_template": "# Python code template for the analysis..."
        }
        
        mock_generate_code.return_value = """
# Sample analysis code
import pandas as pd
import numpy as np

# Load data
expression_df = pd.read_csv(expression_file, index_col=0)
sample_meta_df = pd.read_csv(sample_meta_file, index_col=0)

# Simple differential expression analysis
control_samples = sample_meta_df[sample_meta_df['group'] == 'control'].index
treatment_samples = sample_meta_df[sample_meta_df['group'] == 'treatment'].index

control_mean = expression_df[control_samples].mean(axis=1)
treatment_mean = expression_df[treatment_samples].mean(axis=1)

log2fc = np.log2(treatment_mean / control_mean)
# p-values would normally be calculated with a statistical test
p_values = np.random.uniform(0, 1, len(expression_df))
adj_p_values = p_values * 0.8  # Simple adjustment for testing

# Save results
results = {
    "de_results": {
        "comparison": "treatment vs control",
        "gene_ids": expression_df.index.tolist(),
        "log2_fold_changes": log2fc.tolist(),
        "p_values": p_values.tolist(),
        "adjusted_p_values": adj_p_values.tolist()
    }
}
"""
        
        # Patch _execute_analysis_code to return mock results
        mock_results = {
            "de_results": {
                "comparison": "treatment vs control",
                "gene_ids": ["GENE1", "GENE2", "GENE3", "GENE4", "GENE5"],
                "log2_fold_changes": [1.5, -0.8, 2.3, -1.2, 0.5],
                "p_values": [0.01, 0.04, 0.001, 0.02, 0.07],
                "adjusted_p_values": [0.02, 0.06, 0.003, 0.04, 0.09]
            }
        }
        
        with patch.object(self.service, '_execute_analysis_code', return_value=mock_results):
            # Call the method
            result = self.service.analyze_dataset(self.dataset, self.request)
            
            # Check the result
            self.assertEqual(result.status, AnalysisStatus.COMPLETED)
            self.assertEqual(result.research_question, self.request.research_question)
            self.assertEqual(result.datasets_analyzed, ["GDS123"])
            self.assertEqual(result.omics_types, [OmicsType.TRANSCRIPTOMICS])
            self.assertTrue("GDS123" in result.differential_expression_results)
            
            # Check the differential expression results
            de_result = result.differential_expression_results["GDS123"]
            self.assertEqual(de_result.dataset_id, "GDS123")
            self.assertEqual(de_result.comparison, "treatment vs control")
            self.assertEqual(len(de_result.gene_ids), 5)
            self.assertEqual(len(de_result.log2_fold_changes), 5)
            self.assertEqual(len(de_result.p_values), 5)
            self.assertEqual(len(de_result.adjusted_p_values), 5)
            
            # Check significant genes
            sig_genes = de_result.get_significant_genes(adj_p_threshold=0.05, log2fc_threshold=1.0)
            self.assertEqual(len(sig_genes), 2)  # GENE3 and GENE4 should be significant
    
    @patch.object(LLMService, 'integrate_insights')
    def test_integrate_results(self, mock_integrate_insights):
        """Test integrating results from multiple analyses."""
        # Create sample analysis results
        result1 = AnalysisResult(
            request_id="analysis_123",
            status=AnalysisStatus.COMPLETED,
            research_question="What genes are differentially expressed in cancer?",
            datasets_analyzed=["GDS123"],
            omics_types=[OmicsType.TRANSCRIPTOMICS],
            differential_expression_results={
                "GDS123": DifferentialExpressionResult(
                    dataset_id="GDS123",
                    comparison="tumor vs normal",
                    gene_ids=["GENE1", "GENE2", "GENE3"],
                    log2_fold_changes=[1.5, -0.8, 2.3],
                    p_values=[0.01, 0.04, 0.001],
                    adjusted_p_values=[0.02, 0.06, 0.003]
                )
            }
        )
        
        result2 = AnalysisResult(
            request_id="analysis_456",
            status=AnalysisStatus.COMPLETED,
            research_question="What genes are differentially expressed in cancer?",
            datasets_analyzed=["GDS456"],
            omics_types=[OmicsType.TRANSCRIPTOMICS],
            differential_expression_results={
                "GDS456": DifferentialExpressionResult(
                    dataset_id="GDS456",
                    comparison="tumor vs normal",
                    gene_ids=["GENE1", "GENE3", "GENE4"],
                    log2_fold_changes=[1.2, 2.0, -1.5],
                    p_values=[0.02, 0.003, 0.01],
                    adjusted_p_values=[0.04, 0.006, 0.02]
                )
            }
        )
        
        # Mock the LLM service method
        mock_integrate_insights.return_value = [
            {
                "insight_type": "key_gene",
                "description": "GENE3 is consistently upregulated in tumor samples across datasets",
                "supporting_datasets": ["GDS123", "GDS456"],
                "confidence_score": 0.9,
                "related_entities": ["GENE3", "tumor"]
            },
            {
                "insight_type": "contradiction",
                "description": "GENE1 shows varying fold changes across datasets",
                "supporting_datasets": ["GDS123", "GDS456"],
                "confidence_score": 0.7,
                "related_entities": ["GENE1"]
            }
        ]
        
        # Call the method
        insights = self.service.integrate_results(
            research_question="What genes are differentially expressed in cancer?",
            analysis_results=[result1, result2]
        )
        
        # Check the result
        self.assertEqual(len(insights), 2)
        self.assertEqual(insights[0].insight_type, InsightType.KEY_GENE)
        self.assertEqual(insights[0].description, 
                        "GENE3 is consistently upregulated in tumor samples across datasets")
        self.assertEqual(insights[0].supporting_datasets, ["GDS123", "GDS456"])
        self.assertAlmostEqual(insights[0].confidence_score, 0.9)
        
        self.assertEqual(insights[1].insight_type, InsightType.CONTRADICTION)
        self.assertEqual(insights[1].description, "GENE1 shows varying fold changes across datasets")
        self.assertEqual(insights[1].supporting_datasets, ["GDS123", "GDS456"])
        self.assertAlmostEqual(insights[1].confidence_score, 0.7)
    
    @patch.object(LLMService, 'generate_summary_report')
    def test_generate_summary(self, mock_generate_summary):
        """Test generating a summary."""
        # Create sample insights
        insights = [
            Insight(
                insight_type=InsightType.KEY_GENE,
                description="GENE3 is consistently upregulated in tumor samples across datasets",
                supporting_datasets=["GDS123", "GDS456"],
                confidence_score=0.9,
                related_entities=["GENE3", "tumor"]
            ),
            Insight(
                insight_type=InsightType.CONTRADICTION,
                description="GENE1 shows varying fold changes across datasets",
                supporting_datasets=["GDS123", "GDS456"],
                confidence_score=0.7,
                related_entities=["GENE1"]
            )
        ]
        
        # Mock the LLM service method
        mock_generate_summary.return_value = """
# Analysis Summary

## Research Question
What genes are differentially expressed in cancer?

## Datasets Analyzed
- GDS123
- GDS456

## Key Findings
1. GENE3 is consistently upregulated in tumor samples across datasets (high confidence)
2. GENE1 shows varying fold changes across datasets (medium confidence)

## Implications
These findings suggest that GENE3 may play a crucial role in cancer development and could be a potential biomarker or therapeutic target.
"""
        
        # Call the method
        summary = self.service.generate_summary(
            research_question="What genes are differentially expressed in cancer?",
            datasets_analyzed=["GDS123", "GDS456"],
            insights=insights
        )
        
        # Check the result
        self.assertTrue("GENE3 is consistently upregulated" in summary)
        self.assertTrue("GENE1 shows varying fold changes" in summary)
        self.assertTrue("potential biomarker or therapeutic target" in summary)
    
    def test_cache_results(self):
        """Test caching analysis results."""
        # Create a sample result
        result = AnalysisResult(
            request_id="test_cache",
            status=AnalysisStatus.COMPLETED,
            research_question="What genes are differentially expressed in cancer?",
            datasets_analyzed=["GDS123"],
            omics_types=[OmicsType.TRANSCRIPTOMICS]
        )
        
        # Cache the result
        self.service._cache_results("test_cache", result)
        
        # Check if the file exists
        cache_file = Path(self.temp_dir.name) / "test_cache.json"
        self.assertTrue(cache_file.exists())
        
        # Check the content
        with open(cache_file, 'r') as f:
            cached_data = json.load(f)
            
        self.assertEqual(cached_data["request_id"], "test_cache")
        self.assertEqual(cached_data["research_question"], 
                         "What genes are differentially expressed in cancer?")

if __name__ == '__main__':
    unittest.main()

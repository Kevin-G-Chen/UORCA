"""
Unit tests for the Workflow service.
"""
import unittest
from unittest.mock import patch, MagicMock
import json
import os
import sys
import tempfile
from pathlib import Path

# Add the parent directory to the path so we can import the module
sys.path.append(str(Path(__file__).parent.parent))

from UORCA.services.workflow_service import WorkflowService
from UORCA.services.geo_service import GEOService
from UORCA.services.llm_service import LLMService
from UORCA.services.analysis_service import AnalysisService
from UORCA.models.workflow_models import (
    UORCAWorkflowState,
    WorkflowState,
    DatasetRelevanceScore
)
from UORCA.models.geo_models import (
    GEOSearchQuery,
    GEOSearchResult,
    GEODatasetWithExpression,
    GEODataset,
    GEOPlatform,
    ExpressionMatrix,
    OmicsType
)
from UORCA.models.analysis_models import (
    AnalysisResult,
    AnalysisStatus,
    AnalysisRequest,
    Insight,
    InsightType
)

class TestWorkflowService(unittest.TestCase):
    """Tests for the Workflow service."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.service = WorkflowService()
        
        # Create a temp directory for testing
        self.temp_dir = tempfile.TemporaryDirectory()
        self.service.settings.analysis_cache_dir = self.temp_dir.name
        
        # Mock the services
        self.service.geo_service = MagicMock(spec=GEOService)
        self.service.llm_service = MagicMock(spec=LLMService)
        self.service.analysis_service = MagicMock(spec=AnalysisService)
    
    def tearDown(self):
        """Clean up after tests."""
        self.temp_dir.cleanup()
    
    def test_create_workflow(self):
        """Test creating a new workflow."""
        # Call the method
        workflow_state = self.service.create_workflow("What genes are differentially expressed in cancer?")
        
        # Check the result
        self.assertIsNotNone(workflow_state)
        self.assertTrue(workflow_state.workflow_id.startswith("workflow_"))
        self.assertEqual(workflow_state.research_question, "What genes are differentially expressed in cancer?")
        self.assertEqual(workflow_state.state, WorkflowState.INITIALIZED)
        
        # Check if the file exists
        cache_file = Path(self.temp_dir.name) / f"{workflow_state.workflow_id}_state.json"
        self.assertTrue(cache_file.exists())
    
    def test_search_datasets(self):
        """Test searching for datasets."""
        # Create a workflow state
        workflow_state = UORCAWorkflowState(
            workflow_id="workflow_test",
            state=WorkflowState.INITIALIZED,
            research_question="What genes are differentially expressed in cancer?"
        )
        
        # Mock the LLM service method
        self.service.llm_service.generate_search_queries.return_value = [
            {
                "term": "cancer gene expression",
                "description": "Search for gene expression data in cancer"
            },
            {
                "term": "tumor transcriptomics",
                "description": "Search for transcriptomics data in tumors"
            }
        ]
        
        # Mock the GEO service method
        self.service.geo_service.search_datasets.side_effect = [
            GEOSearchResult(
                query=GEOSearchQuery(term="cancer gene expression", max_results=3),
                total_count=2,
                dataset_ids=["GDS123", "GDS456"],
                series_ids=[]
            ),
            GEOSearchResult(
                query=GEOSearchQuery(term="tumor transcriptomics", max_results=3),
                total_count=1,
                dataset_ids=["GDS789"],
                series_ids=[]
            )
        ]
        
        # Call the method
        updated_state = self.service._search_datasets(workflow_state)
        
        # Check the result
        self.assertEqual(updated_state.state, WorkflowState.SEARCHING_DATASETS)
        self.assertEqual(len(updated_state.search_queries), 2)
        self.assertEqual(len(updated_state.search_results), 2)
        self.assertEqual(updated_state.search_results[0].dataset_ids, ["GDS123", "GDS456"])
        self.assertEqual(updated_state.search_results[1].dataset_ids, ["GDS789"])
    
    def test_evaluate_datasets(self):
        """Test evaluating dataset relevance."""
        # Create a workflow state with search results
        workflow_state = UORCAWorkflowState(
            workflow_id="workflow_test",
            state=WorkflowState.SEARCHING_DATASETS,
            research_question="What genes are differentially expressed in cancer?",
            search_results=[
                GEOSearchResult(
                    query=GEOSearchQuery(term="cancer gene expression", max_results=3),
                    total_count=2,
                    dataset_ids=["GDS123", "GDS456"],
                    series_ids=[]
                )
            ]
        )
        
        # Mock the GEO service methods
        platform = GEOPlatform(
            id="GPL123",
            title="Test Platform",
            technology="Microarray"
        )
        
        dataset1 = GEODataset(
            id="GDS123",
            title="Cancer Gene Expression",
            summary="Gene expression data from cancer samples",
            organism="human",
            platform=platform,
            samples=["GSM123", "GSM456"],
            omics_type=OmicsType.TRANSCRIPTOMICS
        )
        
        dataset2 = GEODataset(
            id="GDS456",
            title="Unrelated Dataset",
            summary="Gene expression data from unrelated condition",
            organism="mouse",
            platform=platform,
            samples=["GSM789", "GSM012"],
            omics_type=OmicsType.TRANSCRIPTOMICS
        )
        
        self.service.geo_service.get_dataset.side_effect = [dataset1, dataset2]
        
        # Mock the LLM service method
        self.service.llm_service.evaluate_dataset_relevance.side_effect = [
            DatasetRelevanceScore(
                dataset_id="GDS123",
                relevance_score=0.9,
                justification="Highly relevant to cancer gene expression"
            ),
            DatasetRelevanceScore(
                dataset_id="GDS456",
                relevance_score=0.3,
                justification="Not directly related to the research question"
            )
        ]
        
        # Call the method
        updated_state = self.service._evaluate_datasets(workflow_state)
        
        # Check the result
        self.assertEqual(updated_state.state, WorkflowState.EVALUATING_RELEVANCE)
        self.assertEqual(len(updated_state.dataset_relevance_scores), 2)
        self.assertEqual(updated_state.dataset_relevance_scores["GDS123"].relevance_score, 0.9)
        self.assertEqual(updated_state.dataset_relevance_scores["GDS456"].relevance_score, 0.3)
        self.assertEqual(len(updated_state.selected_datasets), 1)
        self.assertEqual(updated_state.selected_datasets[0], "GDS123")
    
    def test_fetch_datasets(self):
        """Test fetching datasets."""
        # Create a workflow state with selected datasets
        workflow_state = UORCAWorkflowState(
            workflow_id="workflow_test",
            state=WorkflowState.EVALUATING_RELEVANCE,
            research_question="What genes are differentially expressed in cancer?",
            selected_datasets=["GDS123"]
        )
        
        # Mock the GEO service method
        platform = GEOPlatform(
            id="GPL123",
            title="Test Platform",
            technology="Microarray"
        )
        
        metadata = GEODataset(
            id="GDS123",
            title="Cancer Gene Expression",
            summary="Gene expression data from cancer samples",
            organism="human",
            platform=platform,
            samples=["GSM123", "GSM456"],
            omics_type=OmicsType.TRANSCRIPTOMICS
        )
        
        expression_matrix = ExpressionMatrix(
            gene_ids=["GENE1", "GENE2", "GENE3"],
            sample_ids=["SAMPLE1", "SAMPLE2"],
            values=[
                [1.0, 2.0],
                [3.0, 4.0],
                [5.0, 6.0]
            ]
        )
        
        sample_metadata = {
            "GSM123": {"group": "control", "tissue": "normal"},
            "GSM456": {"group": "treatment", "tissue": "tumor"}
        }
        
        dataset = GEODatasetWithExpression(
            metadata=metadata,
            expression_matrix=expression_matrix,
            sample_metadata=sample_metadata
        )
        
        self.service.geo_service.get_dataset_with_expression.return_value = dataset
        
        # Call the method
        updated_state = self.service._fetch_datasets(workflow_state)
        
        # Check the result
        self.assertEqual(updated_state.state, WorkflowState.FETCHING_DATA)
        self.assertEqual(len(updated_state.fetched_datasets), 1)
        self.assertEqual(updated_state.fetched_datasets["GDS123"].metadata.id, "GDS123")
        self.assertEqual(len(updated_state.fetched_datasets["GDS123"].expression_matrix.gene_ids), 3)
    
    def test_analyze_datasets(self):
        """Test analyzing datasets."""
        # Create a workflow state with fetched datasets
        platform = GEOPlatform(
            id="GPL123",
            title="Test Platform",
            technology="Microarray"
        )
        
        metadata = GEODataset(
            id="GDS123",
            title="Cancer Gene Expression",
            summary="Gene expression data from cancer samples",
            organism="human",
            platform=platform,
            samples=["GSM123", "GSM456"],
            omics_type=OmicsType.TRANSCRIPTOMICS
        )
        
        expression_matrix = ExpressionMatrix(
            gene_ids=["GENE1", "GENE2", "GENE3"],
            sample_ids=["SAMPLE1", "SAMPLE2"],
            values=[
                [1.0, 2.0],
                [3.0, 4.0],
                [5.0, 6.0]
            ]
        )
        
        sample_metadata = {
            "GSM123": {"group": "control", "tissue": "normal"},
            "GSM456": {"group": "treatment", "tissue": "tumor"}
        }
        
        dataset = GEODatasetWithExpression(
            metadata=metadata,
            expression_matrix=expression_matrix,
            sample_metadata=sample_metadata
        )
        
        workflow_state = UORCAWorkflowState(
            workflow_id="workflow_test",
            state=WorkflowState.FETCHING_DATA,
            research_question="What genes are differentially expressed in cancer?",
            fetched_datasets={"GDS123": dataset}
        )
        
        # Mock the analysis service method
        analysis_result = AnalysisResult(
            request_id="analysis_123",
            status=AnalysisStatus.COMPLETED,
            research_question="What genes are differentially expressed in cancer?",
            datasets_analyzed=["GDS123"],
            omics_types=[OmicsType.TRANSCRIPTOMICS]
        )
        
        self.service.analysis_service.analyze_dataset.return_value = analysis_result
        
        # Call the method
        updated_state = self.service._analyze_datasets(workflow_state)
        
        # Check the result
        self.assertEqual(updated_state.state, WorkflowState.ANALYZING_DATA)
        self.assertEqual(len(updated_state.analysis_requests), 1)
        self.assertEqual(len(updated_state.analysis_results), 1)
        self.assertEqual(updated_state.analysis_results["GDS123"].request_id, "analysis_123")
        self.assertEqual(updated_state.analysis_results["GDS123"].status, AnalysisStatus.COMPLETED)
    
    def test_integrate_results(self):
        """Test integrating results."""
        # Create a workflow state with analysis results
        analysis_result1 = AnalysisResult(
            request_id="analysis_123",
            status=AnalysisStatus.COMPLETED,
            research_question="What genes are differentially expressed in cancer?",
            datasets_analyzed=["GDS123"],
            omics_types=[OmicsType.TRANSCRIPTOMICS]
        )
        
        analysis_result2 = AnalysisResult(
            request_id="analysis_456",
            status=AnalysisStatus.COMPLETED,
            research_question="What genes are differentially expressed in cancer?",
            datasets_analyzed=["GDS456"],
            omics_types=[OmicsType.TRANSCRIPTOMICS]
        )
        
        workflow_state = UORCAWorkflowState(
            workflow_id="workflow_test",
            state=WorkflowState.ANALYZING_DATA,
            research_question="What genes are differentially expressed in cancer?",
            analysis_results={
                "GDS123": analysis_result1,
                "GDS456": analysis_result2
            }
        )
        
        # Mock the analysis service method
        insights = [
            Insight(
                insight_type=InsightType.KEY_GENE,
                description="GENE3 is consistently upregulated in cancer",
                supporting_datasets=["GDS123", "GDS456"],
                confidence_score=0.9,
                related_entities=["GENE3", "cancer"]
            ),
            Insight(
                insight_type=InsightType.PATHWAY,
                description="Cancer pathway is enriched in both datasets",
                supporting_datasets=["GDS123", "GDS456"],
                confidence_score=0.8,
                related_entities=["cancer pathway"]
            )
        ]
        
        self.service.analysis_service.integrate_results.return_value = insights
        
        # Call the method
        updated_state = self.service._integrate_results(workflow_state)
        
        # Check the result
        self.assertEqual(updated_state.state, WorkflowState.GENERATING_INSIGHTS)
        self.assertEqual(len(updated_state.integrated_insights), 2)
        self.assertEqual(updated_state.integrated_insights[0].insight_type, InsightType.KEY_GENE)
        self.assertEqual(updated_state.integrated_insights[0].description, 
                         "GENE3 is consistently upregulated in cancer")
        self.assertEqual(updated_state.integrated_insights[1].insight_type, InsightType.PATHWAY)
    
    def test_execute_workflow(self):
        """Test executing the complete workflow."""
        # Create a workflow state
        workflow_state = UORCAWorkflowState(
            workflow_id="workflow_test",
            state=WorkflowState.INITIALIZED,
            research_question="What genes are differentially expressed in cancer?"
        )
        
        # Mock the individual workflow steps
        def mock_search_datasets(state):
            state.state = WorkflowState.SEARCHING_DATASETS
            state.search_queries = [GEOSearchQuery(term="cancer", max_results=3)]
            state.search_results = [
                GEOSearchResult(
                    query=GEOSearchQuery(term="cancer", max_results=3),
                    total_count=1,
                    dataset_ids=["GDS123"],
                    series_ids=[]
                )
            ]
            return state
        
        def mock_evaluate_datasets(state):
            state.state = WorkflowState.EVALUATING_RELEVANCE
            state.dataset_relevance_scores = {
                "GDS123": DatasetRelevanceScore(
                    dataset_id="GDS123",
                    relevance_score=0.9,
                    justification="Highly relevant"
                )
            }
            state.selected_datasets = ["GDS123"]
            return state
        
        def mock_fetch_datasets(state):
            state.state = WorkflowState.FETCHING_DATA
            platform = GEOPlatform(
                id="GPL123",
                title="Test Platform",
                technology="Microarray"
            )
            
            metadata = GEODataset(
                id="GDS123",
                title="Cancer Gene Expression",
                summary="Gene expression data from cancer samples",
                organism="human",
                platform=platform,
                samples=["GSM123", "GSM456"],
                omics_type=OmicsType.TRANSCRIPTOMICS
            )
            
            expression_matrix = ExpressionMatrix(
                gene_ids=["GENE1", "GENE2", "GENE3"],
                sample_ids=["SAMPLE1", "SAMPLE2"],
                values=[
                    [1.0, 2.0],
                    [3.0, 4.0],
                    [5.0, 6.0]
                ]
            )
            
            state.fetched_datasets = {
                "GDS123": GEODatasetWithExpression(
                    metadata=metadata,
                    expression_matrix=expression_matrix,
                    sample_metadata={}
                )
            }
            return state
        
        def mock_analyze_datasets(state):
            state.state = WorkflowState.ANALYZING_DATA
            state.analysis_results = {
                "GDS123": AnalysisResult(
                    request_id="analysis_123",
                    status=AnalysisStatus.COMPLETED,
                    research_question="What genes are differentially expressed in cancer?",
                    datasets_analyzed=["GDS123"],
                    omics_types=[OmicsType.TRANSCRIPTOMICS]
                )
            }
            return state
        
        def mock_integrate_results(state):
            state.state = WorkflowState.INTEGRATING_RESULTS
            state.integrated_insights = [
                Insight(
                    insight_type=InsightType.KEY_GENE,
                    description="GENE3 is upregulated in cancer",
                    supporting_datasets=["GDS123"],
                    confidence_score=0.9,
                    related_entities=["GENE3", "cancer"]
                )
            ]
            return state
        
        # Patch the individual workflow steps
        self.service._search_datasets = MagicMock(side_effect=mock_search_datasets)
        self.service._evaluate_datasets = MagicMock(side_effect=mock_evaluate_datasets)
        self.service._fetch_datasets = MagicMock(side_effect=mock_fetch_datasets)
        self.service._analyze_datasets = MagicMock(side_effect=mock_analyze_datasets)
        self.service._integrate_results = MagicMock(side_effect=mock_integrate_results)
        
        # Call the method
        updated_state = self.service.execute_workflow(workflow_state)
        
        # Check the result
        self.assertEqual(updated_state.state, WorkflowState.COMPLETED)
        self.assertTrue(self.service._search_datasets.called)
        self.assertTrue(self.service._evaluate_datasets.called)
        self.assertTrue(self.service._fetch_datasets.called)
        self.assertTrue(self.service._analyze_datasets.called)
        self.assertTrue(self.service._integrate_results.called)
    
    def test_get_workflow_summary(self):
        """Test generating a workflow summary."""
        # Create a completed workflow state
        workflow_state = UORCAWorkflowState(
            workflow_id="workflow_test",
            state=WorkflowState.COMPLETED,
            research_question="What genes are differentially expressed in cancer?",
            analysis_results={
                "GDS123": AnalysisResult(
                    request_id="analysis_123",
                    status=AnalysisStatus.COMPLETED,
                    research_question="What genes are differentially expressed in cancer?",
                    datasets_analyzed=["GDS123"],
                    omics_types=[OmicsType.TRANSCRIPTOMICS]
                )
            },
            integrated_insights=[
                Insight(
                    insight_type=InsightType.KEY_GENE,
                    description="GENE3 is upregulated in cancer",
                    supporting_datasets=["GDS123"],
                    confidence_score=0.9,
                    related_entities=["GENE3", "cancer"]
                )
            ]
        )
        
        # Mock the analysis service method
        self.service.analysis_service.generate_summary.return_value = """
# Analysis Summary

## Research Question
What genes are differentially expressed in cancer?

## Datasets Analyzed
- GDS123

## Key Findings
1. GENE3 is upregulated in cancer (high confidence)

## Implications
These findings suggest that GENE3 may play a role in cancer development.
"""
        
        # Call the method
        summary = self.service.get_workflow_summary(workflow_state)
        
        # Check the result
        self.assertTrue("GENE3 is upregulated in cancer" in summary)
        self.assertTrue("may play a role in cancer development" in summary)
    
    def test_load_workflow_state(self):
        """Test loading a workflow state from disk."""
        # Create and save a workflow state
        workflow_state = UORCAWorkflowState(
            workflow_id="workflow_test",
            state=WorkflowState.COMPLETED,
            research_question="What genes are differentially expressed in cancer?"
        )
        
        self.service._save_workflow_state(workflow_state)
        
        # Call the method
        loaded_state = self.service.load_workflow_state("workflow_test")
        
        # Check the result
        self.assertIsNotNone(loaded_state)
        self.assertEqual(loaded_state.workflow_id, "workflow_test")
        self.assertEqual(loaded_state.research_question, "What genes are differentially expressed in cancer?")
        self.assertEqual(loaded_state.state, WorkflowState.COMPLETED)

if __name__ == '__main__':
    unittest.main()

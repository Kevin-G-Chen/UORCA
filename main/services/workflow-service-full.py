"""
Service for coordinating the agentic workflow.
"""
import logging
import uuid
from typing import Dict, List, Any, Optional, Union
from pathlib import Path
import json

from ..config import get_settings
from ..models.workflow_models import (
    UORCAWorkflowState,
    WorkflowState,
    DatasetRelevanceScore,
    DatasetSearchTask,
    DatasetEvaluationTask,
    DataAnalysisTask,
    InsightIntegrationTask,
    DatasetSearchAgent,
    DataAnalysisAgent,
    InsightIntegrationAgent
)
from ..models.geo_models import (
    GEOSearchQuery,
    GEOSearchResult,
    GEODatasetWithExpression,
    OmicsType
)
from ..models.analysis_models import (
    AnalysisRequest,
    AnalysisResult,
    AnalysisType,
    AnalysisStatus,
    Insight
)
from ..services.geo_service import GEOService
from ..services.llm_service import LLMService
from ..services.analysis_service import AnalysisService

logger = logging.getLogger(__name__)

class WorkflowService:
    """Service for coordinating the UORCA workflow."""
    
    def __init__(self):
        """Initialize the workflow service."""
        self.settings = get_settings()
        self.geo_service = GEOService()
        self.llm_service = LLMService()
        self.analysis_service = AnalysisService()
        
        # Initialize agents
        self.dataset_search_agent = DatasetSearchAgent()
        self.data_analysis_agent = DataAnalysisAgent()
        self.insight_integration_agent = InsightIntegrationAgent()
    
    def create_workflow(self, research_question: str) -> UORCAWorkflowState:
        """
        Create a new workflow for a research question.
        
        Args:
            research_question: The research question to address
            
        Returns:
            The initialized workflow state
        """
        workflow_id = f"workflow_{uuid.uuid4().hex[:8]}"
        
        # Create the workflow state
        workflow_state = UORCAWorkflowState(
            workflow_id=workflow_id,
            state=WorkflowState.INITIALIZED,
            research_question=research_question
        )
        
        # Save the initial state
        self._save_workflow_state(workflow_state)
        
        logger.info(f"Created new workflow {workflow_id} for research question: {research_question}")
        
        return workflow_state
    
    def execute_workflow(self, workflow_state: UORCAWorkflowState) -> UORCAWorkflowState:
        """
        Execute the complete workflow.
        
        Args:
            workflow_state: The workflow state
            
        Returns:
            The updated workflow state
        """
        logger.info(f"Executing workflow {workflow_state.workflow_id}")
        
        try:
            # Step 1: Search for datasets
            workflow_state = self._search_datasets(workflow_state)
            if workflow_state.state == WorkflowState.FAILED:
                return workflow_state
            
            # Step 2: Evaluate dataset relevance
            workflow_state = self._evaluate_datasets(workflow_state)
            if workflow_state.state == WorkflowState.FAILED:
                return workflow_state
            
            # Step 3: Fetch selected datasets
            workflow_state = self._fetch_datasets(workflow_state)
            if workflow_state.state == WorkflowState.FAILED:
                return workflow_state
            
            # Step 4: Analyze datasets
            workflow_state = self._analyze_datasets(workflow_state)
            if workflow_state.state == WorkflowState.FAILED:
                return workflow_state
            
            # Step 5: Integrate results
            workflow_state = self._integrate_results(workflow_state)
            if workflow_state.state == WorkflowState.FAILED:
                return workflow_state
            
            # Mark as completed
            workflow_state.state = WorkflowState.COMPLETED
            self._save_workflow_state(workflow_state)
            
            logger.info(f"Workflow {workflow_state.workflow_id} completed successfully")
            
            return workflow_state
            
        except Exception as e:
            logger.error(f"Error executing workflow {workflow_state.workflow_id}: {str(e)}")
            workflow_state.add_error(f"Workflow execution error: {str(e)}")
            self._save_workflow_state(workflow_state)
            return workflow_state
    
    def _search_datasets(self, workflow_state: UORCAWorkflowState) -> UORCAWorkflowState:
        """
        Search for datasets relevant to the research question.
        
        Args:
            workflow_state: The workflow state
            
        Returns:
            The updated workflow state
        """
        logger.info(f"Searching for datasets for workflow {workflow_state.workflow_id}")
        
        try:
            workflow_state.state = WorkflowState.SEARCHING_DATASETS
            self._save_workflow_state(workflow_state)
            
            # Generate search queries using LLM
            search_queries_data = self.llm_service.generate_search_queries(
                research_question=workflow_state.research_question,
                max_queries=3  # Limit to 3 queries for testing
            )
            
            # Convert to GEOSearchQuery objects
            for query_data in search_queries_data:
                query = GEOSearchQuery(
                    term=query_data.get("term", ""),
                    max_results=self.settings.max_datasets_to_analyze,
                    omics_type=None,  # We'll let the search term handle this
                    organism=None  # We'll let the search term handle this
                )
                workflow_state.search_queries.append(query)
            
            # Execute each search query
            for query in workflow_state.search_queries:
                try:
                    result = self.geo_service.search_datasets(query)
                    workflow_state.search_results.append(result)
                except Exception as e:
                    logger.error(f"Error executing search query '{query.term}': {str(e)}")
                    workflow_state.add_warning(f"Error executing search query '{query.term}': {str(e)}")
            
            # Check if we found any datasets
            if not any(not result.is_empty() for result in workflow_state.search_results):
                workflow_state.add_warning("No datasets found for any search query")
            
            self._save_workflow_state(workflow_state)
            
            return workflow_state
            
        except Exception as e:
            logger.error(f"Error searching datasets: {str(e)}")
            workflow_state.add_error(f"Dataset search error: {str(e)}")
            self._save_workflow_state(workflow_state)
            return workflow_state
    
    def _evaluate_datasets(self, workflow_state: UORCAWorkflowState) -> UORCAWorkflowState:
        """
        Evaluate the relevance of datasets to the research question.
        
        Args:
            workflow_state: The workflow state
            
        Returns:
            The updated workflow state
        """
        logger.info(f"Evaluating datasets for workflow {workflow_state.workflow_id}")
        
        try:
            workflow_state.state = WorkflowState.EVALUATING_RELEVANCE
            self._save_workflow_state(workflow_state)
            
            # Collect all dataset IDs
            all_dataset_ids = []
            for result in workflow_state.search_results:
                all_dataset_ids.extend(result.dataset_ids)
                # For simplicity, we'll only evaluate datasets (GDS), not series (GSE)
                # all_dataset_ids.extend(result.series_ids)
            
            # Remove duplicates
            all_dataset_ids = list(set(all_dataset_ids))
            
            # Limit the number of datasets to evaluate
            max_datasets = min(len(all_dataset_ids), self.settings.max_datasets_to_analyze)
            datasets_to_evaluate = all_dataset_ids[:max_datasets]
            
            # Fetch metadata for each dataset
            for dataset_id in datasets_to_evaluate:
                try:
                    # Fetch dataset metadata
                    if dataset_id.startswith("GDS"):
                        dataset_metadata = self.geo_service.get_dataset(dataset_id)
                    else:
                        dataset_metadata = self.geo_service.get_series(dataset_id)
                    
                    # Evaluate relevance
                    relevance_score = self.llm_service.evaluate_dataset_relevance(
                        research_question=workflow_state.research_question,
                        dataset_metadata=dataset_metadata.dict()
                    )
                    
                    # Store the relevance score
                    workflow_state.dataset_relevance_scores[dataset_id] = relevance_score
                    
                except Exception as e:
                    logger.error(f"Error evaluating dataset {dataset_id}: {str(e)}")
                    workflow_state.add_warning(f"Error evaluating dataset {dataset_id}: {str(e)}")
            
            # Select datasets with highest relevance scores
            if workflow_state.dataset_relevance_scores:
                # Sort by relevance score in descending order
                sorted_datasets = sorted(
                    workflow_state.dataset_relevance_scores.items(),
                    key=lambda x: x[1].relevance_score,
                    reverse=True
                )
                
                # Select top datasets
                max_to_select = min(len(sorted_datasets), self.settings.max_datasets_to_analyze)
                for i in range(max_to_select):
                    dataset_id, score = sorted_datasets[i]
                    # Only include datasets with relevance score above threshold
                    if score.relevance_score >= 0.5:
                        workflow_state.selected_datasets.append(dataset_id)
            
            # Check if we selected any datasets
            if not workflow_state.selected_datasets:
                workflow_state.add_warning("No relevant datasets found above threshold")
            
            self._save_workflow_state(workflow_state)
            
            return workflow_state
            
        except Exception as e:
            logger.error(f"Error evaluating datasets: {str(e)}")
            workflow_state.add_error(f"Dataset evaluation error: {str(e)}")
            self._save_workflow_state(workflow_state)
            return workflow_state
    
    def _fetch_datasets(self, workflow_state: UORCAWorkflowState) -> UORCAWorkflowState:
        """
        Fetch selected datasets with expression data.
        
        Args:
            workflow_state: The workflow state
            
        Returns:
            The updated workflow state
        """
        logger.info(f"Fetching datasets for workflow {workflow_state.workflow_id}")
        
        try:
            workflow_state.state = WorkflowState.FETCHING_DATA
            self._save_workflow_state(workflow_state)
            
            # Fetch each selected dataset
            for dataset_id in workflow_state.selected_datasets:
                try:
                    # Fetch dataset with expression data
                    if dataset_id.startswith("GDS"):
                        dataset = self.geo_service.get_dataset_with_expression(dataset_id)
                    else:
                        dataset = self.geo_service.get_series_with_expression(dataset_id)
                    
                    # Store the dataset
                    workflow_state.fetched_datasets[dataset_id] = dataset
                    
                except Exception as e:
                    logger.error(f"Error fetching dataset {dataset_id}: {str(e)}")
                    workflow_state.add_warning(f"Error fetching dataset {dataset_id}: {str(e)}")
            
            # Check if we fetched any datasets
            if not workflow_state.fetched_datasets:
                workflow_state.add_error("Failed to fetch any datasets")
                return workflow_state
            
            self._save_workflow_state(workflow_state)
            
            return workflow_state
            
        except Exception as e:
            logger.error(f"Error fetching datasets: {str(e)}")
            workflow_state.add_error(f"Dataset fetching error: {str(e)}")
            self._save_workflow_state(workflow_state)
            return workflow_state
    
    def _analyze_datasets(self, workflow_state: UORCAWorkflowState) -> UORCAWorkflowState:
        """
        Analyze fetched datasets.
        
        Args:
            workflow_state: The workflow state
            
        Returns:
            The updated workflow state
        """
        logger.info(f"Analyzing datasets for workflow {workflow_state.workflow_id}")
        
        try:
            workflow_state.state = WorkflowState.ANALYZING_DATA
            self._save_workflow_state(workflow_state)
            
            # Analyze each fetched dataset
            for dataset_id, dataset in workflow_state.fetched_datasets.items():
                try:
                    # Skip datasets without expression data
                    if not dataset.expression_matrix:
                        workflow_state.add_warning(f"Dataset {dataset_id} has no expression data, skipping analysis")
                        continue
                    
                    # Create analysis request
                    request = AnalysisRequest(
                        research_question=workflow_state.research_question,
                        analysis_type=AnalysisType.CUSTOM,  # Let the LLM determine the appropriate type
                        dataset_ids=[dataset_id],
                        max_datasets=1
                    )
                    
                    # Store the request
                    workflow_state.analysis_requests[dataset_id] = request
                    
                    # Perform the analysis
                    result = self.analysis_service.analyze_dataset(
                        dataset=dataset,
                        request=request
                    )
                    
                    # Store the result
                    workflow_state.analysis_results[dataset_id] = result
                    
                except Exception as e:
                    logger.error(f"Error analyzing dataset {dataset_id}: {str(e)}")
                    workflow_state.add_warning(f"Error analyzing dataset {dataset_id}: {str(e)}")
            
            # Check if we analyzed any datasets
            if not workflow_state.analysis_results:
                workflow_state.add_error("Failed to analyze any datasets")
                return workflow_state
            
            self._save_workflow_state(workflow_state)
            
            return workflow_state
            
        except Exception as e:
            logger.error(f"Error analyzing datasets: {str(e)}")
            workflow_state.add_error(f"Dataset analysis error: {str(e)}")
            self._save_workflow_state(workflow_state)
            return workflow_state
    
    def _integrate_results(self, workflow_state: UORCAWorkflowState) -> UORCAWorkflowState:
        """
        Integrate results from all analyzed datasets.
        
        Args:
            workflow_state: The workflow state
            
        Returns:
            The updated workflow state
        """
        logger.info(f"Integrating results for workflow {workflow_state.workflow_id}")
        
        try:
            workflow_state.state = WorkflowState.INTEGRATING_RESULTS
            self._save_workflow_state(workflow_state)
            
            # Get all successful analysis results
            successful_results = [
                result for result in workflow_state.analysis_results.values()
                if result.status == AnalysisStatus.COMPLETED
            ]
            
            # Skip integration if we have fewer than 2 successful results
            if len(successful_results) < 2:
                workflow_state.add_warning("Not enough successful analyses to perform integration")
                workflow_state.state = WorkflowState.GENERATING_INSIGHTS
                self._save_workflow_state(workflow_state)
                return workflow_state
            
            # Integrate the results
            insights = self.analysis_service.integrate_results(
                research_question=workflow_state.research_question,
                analysis_results=successful_results
            )
            
            # Store the insights
            workflow_state.integrated_insights = insights
            
            workflow_state.state = WorkflowState.GENERATING_INSIGHTS
            self._save_workflow_state(workflow_state)
            
            return workflow_state
            
        except Exception as e:
            logger.error(f"Error integrating results: {str(e)}")
            workflow_state.add_error(f"Result integration error: {str(e)}")
            self._save_workflow_state(workflow_state)
            return workflow_state
    
    def get_workflow_summary(self, workflow_state: UORCAWorkflowState) -> str:
        """
        Generate a human-readable summary of the workflow results.
        
        Args:
            workflow_state: The workflow state
            
        Returns:
            Formatted summary
        """
        # Check if the workflow has completed
        if workflow_state.state != WorkflowState.COMPLETED:
            return f"Workflow is not yet completed. Current state: {workflow_state.state.value}"
        
        # Get datasets analyzed
        datasets_analyzed = list(workflow_state.analysis_results.keys())
        
        # Generate the summary
        return self.analysis_service.generate_summary(
            research_question=workflow_state.research_question,
            datasets_analyzed=datasets_analyzed,
            insights=workflow_state.integrated_insights
        )
    
    def _save_workflow_state(self, workflow_state: UORCAWorkflowState) -> None:
        """
        Save the current workflow state to disk.
        
        Args:
            workflow_state: The workflow state to save
        """
        try:
            # Create a JSON-serializable representation
            state_dict = workflow_state.dict()
            
            # Save to the cache directory
            cache_dir = Path(self.settings.analysis_cache_dir)
            state_file = cache_dir / f"{workflow_state.workflow_id}_state.json"
            
            with open(state_file, 'w') as f:
                json.dump(state_dict, f, indent=2)
                
            logger.debug(f"Saved workflow state to {state_file}")
        except Exception as e:
            logger.error(f"Error saving workflow state: {str(e)}")
    
    def load_workflow_state(self, workflow_id: str) -> Optional[UORCAWorkflowState]:
        """
        Load a workflow state from disk.
        
        Args:
            workflow_id: The workflow ID
            
        Returns:
            The loaded workflow state, or None if not found
        """
        try:
            # Load from the cache directory
            cache_dir = Path(self.settings.analysis_cache_dir)
            state_file = cache_dir / f"{workflow_id}_state.json"
            
            if not state_file.exists():
                logger.warning(f"Workflow state file not found: {state_file}")
                return None
            
            with open(state_file, 'r') as f:
                state_dict = json.load(f)
            
            # Recreate the workflow state
            workflow_state = UORCAWorkflowState(**state_dict)
            
            logger.debug(f"Loaded workflow state from {state_file}")
            return workflow_state
            
        except Exception as e:
            logger.error(f"Error loading workflow state: {str(e)}")
            return None
"""
Pydantic models for the agentic workflow state management.
"""
from typing import Dict, List, Optional, Any, Union
from enum import Enum
from pydantic import BaseModel, Field, validator
from pydanticai import Agent, Task, LLM, Human, Step

from .geo_models import GEOSearchQuery, GEOSearchResult, GEODatasetWithExpression, OmicsType
from .analysis_models import AnalysisRequest, AnalysisResult, AnalysisStatus, Insight


class WorkflowState(str, Enum):
    """State of the UORCA workflow."""
    INITIALIZED = "initialized"
    SEARCHING_DATASETS = "searching_datasets"
    EVALUATING_RELEVANCE = "evaluating_relevance"
    FETCHING_DATA = "fetching_data"
    ANALYZING_DATA = "analyzing_data"
    INTEGRATING_RESULTS = "integrating_results"
    GENERATING_INSIGHTS = "generating_insights"
    COMPLETED = "completed"
    FAILED = "failed"


class DatasetRelevanceScore(BaseModel):
    """Relevance score for a dataset to a research question."""
    dataset_id: str = Field(..., description="Dataset ID")
    relevance_score: float = Field(..., description="Relevance score (0-1)")
    justification: str = Field(..., description="Justification for the relevance score")
    

class UORCAWorkflowState(BaseModel):
    """State of the UORCA workflow."""
    workflow_id: str = Field(..., description="Unique workflow identifier")
    state: WorkflowState = Field(WorkflowState.INITIALIZED, description="Current state of the workflow")
    research_question: str = Field(..., description="Research question being addressed")
    
    # Dataset search and selection
    search_queries: List[GEOSearchQuery] = Field(default_factory=list, description="Search queries generated for the research question")
    search_results: List[GEOSearchResult] = Field(default_factory=list, description="Results of dataset searches")
    dataset_relevance_scores: Dict[str, DatasetRelevanceScore] = Field(default_factory=dict, description="Relevance scores for datasets")
    selected_datasets: List[str] = Field(default_factory=list, description="IDs of selected datasets for analysis")
    
    # Fetched datasets
    fetched_datasets: Dict[str, GEODatasetWithExpression] = Field(default_factory=dict, description="Fetched datasets with expression data")
    
    # Analysis
    analysis_requests: Dict[str, AnalysisRequest] = Field(default_factory=dict, description="Analysis requests by dataset ID")
    analysis_results: Dict[str, AnalysisResult] = Field(default_factory=dict, description="Results of the analyses")
    
    # Integration
    integrated_insights: List[Insight] = Field(default_factory=list, description="Integrated insights from all analyses")
    
    # Errors and warnings
    errors: List[str] = Field(default_factory=list, description="Error messages")
    warnings: List[str] = Field(default_factory=list, description="Warning messages")
    
    def add_error(self, error_message: str) -> None:
        """Add an error message to the workflow state."""
        self.errors.append(error_message)
        if self.state != WorkflowState.FAILED:
            self.state = WorkflowState.FAILED
    
    def add_warning(self, warning_message: str) -> None:
        """Add a warning message to the workflow state."""
        self.warnings.append(warning_message)


# PydanticAI models for the agentic workflow

class DatasetSearchTask(Task):
    """Task to search for relevant datasets based on a research question."""
    research_question: str = Field(..., description="Research question to address")
    target_organism: Optional[str] = Field(None, description="Target organism")
    omics_types: List[OmicsType] = Field(default_factory=list, description="Types of omics data to consider")
    

class DatasetEvaluationTask(Task):
    """Task to evaluate the relevance of datasets to a research question."""
    research_question: str = Field(..., description="Research question to address")
    dataset_metadata: List[Dict[str, Any]] = Field(..., description="Metadata for datasets to evaluate")
    

class DataAnalysisTask(Task):
    """Task to design and run analysis on a dataset."""
    research_question: str = Field(..., description="Research question to address")
    dataset: GEODatasetWithExpression = Field(..., description="Dataset to analyze")
    

class InsightIntegrationTask(Task):
    """Task to integrate insights from multiple dataset analyses."""
    research_question: str = Field(..., description="Research question to address")
    analysis_results: List[AnalysisResult] = Field(..., description="Results of individual dataset analyses")
    

class DatasetSearchAgent(Agent):
    """Agent responsible for searching and selecting relevant datasets."""
    name: str = Field("DatasetSearchAgent", description="Name of the agent")
    description: str = Field(
        "I identify and evaluate datasets relevant to a research question by generating appropriate search queries "
        "and evaluating dataset metadata to determine relevance.",
        description="Agent description"
    )
    
    def search_datasets(self, task: DatasetSearchTask) -> List[GEOSearchQuery]:
        """Generate search queries for the research question."""
        pass
    
    def evaluate_relevance(self, task: DatasetEvaluationTask) -> Dict[str, DatasetRelevanceScore]:
        """Evaluate the relevance of datasets to the research question."""
        pass


class DataAnalysisAgent(Agent):
    """Agent responsible for designing and executing analyses on datasets."""
    name: str = Field("DataAnalysisAgent", description="Name of the agent")
    description: str = Field(
        "I design and execute appropriate analyses on datasets based on their content and the research question. "
        "I can parse metadata, design code for analysis, and produce meaningful results.",
        description="Agent description"
    )
    
    def design_analysis(self, task: DataAnalysisTask) -> AnalysisRequest:
        """Design an appropriate analysis for the dataset."""
        pass
    
    def execute_analysis(self, analysis_request: AnalysisRequest, dataset: GEODatasetWithExpression) -> AnalysisResult:
        """Execute the analysis on the dataset."""
        pass


class InsightIntegrationAgent(Agent):
    """Agent responsible for integrating insights across multiple datasets."""
    name: str = Field("InsightIntegrationAgent", description="Name of the agent")
    description: str = Field(
        "I integrate findings across multiple datasets to derive meaningful insights. "
        "I identify common patterns, contradictions, and novel biological trends.",
        description="Agent description"
    )
    
    def integrate_insights(self, task: InsightIntegrationTask) -> List[Insight]:
        """Integrate insights from multiple dataset analyses."""
        pass


class UORCAWorkflow(BaseModel):
    """The complete UORCA workflow model."""
    workflow_state: UORCAWorkflowState = Field(..., description="Current state of the workflow")
    dataset_search_agent: DatasetSearchAgent = Field(..., description="Agent for dataset search and selection")
    data_analysis_agent: DataAnalysisAgent = Field(..., description="Agent for data analysis")
    insight_integration_agent: InsightIntegrationAgent = Field(..., description="Agent for insight integration")
    
    # Define the workflow steps
    search_datasets_step: Step = Field(
        Step(
            description="Search for relevant datasets based on the research question",
            agent="dataset_search_agent",
            method="search_datasets"
        ),
        description="Step to search for datasets"
    )
    
    evaluate_datasets_step: Step = Field(
        Step(
            description="Evaluate the relevance of datasets to the research question",
            agent="dataset_search_agent",
            method="evaluate_relevance"
        ),
        description="Step to evaluate dataset relevance"
    )
    
    analyze_datasets_step: Step = Field(
        Step(
            description="Design and execute analyses on selected datasets",
            agent="data_analysis_agent",
            method="design_analysis"
        ),
        description="Step to analyze datasets"
    )
    
    integrate_insights_step: Step = Field(
        Step(
            description="Integrate insights from multiple dataset analyses",
            agent="insight_integration_agent",
            method="integrate_insights"
        ),
        description="Step to integrate insights"
    )

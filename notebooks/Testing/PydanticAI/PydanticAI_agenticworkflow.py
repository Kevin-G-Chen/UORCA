# %% Imports
import os
import sys
from typing import List, Dict, Any, Optional, Union, Literal
from pydantic import BaseModel, Field
import pandas as pd
import numpy as np
from enum import Enum

# Add these if you have PydanticAI installed
try:
    from pydantic_ai import Agent, RunContext
    from pydantic_ai.tools import Tool
except ImportError:
    print("PydanticAI not installed. Using mock classes for demonstration.")
    # Mock classes for demonstration
    class Agent:
        def __init__(self, name=None, description=None, deps_type=None, system_prompt=None):
            self.name = name
            self.description = description
            self.deps_type = deps_type
            self.system_prompt = system_prompt

        def tool(self, func=None):
            return func

        def tool_plain(self, func=None):
            return func

        def execute_workflow(self, *args, **kwargs):
            return None

    class RunContext:
        def __init__(self, deps=None):
            self.deps = deps

    class Tool:
        def __init__(self, func=None, name=None, description=None, takes_ctx=True):
            self.func = func
            self.name = name
            self.description = description
            self.takes_ctx = takes_ctx

        def __call__(self, *args, **kwargs):
            if self.func:
                return self.func(*args, **kwargs)
            return None

# %% Define Data Models
class GEODataset(BaseModel):
    """A model representing a GEO dataset."""
    accession: str
    title: str
    summary: str
    organism: str
    samples: int
    platform: str
    relevance_score: Optional[float] = None

class DataFile(BaseModel):
    """A model representing a data file from a GEO dataset."""
    dataset_accession: str
    file_id: str
    file_type: str
    file_url: str
    file_size: Optional[int] = None
    downloaded: bool = False
    file_path: Optional[str] = None

class AnalysisResult(BaseModel):
    """A model representing analysis results."""
    dataset_accession: str
    differentially_expressed_genes: List[str]
    enriched_pathways: List[str]
    visualizations: List[str]
    summary: str

class AnalysisStatus(str, Enum):
    SUCCESS = "success"
    FAILURE = "failure"
    IN_PROGRESS = "in_progress"

class WorkflowStep(str, Enum):
    DATASET_IDENTIFICATION = "dataset_identification"
    DATA_EXTRACTION = "data_extraction"
    DATA_ANALYSIS = "data_analysis"
    COMPLETED = "completed"

class WorkflowState(BaseModel):
    """A model representing the current state of the workflow."""
    query: str
    current_step: WorkflowStep = WorkflowStep.DATASET_IDENTIFICATION
    datasets: List[GEODataset] = []
    selected_dataset: Optional[GEODataset] = None
    data_files: List[DataFile] = []
    analysis_results: Optional[AnalysisResult] = None
    status: AnalysisStatus = AnalysisStatus.IN_PROGRESS
    error_message: Optional[str] = None

class WorkflowDependencies(BaseModel):
    """Dependencies shared across the workflow agents."""
    query: str
    temp_dir: str = "/tmp/rnaseq_data"

    # In a real implementation, this would include shared resources:
    # database_conn: Optional[DatabaseConn] = None  # Database connection
    # api_client: Optional[ApiClient] = None        # API client for external services
    # cache: Dict[str, Any] = Field(default_factory=dict)  # Shared cache across agents

    # You might add helper methods here that all agents could use
    def get_download_path(self, filename: str) -> str:
        """Get the full path for a file to be downloaded."""
        return os.path.join(self.temp_dir, filename)

    def ensure_temp_dir(self) -> None:
        """Ensure the temporary directory exists."""
        os.makedirs(self.temp_dir, exist_ok=True)

# %% Dataset Identification Agent
DatasetIdentificationAgent = Agent[None, List[GEODataset]](
    'openai:gpt-4o',  # Specify a model
    system_prompt=("You are a specialist in identifying relevant genomic datasets from GEO databases for research queries. "
        "Evaluate datasets based on their relevance to the research question and provide ranked results."),
    result_type=List[GEODataset]  # Specify the expected return type
)


@DatasetIdentificationAgent.tool
def query_geo_datasets(ctx: RunContext, query: str, max_results: int = 10) -> List[GEODataset]:
    """
    Query NCBI GEO for datasets matching the given query.

    Args:
        query: The research query.
        max_results: Maximum number of results to return.

    Returns:
        A list of GEODataset objects.
    """
    # In a real implementation, this would make API calls to NCBI GEO
    # For demonstration, return mock data
    mock_datasets = [
        GEODataset(
            accession=f"GSE1000{i}",
            title=f"RNA-seq analysis of {query.split()[0]} in human cells",
            summary=f"This dataset contains RNA-seq data from {query}.",
            organism="Homo sapiens",
            samples=6 + i,
            platform="Illumina HiSeq 2500"
        )
        for i in range(1, max_results+1)
    ]
    return mock_datasets

@DatasetIdentificationAgent.tool
def evaluate_dataset_relevance(ctx: RunContext, datasets: List[GEODataset], query: str) -> List[GEODataset]:
    """
    Evaluate the relevance of each dataset to the research query.

    Args:
        datasets: List of GEODataset objects.
        query: The research query.

    Returns:
        A list of GEODataset objects with relevance scores.
    """
    # In a real implementation, this would use NLP or other techniques to score relevance
    # For demonstration, assign random scores
    import random
    for dataset in datasets:
        dataset.relevance_score = random.uniform(0.5, 1.0)

    # Sort by relevance score
    return sorted(datasets, key=lambda x: x.relevance_score, reverse=True)

@DatasetIdentificationAgent.tool
def select_best_dataset(ctx: RunContext, datasets: List[GEODataset]) -> GEODataset:
    """
    Select the best dataset based on relevance score.

    Args:
        datasets: List of GEODataset objects with relevance scores.

    Returns:
        The best GEODataset object.
    """
    if not datasets:
        raise ValueError("No datasets provided")

    # Select the dataset with the highest relevance score
    return max(datasets, key=lambda x: x.relevance_score or 0)

# %% Data Extraction Agent
DataExtractionAgent = Agent[WorkflowDependencies, Dict[str, Any]](
    'openai:gpt-4o',  # Specify a model
    deps_type=WorkflowDependencies,  # Specify the dependency type
    result_type=Dict[str, Any],  # Specify the expected return type
    system_prompt="You are a specialist in extracting and processing genomic data files from GEO datasets. "
        "Extract FASTQ files, metadata, and prepare the data for analysis."
)

@DataExtractionAgent.tool
def extract_fastq_files(ctx: RunContext, dataset: GEODataset) -> List[DataFile]:
    """
    Extract FASTQ files from the given dataset.

    Args:
        dataset: A GEODataset object.

    Returns:
        A list of DataFile objects representing FASTQ files.
    """
    # In a real implementation, this would make API calls to NCBI GEO/SRA
    # For demonstration, return mock data
    fastq_files = [
        DataFile(
            dataset_accession=dataset.accession,
            file_id=f"{dataset.accession}_sample{i}_R{j}",
            file_type="fastq",
            file_url=f"https://example.com/{dataset.accession}_sample{i}_R{j}.fastq.gz",
            file_size=10000000 + i * 1000000 + j * 100000
        )
        for i in range(1, dataset.samples + 1) for j in range(1, 3)  # Paired-end reads
    ]
    return fastq_files

@DataExtractionAgent.tool
def extract_metadata(ctx: RunContext, dataset: GEODataset) -> Dict[str, Any]:
    """
    Extract metadata from the given dataset.

    Args:
        dataset: A GEODataset object.

    Returns:
        A dictionary containing metadata.
    """
    # In a real implementation, this would parse metadata from GEO
    # For demonstration, return mock data
    metadata = {
        "accession": dataset.accession,
        "title": dataset.title,
        "organism": dataset.organism,
        "platform": dataset.platform,
        "samples": [
            {
                "sample_id": f"{dataset.accession}_sample{i}",
                "condition": "treatment" if i % 2 == 0 else "control",
                "replicate": (i // 2) + 1,
            }
            for i in range(1, dataset.samples + 1)
        ]
    }
    return metadata

@DataExtractionAgent.tool
def download_files(ctx: RunContext, files: List[DataFile], download_dir: str) -> List[DataFile]:
    """
    Download the specified files.

    Args:
        files: List of DataFile objects.
        download_dir: Directory to download files to.

    Returns:
        Updated list of DataFile objects with download status and file paths.
    """
    # In a real implementation, this would download files from URLs
    # For demonstration, update the objects as if they were downloaded
    os.makedirs(download_dir, exist_ok=True)

    for file in files:
        file_name = file.file_url.split("/")[-1]
        file.file_path = os.path.join(download_dir, file_name)
        file.downloaded = True

    return files

# %% Data Analysis Agent
DataAnalysisAgent = Agent[WorkflowDependencies, AnalysisResult](
    'openai:gpt-4o',  # Specify a model
    deps_type=WorkflowDependencies,  # Specify the dependency type
    result_type=AnalysisResult,  # Specify the expected return type
    system_prompt="You are a specialist in analyzing RNA-seq data to identify patterns, differential expression, and biological pathways. "
        "Perform quantification, differential expression analysis, and pathway analysis to generate comprehensive insights."
)

@DataAnalysisAgent.tool
def perform_quantification(ctx: RunContext, fastq_files: List[DataFile], metadata: Dict[str, Any]) -> Dict[str, Any]:
    """
    Perform quantification on RNA-seq data.

    Args:
        fastq_files: List of FASTQ DataFile objects.
        metadata: Metadata dictionary.

    Returns:
        A dictionary containing quantification results.
    """
    # In a real implementation, this would use tools like Salmon, STAR, etc.
    # For demonstration, return mock data
    counts_matrix = {
        "gene1": [100, 120, 90, 110, 95, 105],
        "gene2": [200, 220, 190, 210, 195, 205],
        "gene3": [300, 320, 290, 310, 295, 305],
        # ... more genes
    }

    return {
        "counts_matrix": counts_matrix,
        "samples": [f.file_id.split('_R')[0] for f in fastq_files if '_R1' in f.file_id]
    }

@DataAnalysisAgent.tool
def generate_design_matrix(ctx: RunContext, metadata: Dict[str, Any]) -> pd.DataFrame:
    """
    Generate a design matrix for differential expression analysis.

    Args:
        metadata: Metadata dictionary.

    Returns:
        A pandas DataFrame representing the design matrix.
    """
    # Extract sample information from metadata
    samples = metadata['samples']

    # Create a design matrix
    design_matrix = pd.DataFrame([
        {
            'sample_id': sample['sample_id'],
            'condition': sample['condition'],
            'replicate': sample['replicate']
        }
        for sample in samples
    ])

    return design_matrix

@DataAnalysisAgent.tool
def perform_differential_expression(ctx: RunContext, counts: Dict[str, Any], design: pd.DataFrame) -> Dict[str, Any]:
    """
    Perform differential expression analysis.

    Args:
        counts: Counts dictionary from quantification.
        design: Design matrix DataFrame.

    Returns:
        A dictionary containing differential expression results.
    """
    # In a real implementation, this would use tools like DESeq2, edgeR, etc.
    # For demonstration, return mock data
    de_results = {
        "gene1": {"log2FoldChange": 1.5, "padj": 0.001},
        "gene2": {"log2FoldChange": -0.5, "padj": 0.05},
        "gene3": {"log2FoldChange": 2.0, "padj": 0.0001},
        # ... more genes
    }

    return {
        "de_results": de_results,
        "significant_genes": ["gene1", "gene3"]
    }

@DataAnalysisAgent.tool
def perform_pathway_analysis(ctx: RunContext, de_results: Dict[str, Any]) -> List[str]:
    """
    Perform pathway analysis on differentially expressed genes.

    Args:
        de_results: Differential expression results.

    Returns:
        A list of enriched pathways.
    """
    # In a real implementation, this would use tools like GSEA, Enrichr, etc.
    # For demonstration, return mock data
    significant_genes = de_results.get("significant_genes", [])

    if not significant_genes:
        return []

    # Mock pathway results
    pathways = [
        "Cell Cycle",
        "DNA Repair",
        "Immune Response",
        "Metabolism"
    ]

    return pathways

@DataAnalysisAgent.tool
def generate_results_summary(
    ctx: RunContext,
    dataset: GEODataset,
    de_results: Dict[str, Any],
    pathways: List[str]
) -> AnalysisResult:
    """
    Generate a summary of the analysis results.

    Args:
        dataset: GEODataset object.
        de_results: Differential expression results.
        pathways: List of enriched pathways.

    Returns:
        An AnalysisResult object.
    """
    significant_genes = de_results.get("significant_genes", [])

    summary = f"""
    Analysis of dataset {dataset.accession}: {dataset.title}

    Found {len(significant_genes)} differentially expressed genes.
    Top pathways: {', '.join(pathways[:3])}.

    The results suggest involvement of these pathways in the studied condition.
    """

    visualizations = [
        "volcano_plot.png",
        "heatmap.png",
        "pca_plot.png"
    ]

    return AnalysisResult(
        dataset_accession=dataset.accession,
        differentially_expressed_genes=significant_genes,
        enriched_pathways=pathways,
        visualizations=visualizations,
        summary=summary.strip()
    )

# %% Master Agent
RNAseqResearchAgent = Agent(
    'openai:gpt-4o',  # Specify a model
    system_prompt=(
        "You are a master agent coordinating an RNA-seq research workflow. "
        "Use tools to identify datasets, extract data, and perform analysis. "
        "Examine the user's query, identify relevant RNA-seq datasets, extract data files, "
        "and perform a complete analysis to generate meaningful insights."
    ),
    deps_type=WorkflowDependencies,  # Specify the dependency type
    result_type=WorkflowState  # Specify structured result type
)

@RNAseqResearchAgent.tool
async def process_dataset_identification(ctx: RunContext, query: str) -> List[GEODataset]:
    """
    Execute the dataset identification step of the workflow by delegating to the DatasetIdentificationAgent.

    Args:
        query: The research query.

    Returns:
        A list of relevant GEO datasets.
    """
    try:
        # In a real implementation with PydanticAI:
        result = await DatasetIdentificationAgent.run(
            f"Find relevant datasets for the query: {query}. " +
            "Evaluate their relevance and return a list of the most relevant datasets.",
            usage=ctx.usage
        )
        return result.data
    except Exception as e:
        # For the mock implementation, fall back to direct function calls:
        print(f"Note: In mock mode. Would delegate to DatasetIdentificationAgent. Error: {str(e)}")
        mock_ctx = RunContext(None)
        datasets = await query_geo_datasets(mock_ctx, query)
        relevant_datasets = await evaluate_dataset_relevance(mock_ctx, datasets, query)
        return relevant_datasets

@RNAseqResearchAgent.tool
async def select_best_dataset_for_query(ctx: RunContext, datasets: List[GEODataset]) -> GEODataset:
    """
    Select the best dataset from the list based on relevance by delegating to the DatasetIdentificationAgent.

    Args:
        datasets: List of GEODataset objects with relevance scores.

    Returns:
        The best GEODataset object.
    """
    try:
        # In a real implementation with PydanticAI:
        result = await DatasetIdentificationAgent.run(
            "Select the best dataset from these options based on relevance score.",
            usage=ctx.usage
        )
        return result.data
    except Exception as e:
        # For the mock implementation, fall back to direct function calls:
        print(f"Note: In mock mode. Would delegate to DatasetIdentificationAgent. Error: {str(e)}")
        mock_ctx = RunContext(None)
        return await select_best_dataset(mock_ctx, datasets)

@RNAseqResearchAgent.tool
async def extract_dataset_files(ctx: RunContext, dataset: GEODataset) -> Dict[str, Any]:
    """
    Extract files and metadata from the selected dataset by delegating to the DataExtractionAgent.

    Args:
        dataset: The selected GEO dataset.

    Returns:
        A dictionary containing fastq files and metadata.
    """
    try:
        # In a real implementation with PydanticAI:
        result = await DataExtractionAgent.run(
            f"Extract all files and metadata from dataset {dataset.accession}. " +
            f"The dataset is titled '{dataset.title}' and contains {dataset.samples} samples. " +
            "Download all FASTQ files to the temp directory.",
            usage=ctx.usage
        )
        return result.data
    except Exception as e:
        # For the mock implementation, fall back to direct function calls:
        print(f"Note: In mock mode. Would delegate to DataExtractionAgent. Error: {str(e)}")
        mock_ctx = RunContext(None)
        fastq_files = await extract_fastq_files(mock_ctx, dataset)
        metadata = await extract_metadata(mock_ctx, dataset)

        # Download the files
        download_dir = "/tmp/rnaseq_data"
        updated_files = await download_files(mock_ctx, fastq_files, download_dir)

        return {
            "fastq_files": updated_files,
            "metadata": metadata
        }

@RNAseqResearchAgent.tool
async def analyze_dataset(
    ctx: RunContext,
    dataset: GEODataset,
    files: List[DataFile],
    metadata: Dict[str, Any]
) -> AnalysisResult:
    """
    Perform analysis on the dataset files by delegating to the DataAnalysisAgent.

    Args:
        dataset: The selected GEO dataset.
        files: The extracted data files.
        metadata: The dataset metadata.

    Returns:
        Analysis results.
    """
    try:
        # In a real implementation with PydanticAI:
        # Prepare a detailed prompt for the analysis agent
        analysis_prompt = (
            f"Analyze dataset {dataset.accession} ({dataset.title}). "
            f"Perform RNA-seq analysis including quantification, differential expression, "
            f"and pathway analysis. Generate a summary of the findings."
        )

        result = await DataAnalysisAgent.run(
            analysis_prompt,
            usage=ctx.usage
        )
        return result.data
    except Exception as e:
        # For the mock implementation, fall back to direct function calls:
        print(f"Note: In mock mode. Would delegate to DataAnalysisAgent. Error: {str(e)}")
        mock_ctx = RunContext(None)
        counts = await perform_quantification(mock_ctx, files, metadata)
        design_matrix = await generate_design_matrix(mock_ctx, metadata)
        de_results = await perform_differential_expression(mock_ctx, counts, design_matrix)
        pathways = await perform_pathway_analysis(mock_ctx, de_results)

        return await generate_results_summary(mock_ctx, dataset, de_results, pathways)

@RNAseqResearchAgent.tool
async def execute_workflow(ctx: RunContext[WorkflowDependencies], query: str) -> WorkflowState:
    """
    Execute the complete RNA-seq research workflow.

    Args:
        query: The research query.

    Returns:
        A WorkflowState object representing the final state of the workflow.
    """
    # Ensure the shared temp directory exists
    ctx.deps.ensure_temp_dir()

    # Initialize the workflow state with the query
    workflow_state = WorkflowState(query=query)

    try:
        # Step 1: Dataset Identification
        workflow_state.current_step = WorkflowStep.DATASET_IDENTIFICATION
        print(f"Starting dataset identification for query: {query}")
        relevant_datasets = await process_dataset_identification(ctx, query)
        workflow_state.datasets = relevant_datasets

        if not relevant_datasets:
            workflow_state.status = AnalysisStatus.FAILURE
            workflow_state.error_message = "No relevant datasets found"
            return workflow_state

        selected_dataset = await select_best_dataset_for_query(ctx, relevant_datasets)
        workflow_state.selected_dataset = selected_dataset
        print(f"Selected dataset: {selected_dataset.accession} - {selected_dataset.title}")

        # Step 2: Data Extraction
        workflow_state.current_step = WorkflowStep.DATA_EXTRACTION
        print(f"Extracting data files from dataset {selected_dataset.accession}")
        extraction_result = await extract_dataset_files(ctx, selected_dataset)
        workflow_state.data_files = extraction_result["fastq_files"]
        print(f"Extracted {len(workflow_state.data_files)} data files")

        # Step 3: Data Analysis
        workflow_state.current_step = WorkflowStep.DATA_ANALYSIS
        print(f"Starting analysis of dataset {selected_dataset.accession}")
        analysis_results = await analyze_dataset(
            ctx,
            selected_dataset,
            extraction_result["fastq_files"],
            extraction_result["metadata"]
        )
        workflow_state.analysis_results = analysis_results
        print(f"Analysis completed successfully")

        # Complete the workflow
        workflow_state.current_step = WorkflowStep.COMPLETED
        workflow_state.status = AnalysisStatus.SUCCESS
        print(f"Workflow execution completed successfully")

    except Exception as e:
        workflow_state.status = AnalysisStatus.FAILURE
        workflow_state.error_message = str(e)
        print(f"Workflow execution failed: {str(e)}")

    return workflow_state

# %% Testing the Workflow
async def test_workflow():
    query = "Analyze RNAseq datasets for breast cancer"

    print(f"Executing workflow for query: {query}")

    try:
        # This is how you would use it in a real implementation with PydanticAI:
        # Create workflow dependencies to be shared across agents
        workflow_deps = WorkflowDependencies(query=query)

        # Run the master agent with proper query
        run_result = await RNAseqResearchAgent.run(
            f"Execute RNA-seq workflow for the query: {query}. " +
            "Identify relevant datasets, extract data files, and perform analysis.",
            # In real implementation, you'd pass the usage object to track token usage
            # usage=usage_tracker,
            deps=workflow_deps
        )

        # Get the structured result
        result = run_result.data
        print(f"Workflow completed with status: {result.status}")
    except Exception as e:
        print(f"Unable to run with real PydanticAI agent, using mock implementation: {str(e)}")

        # For the mock implementation, fall back to direct function calls:
        mock_ctx = RunContext(None)
        result = await RNAseqResearchAgent.execute_workflow(mock_ctx, query)
        print(f"Workflow completed with status: {result.status}")

    if result.status == AnalysisStatus.SUCCESS:
        print(f"Selected dataset: {result.selected_dataset.accession} - {result.selected_dataset.title}")
        print(f"Number of data files: {len(result.data_files)}")
        print(f"Analysis results summary:\n{result.analysis_results.summary}")
    else:
        print(f"Workflow failed with error: {result.error_message}")

    return result

def test_workflow_sync():
    """Synchronous wrapper for test_workflow"""
    import asyncio
    return asyncio.run(test_workflow())

# %% Run the test if this script is executed directly
if __name__ == "__main__":
    test_result = test_workflow_sync()

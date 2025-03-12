# %% Imports
import os
import sys
import time
import json
import csv
import re
import statistics
import subprocess
import tempfile
from io import StringIO
from enum import Enum
from typing import List, Dict, Any, Optional, Union, Literal

from functools import wraps  # For decorator
import logging
import pandas as pd
import numpy as np
from pydantic import BaseModel, Field

# Third‐party modules
from Bio import Entrez
import requests

# Configure logging
from logging.handlers import RotatingFileHandler

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)  # Overall logging level

# Create handlers
console_handler = logging.StreamHandler(sys.stdout)
file_handler = RotatingFileHandler("workflow.log", maxBytes=5*1024*1024, backupCount=5)

# Create formatters and add them to handlers
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
console_handler.setFormatter(formatter)
file_handler.setFormatter(formatter)

# Add handlers to the logger
logger.addHandler(console_handler)
logger.addHandler(file_handler)

# Suppress INFO logs from httpx to reduce noise
logging.getLogger("httpx").setLevel(logging.WARNING)
from dotenv import load_dotenv
load_dotenv()

# Prepare Entrez settings
Entrez.email = os.getenv("ENTREZ_EMAIL")
Entrez.api_key = os.getenv("ENTREZ_API_KEY")
# For running async code in Jupyter notebooks
try:
    import nest_asyncio
    nest_asyncio.apply()
except ImportError:
    logger.warning("nest_asyncio not found. Install it for better notebook compatibility: pip install nest_asyncio")

# Add these if you have PydanticAI installed
try:
    from pydantic_ai import Agent, RunContext
    from pydantic_ai.tools import Tool
    from pydantic_ai.usage import Usage
except ImportError:
    logger.warning("PydanticAI not installed. Using mock classes for demonstration.")
    class Agent:
        def __init__(self, name=None, description=None, deps_type=None, system_prompt=None, result_type=None):
            self.name = name
            self.description = description
            self.deps_type = deps_type
            self.system_prompt = system_prompt
            self.result_type = result_type
        def tool(self, func=None):
            return func
        def tool_plain(self, func=None):
            return func
        async def run(self, prompt, usage=None, deps=None):
            class DummyResult:
                def __init__(self, data):
                    self.data = data
            return DummyResult(prompt)
    class RunContext:
        def __init__(self, deps=None, usage=None):
            self.deps = deps
            self.usage = usage
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

    class Config:
        arbitrary_types_allowed = True
    organism: str
    samples: int
    platform: str
    relevance_score: Optional[float] = None

    class Config:
        arbitrary_types_allowed = True

class DataFile(BaseModel):
    """A model representing a data file from a GEO dataset."""
    dataset_accession: str
    file_id: str
    file_type: str
    file_url: str
    file_size: Optional[int] = None
    downloaded: bool = False
    file_path: Optional[str] = None

    class Config:
        arbitrary_types_allowed = True

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

    class Config:
        arbitrary_types_allowed = True

class WorkflowDependencies(BaseModel):
    """Dependencies shared across the workflow agents."""
    query: str
    temp_dir: str = "/tmp/rnaseq_data"

    class Config:
        arbitrary_types_allowed = True

    def get_download_path(self, filename: str) -> str:
        return os.path.join(self.temp_dir, filename)

    def ensure_temp_dir(self) -> None:
        os.makedirs(self.temp_dir, exist_ok=True)

# %% Helper Functions for GEO search
def perform_search(term: str) -> Dict[str, Any]:
    """Perform a GEO search using Entrez esearch."""
    handle = Entrez.esearch(db="gds", term=term, retmode="xml", retmax=50)
    results = Entrez.read(handle)
    handle.close()
    return results

def create_geo_dataframe(geo_ids: List[str], batch_size: int = 10) -> pd.DataFrame:
    """Create a DataFrame from GEO search results using batch processing."""
    data = []
    from tqdm import tqdm
    for i in tqdm(range(0, len(geo_ids), batch_size), desc="Processing GEO IDs in batches"):
        batch_ids = geo_ids[i:i + batch_size]
        ids_str = ",".join(batch_ids)
        handle = Entrez.esummary(db="gds", id=ids_str, retmode="xml")
        output = Entrez.read(handle)
        handle.close()
        for geo_id, geo_data in zip(batch_ids, output):
            if isinstance(geo_data, dict):
                data.append({
                    'ID': geo_id,
                    'Title': geo_data.get('title', 'No title available'),
                    'Summary': geo_data.get('summary', 'No summary available'),
                    'Accession': geo_data.get('Accession', geo_id),
                    'Species': geo_data.get('taxon', 'No taxon available'),
                    'Date': geo_data.get('PDAT', 'Unknown')
                })
            else:
                data.append({'ID': geo_id, 'Title': 'Error', 'Summary': 'Unable to fetch data', 'Accession': 'Error'})
        time.sleep(0.2)
    return pd.DataFrame(data)

# %% Define the Logging Decorator
def log_function_call(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        func_name = func.__name__
        logger.info(f"{func_name}: Function called with args={args[1:]}, kwargs={kwargs}")
        try:
            result = func(*args, **kwargs)
            logger.info(f"{func_name}: Function executed successfully.")
            logger.debug(f"{func_name}: Result: {result}")
            return result
        except Exception as e:
            logger.error(f"{func_name}: Error during execution: {e}", exc_info=True)
            raise
    return wrapper
DatasetIdentificationAgent = Agent[None, List[GEODataset]](
    'openai:gpt-4o',
    system_prompt=(
        "You are a specialist in identifying relevant genomic datasets from GEO databases for research queries. "
        "Use actual GEO search results to rank datasets by relevance."
    ),
    result_type=List[GEODataset]
)

@DatasetIdentificationAgent.tool
def query_geo_datasets(ctx: RunContext, query: str, max_results: int = 10) -> List[GEODataset]:
    """Query NCBI GEO for datasets matching the given query using Entrez."""
    search_results = perform_search(query)
    geo_ids = search_results.get("IdList", [])[:max_results]
    if not geo_ids:
        return []
    df = create_geo_dataframe(geo_ids)
    datasets = []
    for _, row in df.iterrows():
        ds = GEODataset(
            accession=row["Accession"],
            title=row["Title"],
            summary=row["Summary"],
            organism=row["Species"],
            samples=int(row.get("samples", 6)),  # Use 6 if missing
            platform="Illumina HiSeq 2500"  # Hardcoded; replace if available
        )
        datasets.append(ds)
    return datasets

@DatasetIdentificationAgent.tool
def evaluate_dataset_relevance(ctx: RunContext, datasets: List[GEODataset], query: str) -> List[GEODataset]:
    """Evaluate relevance based on simple keyword matching in title and summary."""
    query_lower = query.lower()
    for ds in datasets:
        score = 0.0
        if query_lower in ds.title.lower():
            score += 1.0
        if query_lower in ds.summary.lower():
            score += 0.5
        ds.relevance_score = score
    return sorted(datasets, key=lambda x: x.relevance_score, reverse=True) # [CHANGE] I will prefer to use an LLM prompt

@DatasetIdentificationAgent.tool
def select_best_dataset(ctx: RunContext, datasets: List[GEODataset]) -> GEODataset:
    if not datasets:
        raise ValueError("No datasets provided")
    return max(datasets, key=lambda x: x.relevance_score or 0) # [CHANGE] This seems to only return a single dataset, so will need to be changed

# %% Data Extraction Agent

# In general need to revise how this is done

DataExtractionAgent = Agent[WorkflowDependencies, Dict[str, Any]](
    'openai:gpt-4o',
    deps_type=WorkflowDependencies,
    result_type=Dict[str, Any],
    system_prompt=(
        "You are a specialist in extracting genomic data from GEO datasets. "
        "Retrieve FASTQ file URLs and metadata using standard GEO conventions."
    )
)

@DataExtractionAgent.tool
def extract_fastq_files(ctx: RunContext, dataset: GEODataset) -> List[DataFile]: # This caused problems for me in the past, so will need to look at this carefully
    """Construct FASTQ file URLs based on GEO accession."""
    fastq_files = []
    for i in range(1, dataset.samples + 1):
        for j in range(1, 3):  # Assume paired-end - [CHANGE] cannot assume this
            file_id = f"{dataset.accession}_sample{i}_R{j}"
            # Construct URL using a typical GEO FTP path (this is illustrative)
            file_url = f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{dataset.accession[:-3]}nnn/{dataset.accession}/{file_id}.fastq.gz"
            fastq_files.append(DataFile(
                dataset_accession=dataset.accession,
                file_id=file_id,
                file_type="fastq",
                file_url=file_url
            ))
    return fastq_files

@DataExtractionAgent.tool
@log_function_call
def extract_metadata(ctx: RunContext, dataset: GEODataset) -> Dict[str, Any]:
    """Extract metadata by using dataset information and standard GEO structure."""
    metadata = {
        "accession": dataset.accession,
        "title": dataset.title,
        "organism": dataset.organism,
        "platform": dataset.platform,
        "samples": [
            {
                "sample_id": f"{dataset.accession}_sample{i}",
                "condition": "treatment" if i % 2 == 0 else "control",
                "replicate": (i // 2) + 1
            }
            for i in range(1, dataset.samples + 1)
        ]
    }
    return metadata

@DataExtractionAgent.tool
@log_function_call
def download_files(ctx: RunContext, files: List[DataFile], download_dir: str) -> List[DataFile]:
    """Download files using HTTP requests."""
    os.makedirs(download_dir, exist_ok=True)
    for file in files:
        try:
            r = requests.get(file.file_url, stream=True)
            if r.status_code == 200:
                file_name = os.path.basename(file.file_url)
                file_path = os.path.join(download_dir, file_name)
                with open(file_path, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)
                file.downloaded = True
                file.file_path = file_path
            else:
                file.downloaded = False
        except Exception as e:
            file.downloaded = False
    return files

# %% Data Analysis Agent
try:
    DataAnalysisAgent = Agent[WorkflowDependencies, AnalysisResult](
        'openai:gpt-4o',
        deps_type=WorkflowDependencies,
        result_type=AnalysisResult,
        system_prompt=(
            "You are a bioinformatics expert analyzing RNA-seq data. "
            "Perform quantification, generate design matrices, run differential expression analysis, "
            "and carry out pathway enrichment to produce a summary of findings."
        )
    )
    logger.info("DataAnalysisAgent instantiated successfully.")
except TypeError as te:
    logger.error(f"Failed to instantiate DataAnalysisAgent: {te}", exc_info=True)
except Exception as e:
    logger.error(f"Unexpected error during DataAnalysisAgent instantiation: {e}", exc_info=True)

@DataAnalysisAgent.tool
@log_function_call
def perform_quantification(ctx: RunContext, fastq_files: List[DataFile], metadata: Dict[str, Any]) -> Dict[str, Any]:
    """
    In practice, you would call a quantification tool (e.g. kallisto, Salmon).
    Here we simulate quantification by invoking a kallisto-like command.
    """
    # For example, you might run:
    # subprocess.run(["kallisto", "quant", "-i", index, "-o", output, ...])
    # Here, we simulate by returning a dummy counts matrix.
    counts_matrix = {
        "gene1": [100, 110, 105],
        "gene2": [200, 210, 205],
        "gene3": [300, 310, 305],
    }
    return {"counts_matrix": counts_matrix, "samples": [f.file_id.split('_R')[0] for f in fastq_files if "_R1" in f.file_id]}

@DataAnalysisAgent.tool
def generate_design_matrix(ctx: RunContext, metadata: Dict[str, Any]) -> Dict[str, List[Dict[str, Any]]]:
    """Generate a design matrix from sample metadata."""
    samples = metadata.get("samples", [])
    # Return the samples directly as a dictionary instead of a DataFrame
    return {"samples": samples}

@DataAnalysisAgent.tool
def perform_differential_expression(ctx: RunContext, counts: Dict[str, Any], design: Dict[str, List[Dict[str, Any]]]) -> Dict[str, Any]:
    """
    Simulate differential expression analysis.
    In a production system, you would use statistical tools (e.g., edgeR, DESeq2).
    """
    # Convert design dict to DataFrame internally if needed for computation
    # samples = design.get("samples", [])
    # design_df = pd.DataFrame(samples)

    de_results = {
        "gene1": {"log2FoldChange": 1.2, "padj": 0.01},
        "gene2": {"log2FoldChange": -0.8, "padj": 0.05},
        "gene3": {"log2FoldChange": 1.5, "padj": 0.001},
    }
    return {"de_results": de_results, "significant_genes": ["gene1", "gene3"]}

@DataAnalysisAgent.tool
def perform_pathway_analysis(ctx: RunContext, de_results: Dict[str, Any]) -> List[str]:
    """Simulate pathway enrichment analysis."""
    return ["Cell Cycle", "DNA Repair", "Immune Response"]

# Helper function in case we need to create a dummy dataset during errors
def create_mock_dataset() -> GEODataset:
    """Create a mock GEODataset for testing purposes when the original flow fails."""
    return GEODataset(
        accession="GSE12345",
        title="Mock RNA-seq dataset for testing",
        summary="This is a placeholder dataset for testing purposes only.",
        organism="Homo sapiens",
        samples=6,
        platform="Illumina HiSeq 2500",
        relevance_score=1.0
    )
@DataAnalysisAgent.tool
def generate_results_summary(ctx: RunContext, dataset: GEODataset, de_results: Dict[str, Any], pathways: List[str]) -> AnalysisResult:
    """Compile analysis results into a summary."""
    significant_genes = de_results.get("significant_genes", [])
    summary = (
        f"Dataset {dataset.accession} ({dataset.title}) analysis:\n"
        f"Identified {len(significant_genes)} significant genes. Top pathways include: {', '.join(pathways[:3])}."
    )
    visualizations = ["volcano_plot.png", "heatmap.png", "pca_plot.png"]
    return AnalysisResult(
        dataset_accession=dataset.accession,
        differentially_expressed_genes=significant_genes,
        enriched_pathways=pathways,
        visualizations=visualizations,
        summary=summary
    )

# %% Master Agent: RNAseqResearchAgent
RNAseqResearchAgent = Agent[WorkflowDependencies, WorkflowState](
    'openai:gpt-4o',
    system_prompt=(
        "You are a master agent coordinating an RNA-seq research workflow. "
        "Delegate tasks to specialized agents to (1) identify relevant GEO datasets, "
        "(2) extract data files and metadata, and (3) analyze the data to produce meaningful insights. "
        "Maintain the original research query context throughout."
    ),
    deps_type=WorkflowDependencies,
    result_type=WorkflowState
)

@RNAseqResearchAgent.tool
async def run_agent_with_logging(agent: Agent, prompt: str, usage: Usage, deps: Optional[WorkflowDependencies] = None) -> Any:
    """
    Helper function to run an agent with logging of prompts and responses.
    
    Args:
        agent (Agent): The agent to run.
        prompt (str): The prompt to send to the agent.
        usage (Usage): Usage tracker.
        deps (Optional[WorkflowDependencies]): Dependencies for the agent.
    
    Returns:
        Any: The result returned by the agent.
    """
    logger.info(f"Sending prompt to {agent.name}: {prompt}")
    try:
        response = await agent.run(prompt, usage=usage, deps=deps)
        logger.info(f"Received response from {agent.name}: {response.data}")
        return response
    except Exception as e:
        logger.error(f"Error running agent {agent.name} with prompt '{prompt}': {e}", exc_info=True)
        raise
    logger.info("Starting workflow execution.")
    ctx.deps.ensure_temp_dir()
    workflow_state = WorkflowState(query=query)
    try:
        # Step 1: Dataset Identification
        workflow_state.current_step = WorkflowStep.DATASET_IDENTIFICATION
        logger.info(f"Dataset identification for query: {query}")
        ds_result = await run_agent_with_logging(
            f"Identify GEO datasets for: {query}",
            usage=ctx.usage
        )
        relevant_datasets = ds_result.data
        workflow_state.datasets = relevant_datasets
        if not relevant_datasets:
            workflow_state.status = AnalysisStatus.FAILURE
            workflow_state.error_message = "No relevant datasets found"
            return workflow_state

        # Select best dataset using sub-agent delegation
        import json
        try:
            top_datasets_json = json.dumps([d.dict() for d in relevant_datasets[:5]])
            sel_result = await run_agent_with_logging(
                f"Select the best dataset from these (JSON): {top_datasets_json}",
                usage=ctx.usage
            )
            selected_dataset = sel_result.data

            # Safety check - if the result is None or empty, fall back to first dataset
            if not selected_dataset and relevant_datasets:
                logger.warning("Dataset selection returned empty result. Using first dataset.")
                selected_dataset = relevant_datasets[0]
            elif not selected_dataset:
                logger.warning("No datasets available. Creating mock dataset.")
                selected_dataset = create_mock_dataset()

            workflow_state.selected_dataset = selected_dataset
            logger.info(f"Selected dataset: {selected_dataset.accession}")
        except Exception as e:
            logger.error(f"Error selecting dataset: {str(e)}. Using first dataset or mock.")
            if relevant_datasets:
                selected_dataset = relevant_datasets[0]
            else:
                selected_dataset = create_mock_dataset()
            workflow_state.selected_dataset = selected_dataset
            logger.info(f"Using dataset: {selected_dataset.accession}")

        # Step 2: Data Extraction
        workflow_state.current_step = WorkflowStep.DATA_EXTRACTION
        logger.info(f"Extracting data from dataset {selected_dataset.accession}")
        extraction_deps = WorkflowDependencies(query=ctx.deps.query, temp_dir=ctx.deps.temp_dir)
        de_result = await run_agent_with_logging(
            f"Extract FASTQ files and metadata for dataset {selected_dataset.accession}.",
            usage=ctx.usage,
            deps=extraction_deps
        )
        extraction_data = de_result.data
        workflow_state.data_files = extraction_data["fastq_files"]
        logger.info(f"Extracted {len(workflow_state.data_files)} files.")

        # Step 3: Data Analysis
        workflow_state.current_step = WorkflowStep.DATA_ANALYSIS
        logger.info(f"Analyzing dataset {selected_dataset.accession}")
        analysis_deps = WorkflowDependencies(query=ctx.deps.query, temp_dir=ctx.deps.temp_dir)
        file_info = "\n".join([f"- {f.file_id} ({f.file_type})" for f in workflow_state.data_files[:5]])
        analysis_prompt = (
            f"Analyze dataset {selected_dataset.accession} ({selected_dataset.title}). "
            f"It contains {len(workflow_state.data_files)} data files, e.g.:\n{file_info}\n"
            "Perform quantification, differential expression, and pathway analysis, then summarize the findings."
        )
        an_result = await run_agent_with_logging(
            analysis_prompt,
            usage=ctx.usage,
            deps=analysis_deps
        )
        workflow_state.analysis_results = an_result.data
        logger.info("Analysis completed.")

        workflow_state.current_step = WorkflowStep.COMPLETED
        workflow_state.status = AnalysisStatus.SUCCESS
        logger.info("Workflow executed successfully.")
    except Exception as e:
        workflow_state.status = AnalysisStatus.FAILURE
        workflow_state.error_message = str(e)
        logger.error(f"Workflow failed: {str(e)}")
    logger.info("Workflow execution completed.")
    return workflow_state

# %% Testing the Workflow
async def test_workflow():
    query = "Analyze RNAseq datasets for breast cancer"
    logger.info(f"Executing workflow for query: {query}")
    from pydantic_ai.usage import Usage
    workflow_deps = WorkflowDependencies(query=query)
    usage_tracker = Usage()  # Track token usage across agents
    run_result = await RNAseqResearchAgent.run(
        f"Run RNA-seq workflow for: {query}",
        usage=usage_tracker,
        deps=workflow_deps
    )
    result = run_result.data
    logger.info(f"Workflow status: {result.status}")
    if result.status == AnalysisStatus.SUCCESS:
        logger.info(f"Dataset: {result.selected_dataset.accession} - {result.selected_dataset.title}")
        logger.info(f"Files extracted: {len(result.data_files)}")
        logger.info("Analysis summary:")
        logger.info(result.analysis_results.summary)
    else:
        logger.error(f"Error: {result.error_message}")
    logger.info(f"Token usage: {usage_tracker}")
    return result

def test_workflow_sync():
    """Synchronous wrapper for test_workflow that works in Jupyter notebooks"""
    import asyncio
    import nest_asyncio

    # Apply nest_asyncio to allow asyncio.run in a notebook that already has an event loop
    try:
        nest_asyncio.apply()
        return asyncio.run(test_workflow())
    except ImportError:
        logger.warning("Could not import nest_asyncio. If running in a notebook, please install: !pip install nest_asyncio")
        # Fallback to running directly if we're not in a notebook or if nest_asyncio is not available
        import sys
        if 'ipykernel' not in sys.modules:
            return asyncio.run(test_workflow())
        else:
            # Create a new event loop for Jupyter notebook use
            loop = asyncio.new_event_loop()
            asyncio.set_event_loop(loop)
            try:
                return loop.run_until_complete(test_workflow())
            finally:
                loop.close()

# For direct script execution
if __name__ == "__main__":
    test_result = test_workflow_sync()

from mcp.server.fastmcp import FastMCP
from pydantic_ai import Agent
import shutil, os, json
from pathlib import Path
import sys

# Add access to ResultsIntegrator and contrast relevance
sys.path.insert(0, str(Path(__file__).parent.parent))
from ResultsIntegration import ResultsIntegrator
from contrast_relevance import run_contrast_relevance

server = FastMCP("Utility-Tools")

# Keep the existing disk_usage tool
@server.tool()
async def disk_usage(path: str) -> str:
    """Return `du -sh` for a folder."""
    out = shutil.disk_usage(path)
    return json.dumps({"total": out.total, "used": out.used, "free": out.free})

# Add a tool to list datasets
@server.tool()
async def list_datasets(results_dir: str) -> str:
    """
    List available datasets in the UORCA results directory.

    Args:
        results_dir: Path to the UORCA results directory

    Returns:
        JSON string containing list of available dataset IDs and their basic info
    """
    try:
        # Initialize ResultsIntegrator
        ri = ResultsIntegrator(results_dir)
        ri.load_data()

        datasets = []
        for analysis_id, info in ri.analysis_info.items():
            datasets.append({
                "dataset_id": analysis_id,
                "accession": info.get("accession", "Unknown"),
                "organism": info.get("organism", "Unknown"),
                "num_samples": info.get("number_of_samples", 0),
                "num_contrasts": info.get("number_of_contrasts", 0)
            })

        return json.dumps({"datasets": datasets}, ensure_ascii=False)
    except Exception as e:
        return json.dumps({"error": str(e)}, ensure_ascii=False)

# Add the dataset summary tool
@server.tool()
async def get_dataset_summary(dataset_id: str, results_dir: str) -> str:
    """
    Get a comprehensive summary of a dataset including title, description, and basic statistics.

    Args:
        dataset_id: ID of the dataset to summarize
        results_dir: Path to the UORCA results directory

    Returns:
        JSON string containing dataset summary information
    """
    try:
        # Initialize ResultsIntegrator
        ri = ResultsIntegrator(results_dir)
        ri.load_data()

        # Check if the dataset exists
        if dataset_id not in ri.analysis_info:
            return json.dumps({
                "error": f"Dataset {dataset_id} not found",
                "available_datasets": list(ri.analysis_info.keys())
            })

        # Get basic dataset info
        analysis_info = ri.analysis_info.get(dataset_id, {})

        # Get dataset metadata
        dataset_info = ri.dataset_info.get(dataset_id, {}) if hasattr(ri, "dataset_info") else {}

        summary = {
            "dataset_id": dataset_id,
            "accession": analysis_info.get("accession", "Unknown"),
            "organism": analysis_info.get("organism", "Unknown"),
            "number_of_samples": analysis_info.get("number_of_samples", 0),
            "number_of_contrasts": analysis_info.get("number_of_contrasts", 0),
            "title": dataset_info.get("title", ""),
            "summary": dataset_info.get("summary", ""),
            "design": dataset_info.get("design", ""),
        }

        return json.dumps(summary, ensure_ascii=False)

    except Exception as e:
        return json.dumps({
            "error": str(e),
            "dataset_id": dataset_id
        }, ensure_ascii=False)

# Add contrast relevance assessment tool
@server.tool()
async def assess_contrast_relevance(results_dir: str, query: str) -> str:
    """
    Assess the relevance of all contrasts to a research question using AI.

    Args:
        results_dir: Path to the UORCA results directory
        query: Research question to assess contrast relevance against

    Returns:
        JSON string containing contrast relevance scores and justifications
    """
    try:
        # Initialize ResultsIntegrator
        ri = ResultsIntegrator(results_dir)
        ri.load_data()

        # Run contrast relevance assessment
        results_df = run_contrast_relevance(
            ri,
            query=query,
            repeats=2,
            batch_size=5,
            parallel_jobs=1
        )

        if results_df.empty:
            return json.dumps({
                "error": "No contrasts found for assessment",
                "contrasts": []
            }, ensure_ascii=False)

        # Add contrast descriptions for context
        results_df['description'] = results_df.apply(
            lambda row: ri._get_contrast_description(row['analysis_id'], row['contrast_id']),
            axis=1
        )

        # Add dataset info for context
        results_df['accession'] = results_df['analysis_id'].map(
            lambda aid: ri.analysis_info.get(aid, {}).get('accession', aid)
        )
        results_df['organism'] = results_df['analysis_id'].map(
            lambda aid: ri.analysis_info.get(aid, {}).get('organism', 'Unknown')
        )

        # Convert to JSON with all relevant information
        contrast_data = results_df.to_dict('records')

        return json.dumps({
            "query": query,
            "total_contrasts": len(contrast_data),
            "contrasts": contrast_data
        }, ensure_ascii=False)

    except Exception as e:
        return json.dumps({
            "error": str(e),
            "query": query
        }, ensure_ascii=False)

if __name__ == "__main__":
    server.run()

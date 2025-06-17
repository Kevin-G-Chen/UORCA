from mcp.server.fastmcp import FastMCP
from pydantic_ai import Agent
import shutil, os, json
from pathlib import Path
import sys

# Add access to ResultsIntegrator
sys.path.insert(0, str(Path(__file__).parent.parent))
from ResultsIntegration import ResultsIntegrator

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

if __name__ == "__main__":
    server.run()

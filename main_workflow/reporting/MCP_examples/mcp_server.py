from mcp.server.fastmcp import FastMCP
from pydantic_ai import Agent
import shutil, os, json
from pathlib import Path
import sys
import pandas as pd

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

# Add new tools for gene analysis
@server.tool()
async def filter_degs(
    results_dir: str,
    dataset_id: str,
    contrast_id: str,
    p_value_thresh: float,
    lfc_thresh: float
) -> str:
    """
    Filter genes from a contrast's DEG.csv file based on p-value and log fold change thresholds.

    Args:
        results_dir: Path to the UORCA results directory
        dataset_id: ID of the dataset/analysis
        contrast_id: ID of the contrast
        p_value_thresh: P-value threshold (genes with p < threshold are selected)
        lfc_thresh: Log fold change threshold (genes with |logFC| > threshold are selected)

    Returns:
        JSON string containing filtered gene list with statistics
    """
    try:
        # Find the DEG file
        deg_file = os.path.join(results_dir, dataset_id, "RNAseqAnalysis", contrast_id, "DEG.csv")

        if not os.path.exists(deg_file):
            return json.dumps({
                "error": f"DEG file not found: {deg_file}",
                "genes": []
            }, ensure_ascii=False)

        # Load DEG data
        deg_df = pd.read_csv(deg_file)

        # Ensure required columns exist
        if 'Gene' not in deg_df.columns:
            return json.dumps({
                "error": "Gene column not found in DEG file",
                "genes": []
            }, ensure_ascii=False)

        # Find p-value column (prefer adjusted)
        p_col = None
        if 'adj.P.Val' in deg_df.columns:
            p_col = 'adj.P.Val'
        elif 'P.Value' in deg_df.columns:
            p_col = 'P.Value'

        # Find log fold change column
        lfc_col = None
        if 'logFC' in deg_df.columns:
            lfc_col = 'logFC'
        elif 'log2FoldChange' in deg_df.columns:
            lfc_col = 'log2FoldChange'

        if not p_col or not lfc_col:
            return json.dumps({
                "error": f"Required columns not found. Available: {list(deg_df.columns)}",
                "genes": []
            }, ensure_ascii=False)

        # Apply filters
        filtered_df = deg_df[
            (deg_df[p_col] < p_value_thresh) &
            (abs(deg_df[lfc_col]) > lfc_thresh)
        ].copy()

        # Prepare output
        genes_list = []
        for _, row in filtered_df.iterrows():
            gene_info = {
                "gene": row['Gene'],
                "logFC": float(row[lfc_col]),
                "pvalue": float(row[p_col])
            }
            # Add adjusted p-value if available
            if 'adj.P.Val' in deg_df.columns and p_col != 'adj.P.Val':
                gene_info["adj_pvalue"] = float(row['adj.P.Val'])
            genes_list.append(gene_info)

        result = {
            "dataset_id": dataset_id,
            "contrast_id": contrast_id,
            "filters_applied": {
                "p_value_thresh": p_value_thresh,
                "lfc_thresh": lfc_thresh
            },
            "total_genes_in_file": len(deg_df),
            "filtered_genes_count": len(genes_list),
            "genes": genes_list
        }

        return json.dumps(result, ensure_ascii=False)

    except Exception as e:
        return json.dumps({
            "error": str(e),
            "dataset_id": dataset_id,
            "contrast_id": contrast_id
        }, ensure_ascii=False)

@server.tool()
async def get_gene_stats(
    results_dir: str,
    dataset_id: str,
    contrast_id: str,
    gene: str
) -> str:
    """
    Get statistics for a specific gene in a specific contrast's DEG results.

    Args:
        results_dir: Path to the UORCA results directory
        dataset_id: ID of the dataset/analysis
        contrast_id: ID of the contrast
        gene: Gene symbol to look up

    Returns:
        JSON string containing gene statistics
    """
    try:
        # Find the DEG file
        deg_file = os.path.join(results_dir, dataset_id, "RNAseqAnalysis", contrast_id, "DEG.csv")

        if not os.path.exists(deg_file):
            return json.dumps({
                "error": f"DEG file not found: {deg_file}",
                "found": False
            }, ensure_ascii=False)

        # Load DEG data
        deg_df = pd.read_csv(deg_file)

        # Find the gene
        if 'Gene' not in deg_df.columns:
            return json.dumps({
                "error": "Gene column not found in DEG file",
                "found": False
            }, ensure_ascii=False)

        gene_row = deg_df[deg_df['Gene'] == gene]

        if gene_row.empty:
            return json.dumps({
                "dataset_id": dataset_id,
                "contrast_id": contrast_id,
                "gene": gene,
                "found": False,
                "message": f"Gene {gene} not found in this contrast"
            }, ensure_ascii=False)

        # Extract statistics
        row = gene_row.iloc[0]
        stats = {
            "dataset_id": dataset_id,
            "contrast_id": contrast_id,
            "gene": gene,
            "found": True
        }

        # Add available statistics
        for col in ['logFC', 'log2FoldChange', 'AveExpr', 'P.Value', 'adj.P.Val', 't', 'B']:
            if col in deg_df.columns:
                stats[col] = float(row[col]) if pd.notna(row[col]) else None

        return json.dumps(stats, ensure_ascii=False)

    except Exception as e:
        return json.dumps({
            "error": str(e),
            "dataset_id": dataset_id,
            "contrast_id": contrast_id,
            "gene": gene,
            "found": False
        }, ensure_ascii=False)

@server.tool()
async def count_gene_occurrences(
    results_dir: str,
    selections: str,
    p_value_thresh: float = 0.05,
    lfc_thresh: float = 1.0
) -> str:
    """
    Count how many times each gene appears as differentially expressed across selected contrasts.

    Args:
        results_dir: Path to the UORCA results directory
        selections: JSON string of list of dicts with 'dataset_id' and 'contrast_id' keys
        p_value_thresh: P-value threshold for significance
        lfc_thresh: Log fold change threshold for significance

    Returns:
        JSON string containing gene occurrence counts and statistics
    """
    try:
        # Parse selections
        import json as json_lib
        selection_list = json_lib.loads(selections)

        gene_counts = {}
        contrast_info = []
        total_contrasts = len(selection_list)

        for selection in selection_list:
            dataset_id = selection['dataset_id']
            contrast_id = selection['contrast_id']

            # Find the DEG file
            deg_file = os.path.join(results_dir, dataset_id, "RNAseqAnalysis", contrast_id, "DEG.csv")

            if not os.path.exists(deg_file):
                contrast_info.append({
                    "dataset_id": dataset_id,
                    "contrast_id": contrast_id,
                    "status": "file_not_found",
                    "deg_count": 0
                })
                continue

            # Load and filter DEG data
            deg_df = pd.read_csv(deg_file)

            if 'Gene' not in deg_df.columns:
                contrast_info.append({
                    "dataset_id": dataset_id,
                    "contrast_id": contrast_id,
                    "status": "no_gene_column",
                    "deg_count": 0
                })
                continue

            # Find p-value and log fold change columns
            p_col = 'adj.P.Val' if 'adj.P.Val' in deg_df.columns else 'P.Value'
            lfc_col = 'logFC' if 'logFC' in deg_df.columns else 'log2FoldChange'

            if p_col not in deg_df.columns or lfc_col not in deg_df.columns:
                contrast_info.append({
                    "dataset_id": dataset_id,
                    "contrast_id": contrast_id,
                    "status": "missing_required_columns",
                    "deg_count": 0
                })
                continue

            # Filter significant genes
            sig_genes = deg_df[
                (deg_df[p_col] < p_value_thresh) &
                (abs(deg_df[lfc_col]) > lfc_thresh)
            ]['Gene'].tolist()

            # Count occurrences
            for gene in sig_genes:
                gene_counts[gene] = gene_counts.get(gene, 0) + 1

            contrast_info.append({
                "dataset_id": dataset_id,
                "contrast_id": contrast_id,
                "status": "success",
                "deg_count": len(sig_genes)
            })

        # Sort genes by occurrence count
        sorted_genes = sorted(gene_counts.items(), key=lambda x: x[1], reverse=True)

        result = {
            "total_contrasts_analyzed": total_contrasts,
            "filters_applied": {
                "p_value_thresh": p_value_thresh,
                "lfc_thresh": lfc_thresh
            },
            "unique_genes_found": len(gene_counts),
            "gene_counts": dict(sorted_genes),
            "top_recurring_genes": sorted_genes[:20],  # Top 20 most frequent
            "contrast_details": contrast_info
        }

        return json.dumps(result, ensure_ascii=False)

    except Exception as e:
        return json.dumps({
            "error": str(e),
            "selections": selections
        }, ensure_ascii=False)

if __name__ == "__main__":
    server.run()

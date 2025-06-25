import os
import json
import pandas as pd
from datetime import datetime
from mcp.server.fastmcp import FastMCP
from ResultsIntegration import ResultsIntegrator
from streamlit_tabs.helpers.ai_agent_tool_logger import log_ai_agent_tool


server = FastMCP("UORCA-Tools")

def load_long_table(results_dir: str) -> pd.DataFrame:
    """
    Build a DataFrame with one row per (analysis_id, contrast_id, gene)
    with columns: analysis_id, contrast_id, Gene, logFC, pvalue.
    """
    ri = ResultsIntegrator(results_dir)
    ri.load_data()
    pieces = []
    for aid, contrasts in ri.deg_data.items():
        for cid, df in contrasts.items():
            # normalize p-value column - prefer adjusted p-value if available
            df2 = df.copy()
            if 'adj.P.Val' in df.columns:
                df2 = df2.rename(columns={'adj.P.Val': 'pvalue'})
            elif 'P.Value' in df.columns:
                df2 = df2.rename(columns={'P.Value': 'pvalue'})

            # Ensure we have the required columns
            if 'Gene' not in df2.columns or 'logFC' not in df2.columns or 'pvalue' not in df2.columns:
                continue

            df_sub = df2[['Gene','logFC','pvalue']].copy()
            df_sub['analysis_id'] = aid
            df_sub['contrast_id'] = cid
            pieces.append(df_sub[['analysis_id','contrast_id','Gene','logFC','pvalue']])

    if not pieces:
        # Return empty DataFrame with correct structure if no data
        return pd.DataFrame(columns=['analysis_id','contrast_id','Gene','logFC','pvalue'])

    return pd.concat(pieces, ignore_index=True)

# Load upon import
RESULTS_DIR = os.getenv('RESULTS_DIR', '/UORCA_results')
try:
    FULL_DF = load_long_table(RESULTS_DIR)
    print(f"Loaded full dataframe: {len(FULL_DF)} rows")
except Exception as e:
    print(f"Warning: Could not load long table from {RESULTS_DIR}: {e}")
    # Create empty DataFrame with correct structure as fallback
    FULL_DF = pd.DataFrame(columns=['analysis_id','contrast_id','Gene','logFC','pvalue'])

def get_filtered_dataframe() -> pd.DataFrame:
    """
    Get the dataframe filtered to selected contrasts, checking environment variable dynamically.
    """
    global FULL_DF

    # Check for selected contrasts environment variable
    selected_contrasts_json = os.getenv('SELECTED_CONTRASTS_FOR_AI')
    if selected_contrasts_json:
        try:
            selected_contrasts = json.loads(selected_contrasts_json)
            print(f"Filtering to {len(selected_contrasts)} selected contrasts")

            # Create filter condition
            contrast_tuples = [(sc['analysis_id'], sc['contrast_id']) for sc in selected_contrasts]
            mask = FULL_DF.apply(lambda row: (row['analysis_id'], row['contrast_id']) in contrast_tuples, axis=1)
            filtered_df = FULL_DF[mask].copy()

            print(f"Filtered dataframe: {len(filtered_df)} rows from {len(FULL_DF)} total rows")
            return filtered_df
        except json.JSONDecodeError:
            print("Warning: Could not parse SELECTED_CONTRASTS_FOR_AI, using full dataframe")
            return FULL_DF
    else:
        print("No selected contrasts specified, using full dataframe")
        return FULL_DF

# 1) Top recurring DEGs
@server.tool()
@log_ai_agent_tool
async def get_most_common_genes(lfc_thresh: float, p_thresh: float, top_n: int) -> str:
    """
    Find the most commonly differentially expressed genes across contrasts. Useful to get an overview of frequently observed DEGs.

    Args:
        lfc_thresh: Minimum absolute log2 fold change threshold
        p_thresh: Maximum p-value threshold
        top_n: Number of top genes to return

    Returns:
        JSON string with list of genes and their occurrence counts
    """
    df = get_filtered_dataframe()
    df = df[(df.logFC.abs() >= lfc_thresh) & (df.pvalue < p_thresh)]
    counts = df.groupby("Gene").size().nlargest(top_n)
    result = json.dumps([{"gene": g, "count": int(c)} for g, c in counts.items()])

    return result

# 2) Per-gene + contrast stats
@server.tool()
@log_ai_agent_tool
async def get_gene_contrast_stats(gene: str, contrast_id: str = None) -> str:
    """
    Get statistics for a specific gene across contrasts. Useful to assess a specific gene to get more information about its expression changes across different conditions.

    Args:
        gene: Gene symbol to look up
        contrast_id: Optional specific contrast to filter to

    Returns:
        JSON string with gene statistics across contrasts
    """
    df = get_filtered_dataframe()
    df = df[df.Gene == gene]
    if contrast_id:
        df = df[df.contrast_id == contrast_id]
    result = json.dumps(df[["contrast_id","logFC","pvalue"]].to_dict("records"))

    return result

# 3) Filter by set A vs set B
@server.tool()
@log_ai_agent_tool
async def filter_genes_by_contrast_sets(set_a: list, set_b: list, lfc_thresh: float, p_thresh: float) -> str:
    """
    Find genes that are significant in contrast set A but not in set B. Useful for comparing two sets of contrasts, for example when wanting to find genes unique to a specific condition, and not (for example) a control condition.

    Args:
        set_a: List of contrast IDs for set A
        set_b: List of contrast IDs for set B
        lfc_thresh: Minimum absolute log2 fold change threshold
        p_thresh: Maximum p-value threshold

    Returns:
        JSON string with genes unique to set A and summary counts
    """
    df = get_filtered_dataframe()

    dfA = df[
        df.contrast_id.isin(set_a) &
        (df.logFC.abs() >= lfc_thresh) &
        (df.pvalue < p_thresh)
    ]
    genesA = set(dfA.Gene)

    dfB = df[
        df.contrast_id.isin(set_b) &
        (df.logFC.abs() >= lfc_thresh) &
        (df.pvalue < p_thresh)
    ]
    genesB = set(dfB.Gene)

    result_genes = sorted(genesA - genesB)
    result = json.dumps({
        "genes": result_genes,
        "counts": {
            "in_A": len(result_genes),
            "total_A": len(genesA),
            "total_B": len(genesB)
        }
    })

    return result

# 4) Contrast mini-summary
@server.tool()
@log_ai_agent_tool
async def summarize_contrast(contrast_id: str, lfc_thresh: float, p_thresh: float, max_genes: int = 10) -> str:
    """
    Summarize a specific contrast with top DEGs and statistics.

    Args:
        contrast_id: Contrast to summarize
        lfc_thresh: Minimum absolute log2 fold change threshold
        p_thresh: Maximum p-value threshold
        max_genes: Maximum number of top genes to include

    Returns:
        JSON string with contrast summary and top genes
    """
    df = get_filtered_dataframe()
    df = df[df.contrast_id == contrast_id]
    df = df[(df.logFC.abs() >= lfc_thresh) & (df.pvalue < p_thresh)]

    if df.empty:
        result = json.dumps({
            "total_DEGs": 0,
            "top_genes": [],
            "mean_logFC": None,
            "median_logFC": None
        })
    else:
        top = df.reindex(df.logFC.abs().sort_values(ascending=False).index).head(max_genes)
        result = json.dumps({
            "total_DEGs": int(len(df)),
            "top_genes": [{"gene": r.Gene, "logFC": float(r.logFC)} for _, r in top.iterrows()],
            "mean_logFC": float(df.logFC.mean()),
            "median_logFC": float(df.logFC.median())
        })

    return result



if __name__ == "__main__":
    server.run()

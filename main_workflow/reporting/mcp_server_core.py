import os
import json
import pandas as pd
from mcp.server.fastmcp import FastMCP
from ResultsIntegration import ResultsIntegrator
from streamlit_tabs.helpers import log_streamlit_function

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
    LONG_DF = load_long_table(RESULTS_DIR)
except Exception as e:
    print(f"Warning: Could not load long table from {RESULTS_DIR}: {e}")
    # Create empty DataFrame with correct structure as fallback
    LONG_DF = pd.DataFrame(columns=['analysis_id','contrast_id','Gene','logFC','pvalue'])

# 1) Top recurring DEGs
@server.tool()
@log_streamlit_function
async def get_most_common_genes(lfc_thresh: float, p_thresh: float, top_n: int) -> str:
    """
    Find the most commonly differentially expressed genes across contrasts.

    Args:
        lfc_thresh: Minimum absolute log2 fold change threshold
        p_thresh: Maximum p-value threshold
        top_n: Number of top genes to return

    Returns:
        JSON string with list of genes and their occurrence counts
    """
    df = LONG_DF[(LONG_DF.logFC.abs() >= lfc_thresh) & (LONG_DF.pvalue < p_thresh)]
    counts = df.groupby("Gene").size().nlargest(top_n)
    return json.dumps([{"gene": g, "count": int(c)} for g, c in counts.items()])

# 2) Per-gene + contrast stats
@server.tool()
@log_streamlit_function
async def get_gene_contrast_stats(gene: str, contrast_id: str = None) -> str:
    """
    Get statistics for a specific gene across contrasts.

    Args:
        gene: Gene symbol to look up
        contrast_id: Optional specific contrast to filter to

    Returns:
        JSON string with gene statistics across contrasts
    """
    df = LONG_DF[LONG_DF.Gene == gene]
    if contrast_id:
        df = df[df.contrast_id == contrast_id]
    return json.dumps(df[["contrast_id","logFC","pvalue"]].to_dict("records"))

# 3) Filter by set A vs set B
@server.tool()
@log_streamlit_function
async def filter_genes_by_contrast_sets(set_a: list, set_b: list, lfc_thresh: float, p_thresh: float) -> str:
    """
    Find genes that are significant in contrast set A but not in set B.

    Args:
        set_a: List of contrast IDs for set A
        set_b: List of contrast IDs for set B
        lfc_thresh: Minimum absolute log2 fold change threshold
        p_thresh: Maximum p-value threshold

    Returns:
        JSON string with genes unique to set A and summary counts
    """
    dfA = LONG_DF[
        LONG_DF.contrast_id.isin(set_a) &
        (LONG_DF.logFC.abs() >= lfc_thresh) &
        (LONG_DF.pvalue < p_thresh)
    ]
    genesA = set(dfA.Gene)

    dfB = LONG_DF[
        LONG_DF.contrast_id.isin(set_b) &
        (LONG_DF.logFC.abs() >= lfc_thresh) &
        (LONG_DF.pvalue < p_thresh)
    ]
    genesB = set(dfB.Gene)

    result = sorted(genesA - genesB)
    return json.dumps({
        "genes": result,
        "counts": {
            "in_A": len(result),
            "total_A": len(genesA),
            "total_B": len(genesB)
        }
    })

# 4) Contrast mini-summary
@server.tool()
@log_streamlit_function
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
    df = LONG_DF[LONG_DF.contrast_id == contrast_id]
    df = df[(df.logFC.abs() >= lfc_thresh) & (df.pvalue < p_thresh)]

    if df.empty:
        return json.dumps({
            "total_DEGs": 0,
            "top_genes": [],
            "mean_logFC": None,
            "median_logFC": None
        })

    top = df.reindex(df.logFC.abs().sort_values(ascending=False).index).head(max_genes)
    return json.dumps({
        "total_DEGs": int(len(df)),
        "top_genes": [{"gene": r.Gene, "logFC": float(r.logFC)} for _, r in top.iterrows()],
        "mean_logFC": float(df.logFC.mean()),
        "median_logFC": float(df.logFC.median())
    })

if __name__ == "__main__":
    server.run()

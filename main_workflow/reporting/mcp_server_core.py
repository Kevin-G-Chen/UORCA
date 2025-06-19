import os
import json
import pandas as pd
from datetime import datetime
from mcp.server.fastmcp import FastMCP
from ResultsIntegration import ResultsIntegrator
from streamlit_tabs.helpers import log_streamlit_function


server = FastMCP("UORCA-Tools")

# Global tool call log for tracking AI tool usage
tool_call_log = []

def log_tool_call(tool_name: str, params: dict, result: str):
    """Log tool calls for transparency in AI analysis."""
    tool_call_log.append({
        'tool': tool_name,
        'parameters': params,
        'result': result[:500] + "..." if len(result) > 500 else result,  # Truncate long results
        'timestamp': datetime.now().isoformat()
    })

def clear_tool_call_log():
    """Clear the tool call log."""
    global tool_call_log
    tool_call_log = []

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
    Find the most commonly differentially expressed genes across contrasts. Useful to get an overview of frequently observed DEGs.

    Args:
        lfc_thresh: Minimum absolute log2 fold change threshold
        p_thresh: Maximum p-value threshold
        top_n: Number of top genes to return

    Returns:
        JSON string with list of genes and their occurrence counts
    """
    params = {'lfc_thresh': lfc_thresh, 'p_thresh': p_thresh, 'top_n': top_n}

    df = LONG_DF[(LONG_DF.logFC.abs() >= lfc_thresh) & (LONG_DF.pvalue < p_thresh)]
    counts = df.groupby("Gene").size().nlargest(top_n)
    result = json.dumps([{"gene": g, "count": int(c)} for g, c in counts.items()])

    log_tool_call('get_most_common_genes', params, result)
    return result

# 2) Per-gene + contrast stats
@server.tool()
@log_streamlit_function
async def get_gene_contrast_stats(gene: str, contrast_id: str = None) -> str:
    """
    Get statistics for a specific gene across contrasts. Useful to assess a specific gene to get more information about its expression changes across different conditions.

    Args:
        gene: Gene symbol to look up
        contrast_id: Optional specific contrast to filter to

    Returns:
        JSON string with gene statistics across contrasts
    """
    params = {'gene': gene, 'contrast_id': contrast_id}

    df = LONG_DF[LONG_DF.Gene == gene]
    if contrast_id:
        df = df[df.contrast_id == contrast_id]
    result = json.dumps(df[["contrast_id","logFC","pvalue"]].to_dict("records"))

    log_tool_call('get_gene_contrast_stats', params, result)
    return result

# 3) Filter by set A vs set B
@server.tool()
@log_streamlit_function
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
    params = {'set_a': set_a, 'set_b': set_b, 'lfc_thresh': lfc_thresh, 'p_thresh': p_thresh}

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

    result_genes = sorted(genesA - genesB)
    result = json.dumps({
        "genes": result_genes,
        "counts": {
            "in_A": len(result_genes),
            "total_A": len(genesA),
            "total_B": len(genesB)
        }
    })

    log_tool_call('filter_genes_by_contrast_sets', params, result)
    return result

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
    params = {'contrast_id': contrast_id, 'lfc_thresh': lfc_thresh, 'p_thresh': p_thresh, 'max_genes': max_genes}

    df = LONG_DF[LONG_DF.contrast_id == contrast_id]
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

    log_tool_call('summarize_contrast', params, result)
    return result

# Tool to retrieve and manage tool call logs
@server.tool()
@log_streamlit_function
async def get_tool_call_log() -> str:
    """
    Get the current tool call log and clear it for next analysis.

    Returns:
        JSON string with list of tool calls made during analysis
    """
    global tool_call_log
    result = json.dumps(tool_call_log)
    return result

@server.tool()
@log_streamlit_function
async def clear_tool_log() -> str:
    """
    Clear the tool call log.

    Returns:
        Confirmation message
    """
    clear_tool_call_log()
    return json.dumps({"status": "cleared"})

if __name__ == "__main__":
    server.run()

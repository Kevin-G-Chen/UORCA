import os
from typing import Optional, Dict

import logging

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist

logger = logging.getLogger(__name__)


def load_deg(path: str) -> pd.DataFrame:
    """Load differential expression results from a CSV file."""
    return pd.read_csv(path)


def load_cpm(path: str) -> pd.DataFrame:
    """Load CPM expression values and set 'Gene' as index."""
    df = pd.read_csv(path)
    if "Gene" in df.columns:
        df = df.set_index("Gene")
    return df


def _first_existing(paths):
    for p in paths:
        if os.path.isfile(p):
            return p
    return None


def load_sample_groups(results_dir: str, analysis_id: str, cpm_df: pd.DataFrame, analysis_info: Optional[Dict] = None) -> Optional[pd.Series]:
    """Return a mapping of sample names to group labels if available."""
    analysis_dir = os.path.join(results_dir, analysis_id)
    sample_file_locations = [
        os.path.join(analysis_dir, "metadata", "edger_analysis_samples.csv"),
        os.path.join(analysis_dir, "edger_analysis_samples.csv"),
        os.path.join(results_dir, "edger_analysis_samples.csv"),
        os.path.join(os.path.dirname(results_dir), "edger_analysis_samples.csv"),
    ]
    sample_file = _first_existing(sample_file_locations)
    sample_cols = [c for c in cpm_df.columns if c != "Gene"]
    if not sample_file:
        return pd.Series({c: c for c in sample_cols})
    try:
        tmp = pd.read_csv(sample_file, nrows=0)
        first_col = tmp.columns[0]
        has_index = first_col == "" or first_col.startswith("Unnamed")
        sample_df = pd.read_csv(sample_file, index_col=0 if has_index else None)
        group_col = None
        if analysis_info and "analysis_column" in analysis_info:
            group_col = analysis_info["analysis_column"]
        elif "merged_analysis_group" in sample_df.columns:
            group_col = "merged_analysis_group"
        else:
            for col in sample_df.columns:
                if "group" in col.lower() or "condition" in col.lower():
                    group_col = col
                    break
        if group_col and group_col in sample_df.columns:
            sample_pattern = [c for c in sample_cols if c.startswith("Sample")]
            mapping = {}
            if sample_pattern:
                for i in range(min(len(sample_pattern), len(sample_df))):
                    mapping[sample_pattern[i]] = sample_df.iloc[i][group_col]
            else:
                for i in range(min(len(sample_cols), len(sample_df))):
                    mapping[sample_cols[i]] = sample_df.iloc[i][group_col]
            return pd.Series(mapping)
    except Exception:
        pass
    return pd.Series({c: c for c in sample_cols})


def create_volcano_plot(
    deg: pd.DataFrame,
    lfc_col: str = "logFC",
    p_col: str = "adj.P.Val",
    lfc_thresh: float = 1.0,
    p_thresh: float = 0.05,
) -> go.Figure:
    """Create an interactive volcano plot."""
    df = deg.dropna(subset=[lfc_col, p_col]).copy()
    df["minus_log10_p"] = -np.log10(df[p_col])

    # Categorize genes by significance and direction of change
    conditions = [
        (df[p_col] < p_thresh) & (df[lfc_col] > lfc_thresh),
        (df[p_col] < p_thresh) & (df[lfc_col] < -lfc_thresh),
    ]
    choices = ["Up", "Down"]
    df["category"] = np.select(conditions, choices, default="Not significant")

    fig = px.scatter(
        df,
        x=lfc_col,
        y="minus_log10_p",
        color="category",
        hover_name="Gene",
        hover_data={p_col: True, "minus_log10_p": False},
        color_discrete_map={
            "Up": "red",
            "Down": "blue",
            "Not significant": "gray",
        },
    )
    fig.update_layout(
        xaxis_title="log2 fold change",
        yaxis_title="-log10(p value)",
        legend_title="Gene category",
    )
    return fig


def create_ma_plot(
    deg: pd.DataFrame,
    avg_col: str = "AveExpr",
    lfc_col: str = "logFC",
    p_col: str = "adj.P.Val",
    lfc_thresh: float = 1.0,
    p_thresh: float = 0.05,
) -> go.Figure:
    """Create an interactive MA plot."""
    df = deg.dropna(subset=[avg_col, lfc_col, p_col]).copy()

    conditions = [
        (df[p_col] < p_thresh) & (df[lfc_col] > lfc_thresh),
        (df[p_col] < p_thresh) & (df[lfc_col] < -lfc_thresh),
    ]
    choices = ["Up", "Down"]
    df["category"] = np.select(conditions, choices, default="Not significant")

    fig = px.scatter(
        df,
        x=avg_col,
        y=lfc_col,
        color="category",
        hover_name="Gene",
        color_discrete_map={"Up": "red", "Down": "blue", "Not significant": "gray"},
    )
    fig.update_layout(
        xaxis_title="Average log2 expression",
        yaxis_title="log2 fold change",
        legend_title="Gene category",
    )
    return fig


def create_deg_heatmap(
    cpm: pd.DataFrame,
    deg: pd.DataFrame,
    group_labels: Optional[pd.Series] = None,
    p_col: str = "adj.P.Val",
    top_n: int = 50,
) -> go.Figure:
    """Interactive heatmap of logCPM values for the top N DE genes."""
    logger.debug(
        "Creating DEG heatmap using p-value column '%s' with top_n=%s",
        p_col,
        top_n,
    )

    top_genes = deg.sort_values(p_col).head(top_n)["Gene"].tolist()
    logger.debug("Top genes selected: %s", top_genes)

    genes_in_cpm = [g for g in top_genes if g in cpm.index]
    missing = set(top_genes) - set(genes_in_cpm)
    if missing:
        logger.debug("%d genes missing from CPM data: %s", len(missing), list(missing))

    heat_df = cpm.loc[genes_in_cpm]
    logger.debug("Heatmap dataframe shape before clustering: %s", heat_df.shape)

    # Cluster genes and samples for improved visualization
    if heat_df.shape[0] > 1:
        gene_order = leaves_list(linkage(pdist(heat_df), method="average"))
        heat_df = heat_df.iloc[gene_order]
        logger.debug("Gene clustering order: %s", gene_order)
    if heat_df.shape[1] > 1:
        sample_order = leaves_list(linkage(pdist(heat_df.T), method="average"))
        heat_df = heat_df.iloc[:, sample_order]
        logger.debug("Sample clustering order: %s", sample_order)

    fig = px.imshow(
        heat_df,
        labels={"x": "Sample", "y": "Gene", "color": "log2CPM"},
        zmin=heat_df.values.min(),
        zmax=heat_df.values.max(),
        color_continuous_scale="RdBu_r",
        aspect="auto",
    )
    if group_labels is not None:
        labels = group_labels.reindex(heat_df.columns)
        fig.update_xaxes(
            ticktext=labels.fillna(heat_df.columns).values,
            tickvals=list(range(len(heat_df.columns))),
        )
        logger.debug("Applied group labels to heatmap x-axis")
    fig.update_layout(xaxis_title="Sample", yaxis_title="Gene")
    return fig


def create_pca_plot(cpm: pd.DataFrame, group_labels: Optional[pd.Series] = None,
                    n_components: int = 2) -> go.Figure:
    """Perform PCA on logCPM expression values."""
    if "Gene" in cpm.columns:
        cpm = cpm.set_index("Gene")
    scaler = StandardScaler()
    scaled = scaler.fit_transform(cpm.T)
    pca = PCA(n_components=n_components)
    coords = pca.fit_transform(scaled)
    pc_df = pd.DataFrame(coords, columns=[f"PC{i+1}" for i in range(n_components)])
    pc_df.index = cpm.columns
    if group_labels is not None:
        pc_df["group"] = pc_df.index.map(group_labels).fillna(pc_df.index)
        fig = px.scatter(pc_df, x="PC1", y="PC2", color="group", hover_name=pc_df.index)
    else:
        fig = px.scatter(pc_df, x="PC1", y="PC2", hover_name=pc_df.index)
    fig.update_layout(xaxis_title="PC1", yaxis_title="PC2")
    return fig

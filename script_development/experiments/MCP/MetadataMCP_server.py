# mcp_servers/metadata_server.py
from __future__ import annotations
from mcp.server.fastmcp import FastMCP                 # <- FastMCP!
import pandas as pd, json, re, os
from typing import List, Dict, Any
from dataclasses import dataclass
from unidecode import unidecode
from Bio import Entrez

# ---------- optional context object -----------------
@dataclass
class ServerCtx:
    geo_email: str = "kevin.chen@telethonkids.org.au"

server = FastMCP("RNA-seq Metadata Tools 🚀", context_type=ServerCtx)

# ---------- helper ----------------------------------
def _clean_string(s: str) -> str:
    s = "NA" if pd.isna(s) else str(s).strip()
    s = unidecode(s).replace(" ", "_")
    return re.sub(r"[^\w]", "", s)

# ---------- TOOLS -----------------------------------
@server.tool()
def process_metadata(metadata_path: str) -> Dict[str, Any]:
    """
    Load and clean an RNA-seq metadata file (csv/tsv/txt).
    """
    if metadata_path.endswith(".csv"):
        df = pd.read_csv(metadata_path)
    else:
        df = pd.read_csv(metadata_path, sep=None, engine="python")

    df = df.loc[:, ~df.columns.str.contains("Run|Experiment", case=False)]
    df = df.drop_duplicates()
    df = df.loc[:, df.nunique() < df.shape[0]]

    df.rename(columns={c: _clean_string(c) for c in df.columns}, inplace=True)
    for c in df.columns:
        df[c] = df[c].apply(_clean_string)

    column_stats = {
        c: {"unique_count": df[c].nunique(),
            "values": df[c].unique()[:20].tolist()}
        for c in df.columns
    }
    return {
        "shape": df.shape,
        "columns": list(df.columns),
        "column_stats": column_stats,
    }

@server.tool()
def merge_analysis_columns(
    metadata_path: str,
    columns: List[str],
) -> Dict[str, str]:
    """
    Merge columns into a single `merged_analysis_group` if >1 supplied.
    """
    df = pd.read_csv(metadata_path, sep=None, engine="python")
    cols = [c for c in columns if c in df.columns]
    if not cols:
        raise ValueError("None of the requested columns exist")

    if len(cols) == 1:
        return {"merged_column": cols[0]}

    merged = df[cols].astype(str).agg("_".join, axis=1)
    out_col = "merged_analysis_group"
    df[out_col] = merged
    df.to_csv(metadata_path, index=False)              # persist change
    return {"merged_column": out_col}

@server.tool()
def extract_unique_values(metadata_path: str, column: str) -> Dict[str, Any]:
    """Return sorted unique values for `column`."""
    df = pd.read_csv(metadata_path, sep=None, engine="python")
    if column not in df.columns:
        raise ValueError(f"{column} not in metadata")
    vals = sorted(df[column].dropna().unique())
    return {"unique_values": vals, "count": len(vals)}

@server.tool()
def fetch_geo_summary(metadata_path: str, ctx: ServerCtx) -> Dict[str, str]:
    """
    Pull the GEO summary based on an accession inferred from file name.
    """
    m = re.search(r"(GSE\d+)", metadata_path)
    if not m:
        return {"message": "No GEO accession in file name"}
    acc = m.group(1)
    Entrez.email = ctx.geo_email
    search = Entrez.esearch(db="gds", term=f"{acc}[ACCN]", retmode="xml")
    ids = Entrez.read(search)["IdList"]
    if not ids:
        return {"message": f"No GEO record for {acc}"}
    summary = Entrez.read(Entrez.esummary(db="gds", id=ids[0], retmode="xml"))[0]["summary"]
    return {"geo_accession": acc, "summary": summary}

# ---------- entry-point -----------------------------
if __name__ == "__main__":
    # Explicitly specify HTTP mode and port to make debugging easier
    server.run(transport="sse", host="127.0.0.1", port=8000)

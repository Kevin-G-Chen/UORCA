"""
Helpers to create an R reproduction bundle for a single dataset analysis.

This produces a ZIP archive with:
- RNAseq.R (single, self-contained script with embedded parameters)
- metadata/edger_analysis_samples.csv (with corrected relative abundance_file paths)
- metadata/contrasts.csv (copied or synthesized when available)
- t2g.txt (included if resolvable from analysis_info or local indices)
- README.txt (usage instructions)

The resulting bundle can be executed directly via:
  Rscript RNAseq.R
"""

from __future__ import annotations

import io
import os
import json
import zipfile
from pathlib import Path
from typing import Optional, Dict, Any, List, Tuple

import pandas as pd

from ResultsIntegration import ResultsIntegrator

# Prefer reusing the same directory-zip logic as the general download package
try:
    # When imported as part of the package
    from . import add_directory_to_zip as _add_dir_to_zip
except Exception:
    # Fallback for direct import in tests
    def _add_dir_to_zip(zip_file, source_dir: str, archive_dir: str):
        if not os.path.exists(source_dir):
            return
        for root, _, files in os.walk(source_dir):
            for file in files:
                file_path = os.path.join(root, file)
                arcname = os.path.join(archive_dir, os.path.relpath(file_path, source_dir))
                try:
                    with open(file_path, 'rb') as f:
                        zip_file.writestr(arcname, f.read())
                except Exception:
                    pass


def _repo_paths() -> Dict[str, Path]:
    here = Path(__file__).resolve()
    # helpers -> streamlit_tabs -> reporting -> main_workflow
    reporting_dir = here.parents[2]
    main_workflow_dir = reporting_dir.parent
    addl_scripts = main_workflow_dir / "additional_scripts"
    repo_root = main_workflow_dir.parent
    return {
        "reporting_dir": reporting_dir,
        "main_workflow_dir": main_workflow_dir,
        "rnaseq_r": addl_scripts / "RNAseq.R",
        "repo_root": repo_root,
    }


def _first_existing(paths: List[Path]) -> Optional[Path]:
    for p in paths:
        if p and p.exists():
            return p
    return None


def _load_sample_mapping(results_dir: str, dataset_id: str) -> Optional[pd.DataFrame]:
    analysis_dir = Path(results_dir) / dataset_id
    candidate_paths = [
        analysis_dir / "metadata" / "edger_analysis_samples.csv",
        analysis_dir / "edger_analysis_samples.csv",
        Path(results_dir) / "edger_analysis_samples.csv",
        Path(results_dir).parent / "edger_analysis_samples.csv",
    ]
    f = _first_existing(candidate_paths)
    if f is None:
        return None
    # Try robust load similar to ResultsIntegration
    try:
        df_head = pd.read_csv(f, nrows=0)
        first_col = df_head.columns[0] if len(df_head.columns) else None
        index_col = 0 if (first_col == '' or (isinstance(first_col, str) and first_col.startswith('Unnamed'))) else None
        df = pd.read_csv(f, index_col=index_col)
        return df
    except Exception:
        # Fallback simple read
        return pd.read_csv(f)


def _ensure_absolute_abundance_paths(df: pd.DataFrame, results_dir: str, dataset_id: str) -> pd.DataFrame:
    df = df.copy()
    if "abundance_file" not in df.columns:
        return df
    base = Path(results_dir) / dataset_id
    abs_paths: List[str] = []
    for val in df["abundance_file"].astype(str).tolist():
        p = Path(val)
        if p.is_absolute():
            abs_paths.append(str(p))
        else:
            # Treat as relative to dataset root
            abs_paths.append(str((base / p).resolve()))
    df["abundance_file"] = abs_paths
    return df


def _synthesize_contrasts_from_analysis_info(ri: ResultsIntegrator, dataset_id: str) -> Optional[pd.DataFrame]:
    try:
        if hasattr(ri, "analysis_info") and dataset_id in ri.analysis_info:
            info = ri.analysis_info[dataset_id]
            contr = info.get("contrasts")
            if isinstance(contr, list) and contr:
                # Normalize to DataFrame with at least name/expression/description
                rows = []
                for c in contr:
                    if isinstance(c, dict) and ("name" in c and "expression" in c):
                        rows.append({
                            "name": c.get("name"),
                            "expression": c.get("expression"),
                            "description": c.get("description", ""),
                        })
                if rows:
                    return pd.DataFrame(rows)
    except Exception:
        pass
    return None


def _load_analysis_info(results_dir: str, dataset_id: str) -> Optional[Dict[str, Any]]:
    analysis_dir = Path(results_dir) / dataset_id
    candidate_paths = [
        analysis_dir / "metadata" / "analysis_info.json",
        analysis_dir / "analysis_info.json",
    ]
    f = _first_existing(candidate_paths)
    if f is None:
        return None
    try:
        with open(f, 'r') as fh:
            return json.load(fh)
    except Exception:
        return None


def _iter_files(dir_path: Path):
    for root, _, files in os.walk(dir_path):
        for fn in files:
            p = Path(root) / fn
            yield p


def _unused_placeholder() -> None:
    # Historically, a wrapper script was generated here. We now create a single
    # self-contained RNAseq.R, so this placeholder remains for backward diffs.
    return None


def create_r_reproduction_bundle(
    ri: ResultsIntegrator,
    results_dir: str,
    dataset_id: str,
    dataset_accession: str,
) -> Optional[bytes]:
    """
    Create a minimal R reproduction bundle ZIP for the given dataset.

    The bundle includes RNAseq.R, analysis_info.json, metadata CSVs, and abundance files.
    abundance_file paths are rewritten to be relative (abundance/<sample>/abundance.tsv).
    """
    try:
        paths = _repo_paths()
        rnaseq_r = paths["rnaseq_r"]

        if not rnaseq_r.exists():
            return None

        dataset_path = Path(results_dir) / dataset_id
        metadata_dir = dataset_path / "metadata"

        # Prepare in-memory zip
        zip_buffer = io.BytesIO()
        with zipfile.ZipFile(zip_buffer, mode='w', compression=zipfile.ZIP_DEFLATED) as zf:

            # Include analysis_info.json if present for transparency and primary source
            ai = _load_analysis_info(results_dir, dataset_id)
            if ai is not None:
                try:
                    zf.writestr("metadata/analysis_info.json", json.dumps(ai, indent=2))
                except Exception:
                    pass

            # Sample mapping
            sample_df = _load_sample_mapping(results_dir, dataset_id)
            if sample_df is None:
                # Try to synthesize from abundance directory
                abund_dir = dataset_path / "abundance"
                if not abund_dir.exists():
                    # Cannot proceed without sample mapping or abundance
                    return None
                rows: List[Dict[str, Any]] = []
                for sub in sorted([p for p in abund_dir.iterdir() if p.is_dir()]):
                    tsv = sub / "abundance.tsv"
                    if tsv.exists():
                        rows.append({
                            "sample": sub.name,
                            "abundance_file": str(tsv.resolve()),
                            "Group": sub.name.split("_")[0] if "_" in sub.name else "Group1",
                        })
                if not rows:
                    return None
                sample_df = pd.DataFrame(rows)

            # Ensure absolute paths for abundance_file
            sample_df = _ensure_absolute_abundance_paths(sample_df, results_dir, dataset_id)

            # Copy abundance directory using the same helper as the main package
            abund_dir = dataset_path / "abundance"
            if abund_dir.exists():
                _add_dir_to_zip(zf, str(abund_dir), "abundance")
                # Overwrite abundance_file column with relative paths based on parent dir name
                new_paths: List[str] = []
                for _, row in sample_df.iterrows():
                    abspath = Path(str(row.get("abundance_file", "")))
                    sample_dir = abspath.parent.name if abspath.parent.name else None
                    if not sample_dir:
                        sample_dir = abspath.parent.name if abspath.exists() else None
                    relp = f"abundance/{sample_dir or 'UNKNOWN'}/abundance.tsv"
                    new_paths.append(relp)
                if new_paths:
                    sample_df["abundance_file"] = new_paths
            else:
                # Fallback: copy directly from each row's abundance_file parent
                new_paths = []
                for _, row in sample_df.iterrows():
                    abspath = Path(str(row.get("abundance_file", "")))
                    src_dir = abspath.parent
                    if src_dir.exists():
                        for p in _iter_files(src_dir):
                            arcname = Path("abundance") / src_dir.name / p.relative_to(src_dir)
                            try:
                                with open(p, 'rb') as f:
                                    zf.writestr(str(arcname), f.read())
                                copied_any = True
                            except Exception:
                                pass
                        new_paths.append(f"abundance/{src_dir.name}/abundance.tsv")
                    else:
                        new_paths.append("abundance/UNKNOWN/abundance.tsv")
                if new_paths:
                    sample_df["abundance_file"] = new_paths

            # Write sample mapping into zip (with relative abundance paths)
            sample_csv_bytes = sample_df.to_csv(index=False)
            zf.writestr("metadata/edger_analysis_samples.csv", sample_csv_bytes)

            # Contrasts: copy if present, else synthesize from analysis_info
            contr_file = _first_existing([
                dataset_path / "metadata" / "contrasts.csv",
                dataset_path / "contrasts.csv",
            ])
            if contr_file and contr_file.exists():
                with open(contr_file, 'rb') as f:
                    zf.writestr("metadata/contrasts.csv", f.read())
                contr_note = f"Using contrasts from original analysis: metadata/contrasts.csv"
            else:
                # Prefer analysis_info.json if available; else fall back to in-memory ri
                contr_df = None
                if ai and isinstance(ai.get("contrasts"), list) and ai["contrasts"]:
                    rows = []
                    for c in ai["contrasts"]:
                        if isinstance(c, dict) and ("name" in c and "expression" in c):
                            rows.append({
                                "name": c.get("name"),
                                "expression": c.get("expression"),
                                "description": c.get("description", ""),
                            })
                    if rows:
                        contr_df = pd.DataFrame(rows)
                if contr_df is None:
                    contr_df = _synthesize_contrasts_from_analysis_info(ri, dataset_id)
                if contr_df is not None and len(contr_df) > 0 and {"name", "expression"}.issubset(contr_df.columns):
                    zf.writestr("metadata/contrasts.csv", contr_df.to_csv(index=False))
                    contr_note = "Synthesized contrasts from analysis_info.json"
                else:
                    contr_note = "No explicit contrasts included; the script will generate simple pairwise if feasible."

            # Try to include a t2g file (prefer analysis_info path, else organism fallback)
            t2g_note = ""
            try:
                # 1) analysis_info explicit path
                explicit_t2g = None
                if ai and isinstance(ai.get("tx2gene_file_used"), str):
                    explicit_t2g = Path(ai["tx2gene_file_used"]).resolve()
                    if explicit_t2g.exists():
                        with open(explicit_t2g, 'rb') as f:
                            zf.writestr("t2g.txt", f.read())
                        t2g_note = "Included t2g.txt from analysis_info.json (tx2gene_file_used)"

                # 2) fallback by organism
                if not t2g_note:
                    org = None
                    if ai and isinstance(ai.get("organism"), str):
                        org = ai["organism"].strip()
                    org_map = {
                        "Homo sapiens": "human",
                        "Mus musculus": "mouse",
                        "Canis lupus familiaris": "dog",
                        "Macaca mulatta": "monkey",
                        "Danio rerio": "zebrafish",
                    }
                    species = org_map.get(org, None)
                    if species:
                        t2g_path = paths["repo_root"] / "data" / "kallisto_indices" / species / "t2g.txt"
                        if t2g_path.exists():
                            with open(t2g_path, 'rb') as f:
                                zf.writestr("t2g.txt", f.read())
                            t2g_note = "Included t2g.txt inferred from organism"
                if not t2g_note:
                    t2g_note = "No t2g.txt found locally; script will run with txOut=TRUE"
            except Exception:
                t2g_note = "No t2g.txt found locally; script will run with txOut=TRUE"

            # Render single self-contained R script by patching template with fixed parameters
            group_col = None
            if ai and isinstance(ai.get("merged_column"), str) and ai["merged_column"] in sample_df.columns:
                group_col = ai["merged_column"]
            else:
                # Prefer a 'merged_column' column if present in samples
                if "merged_column" in sample_df.columns:
                    group_col = "merged_column"
                else:
                    # Heuristic fallback: choose a non-file categorical column with few uniques
                    for cn in sample_df.columns:
                        if cn.lower().endswith("_file") or cn.lower().endswith("_path") or cn in ("abundance_file", "sample", "geo_accession", "title"):
                            continue
                        vals = sample_df[cn]
                        try:
                            uniq = len(pd.unique(vals))
                        except Exception:
                            continue
                        if 2 <= uniq <= min(6, len(sample_df) - 1):
                            group_col = cn
                            break
            if not group_col:
                group_col = "Group"

            # Determine contrasts file usage in script
            use_contrasts_path = "metadata/contrasts.csv" if _first_existing([dataset_path / "metadata" / "contrasts.csv", dataset_path / "contrasts.csv"]) or (ai and ai.get("contrasts")) else None

            # Read RNAseq.R template and patch argument parsing into fixed assignments
            template_text = rnaseq_r.read_text(encoding='utf-8')
            fixed_block = (
                f"metadata_file <- \"metadata/edger_analysis_samples.csv\"\n"
                f"merged_group <- \"{group_col}\"\n"
                f"output_dir <- \".\"\n"
                f"tx2gene_file <- if (file.exists(\"t2g.txt\")) \"t2g.txt\" else \"NA\"\n"
                + ("contrasts_file <- \"metadata/contrasts.csv\"\n" if use_contrasts_path else "contrasts_file <- NULL\n")
            )

            # Replace the args parsing section with our fixed assignments
            # Use simple string replacement based on known lines in template
            patched_text = template_text
            # Remove/replace the lines where args are defined
            patched_text = patched_text.replace(
                "args <- commandArgs(trailingOnly = TRUE)\n\nmetadata_file <- args[1]\nmerged_group <- args[2]\noutput_dir <- args[3]\ntx2gene_file <- args[4]\ncontrasts_file <- if(length(args) > 4) args[5] else NULL\n",
                fixed_block,
            )
            # Fallback: if the exact block isn't matched, append fixed block after library section
            if fixed_block not in patched_text:
                insert_after = "suppressPackageStartupMessages(library(gplots))\n"
                if insert_after in patched_text:
                    patched_text = patched_text.replace(insert_after, insert_after + "\n" + fixed_block + "\n")
                else:
                    patched_text = fixed_block + "\n" + patched_text

            # Write single R script into zip (named RNAseq.R for direct execution)
            zf.writestr("RNAseq.R", patched_text)

            # README
            readme = (
                f"UORCA R Reproduction Bundle for {dataset_accession}\n"
                f"===============================\n\n"
                f"Contents:\n"
                f"- RNAseq.R: Single R script with embedded parameters (group/contrasts)\n"
                f"- metadata/analysis_info.json: Analysis metadata and contrasts used in original run (if available)\n"
                f"- metadata/edger_analysis_samples.csv: Sample mapping with abundance_file paths (relative)\n"
                f"- metadata/contrasts.csv: Contrast definitions (if available)\n"
                f"- abundance/<sample>/: Kallisto outputs for each sample\n\n"
                f"Usage:\n"
                f"1) Ensure R and required packages are installed: edgeR, limma, tximport, ComplexHeatmap, gplots\n"
                f"2) From this folder, run: Rscript RNAseq.R\n\n"
                f"Notes:\n"
                f"- {contr_note}\n"
                f"- {t2g_note}.\n"
                f"- All abundance_file paths are relative to this folder (abundance/<sample>/abundance.tsv).\n"
            )
            zf.writestr("README.txt", readme)

        zip_buffer.seek(0)
        return zip_buffer.getvalue()

    except Exception:
        return None

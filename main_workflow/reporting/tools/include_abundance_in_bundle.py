#!/usr/bin/env python3
"""
Include abundance directories into an existing R reproduction bundle and
rewrite metadata/edger_analysis_samples.csv to use relative paths.

Usage:
  uv run python main_workflow/reporting/tools/include_abundance_in_bundle.py \
      /path/to/<ACCESSION>_R_repro_bundle \
      --source-root /UORCA_results

Notes:
- The script tries two ways to find source sample dirs for each row:
  1) Use the absolute path in 'abundance_file' if it exists (take its parent dir)
  2) Fall back to <source_root>/<accession>/abundance/<sample>
     where <sample> is the parent directory name of the abundance_file path.
"""

from __future__ import annotations

import argparse
import os
import shutil
from pathlib import Path
import pandas as pd


def detect_accession_from_bundle(bundle_dir: Path) -> str | None:
    # Try to infer from analysis_info.json if present
    ai = bundle_dir / "metadata" / "analysis_info.json"
    if ai.exists():
        try:
            import json
            with open(ai, 'r') as f:
                data = json.load(f)
            acc = data.get("accession") or data.get("dataset_accession")
            if isinstance(acc, str) and acc:
                return acc
        except Exception:
            pass
    # Try from repro directory name (ends with _R_repro_bundle)
    stem = bundle_dir.name
    if stem.endswith("_R_repro_bundle"):
        return stem.replace("_R_repro_bundle", "")
    return None


def main():
    ap = argparse.ArgumentParser(description="Include abundance dirs into R repro bundle")
    ap.add_argument("bundle_dir", type=str, help="Path to <ACCESSION>_R_repro_bundle directory")
    ap.add_argument("--source-root", type=str, required=False, default=None,
                    help="Root of original results (e.g., /UORCA_results). If omitted, absolute paths in CSV must exist.")
    args = ap.parse_args()

    bundle_dir = Path(args.bundle_dir).resolve()
    if not bundle_dir.exists():
        raise SystemExit(f"Bundle directory not found: {bundle_dir}")

    meta_csv = bundle_dir / "metadata" / "edger_analysis_samples.csv"
    if not meta_csv.exists():
        raise SystemExit(f"Bundle missing metadata/edger_analysis_samples.csv: {meta_csv}")

    df = pd.read_csv(meta_csv)
    if "abundance_file" not in df.columns:
        raise SystemExit("CSV is missing 'abundance_file' column")

    accession = detect_accession_from_bundle(bundle_dir)
    source_root = Path(args.source_root) if args.source_root else None

    target_root = bundle_dir / "abundance"
    target_root.mkdir(parents=True, exist_ok=True)

    new_paths: list[str] = []
    copied = 0
    missing = 0

    for _, row in df.iterrows():
        src_file = Path(str(row["abundance_file"]))
        sample_dir_name = src_file.parent.name if src_file.parent.name else None
        # Candidate source directories
        candidates: list[Path] = []
        if src_file.exists():
            candidates.append(src_file.parent)
        if accession and source_root and sample_dir_name:
            candidates.append(source_root / accession / "abundance" / sample_dir_name)

        src_dir = next((p for p in candidates if p.exists()), None)
        if not src_dir:
            missing += 1
            new_paths.append(f"abundance/{sample_dir_name or 'UNKNOWN'}/abundance.tsv")
            continue

        dest_dir = target_root / src_dir.name
        # Copy directory if not already copied
        if not dest_dir.exists():
            shutil.copytree(src_dir, dest_dir)
            copied += 1

        new_paths.append(f"abundance/{dest_dir.name}/abundance.tsv")

    # Update CSV with relative paths
    df["abundance_file"] = new_paths
    df.to_csv(meta_csv, index=False)

    print(f"Copied {copied} sample directories into {target_root}")
    if missing:
        print(f"Warning: {missing} samples could not be located; paths set to abundance/UNKNOWN/abundance.tsv")
        print("Tip: re-run with --source-root /path/to/UORCA_results if absolute paths are not accessible.")


if __name__ == "__main__":
    main()


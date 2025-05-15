import os
import glob
import argparse
import pandas as pd
import logging

# Set up basic logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger(__name__)

def find_contrast_csvs(root_dir: str) -> list[str]:
    """Recursively find all contrasts.csv under root_dir."""
    pattern = os.path.join(root_dir, "**", "contrasts.csv")
    matches = glob.glob(pattern, recursive=True)
    logger.info(f"Found {len(matches)} contrast CSV files")
    return matches

def read_contrast_names(contrast_csv: str) -> list[str]:
    """Read the 'name' column from a contrasts.csv."""
    df = pd.read_csv(contrast_csv)
    names = df["name"].dropna().astype(str).tolist()
    logger.debug(f"Read {len(names)} contrast names from {contrast_csv}")
    return names

def find_deg_file(root_dir: str, contrast_name: str) -> str | None:
    """Find a file named deg_{contrast_name}.csv anywhere under root_dir."""
    pattern = os.path.join(root_dir, "**", f"deg_{contrast_name}.csv")
    matches = glob.glob(pattern, recursive=True)
    return matches[0] if matches else None

def build_binary_matrix(root_dir: str, lfc_threshold: float = 0, padj_threshold: float = 1.0) -> pd.DataFrame:
    logger.info(f"Building binary matrix with |logFC| >= {lfc_threshold} and adj.P.Val <= {padj_threshold}")

    # 1) Collect all contrast names
    contrast_csvs = find_contrast_csvs(root_dir)
    contrast_names: list[str] = []
    for cfile in contrast_csvs:
        contrast_names += read_contrast_names(cfile)
    contrast_names = sorted(set(contrast_names))
    logger.info(f"Found {len(contrast_names)} unique contrast names")

    # 2) For each contrast, load the DEG file and record its genes that pass thresholds
    gene_sets: dict[str, set[str]] = {}
    all_genes: set[str] = set()
    total_genes_before_filtering = 0
    found_deg_files = 0

    for name in contrast_names:
        deg_path = find_deg_file(root_dir, name)
        if not deg_path:
            logger.warning(f"No DEG file found for contrast: {name}")
            continue

        found_deg_files += 1
        df = pd.read_csv(deg_path)
        if "Gene" not in df.columns:
            logger.error(f"'Gene' column not found in {deg_path}")
            raise KeyError(f"'Gene' column not found in {deg_path}")

        all_genes_in_file = set(df["Gene"].astype(str))
        total_genes_before_filtering += len(all_genes_in_file)

        # Filter genes by logFC and adjusted p-value thresholds
        if "logFC" in df.columns and "adj.P.Val" in df.columns:
            filtered_df = df[(abs(df["logFC"]) >= lfc_threshold) &
                            (df["adj.P.Val"] <= padj_threshold)]
            genes = set(filtered_df["Gene"].astype(str))

            logger.info(f"Contrast '{name}': {len(genes)}/{len(df)} genes pass thresholds "
                      f"(|logFC| >= {lfc_threshold}, adj.P.Val <= {padj_threshold})")
        else:
            # Fall back to old behavior if expected columns not found
            logger.warning(f"Missing logFC or adj.P.Val columns in {deg_path}. Using all genes.")
            genes = all_genes_in_file

        gene_sets[name] = genes
        all_genes |= genes

    logger.info(f"Processed {found_deg_files}/{len(contrast_names)} DEG files")
    logger.info(f"Total genes before filtering: {total_genes_before_filtering}")
    logger.info(f"Total unique genes after filtering: {len(all_genes)}")

    # 3) Build the binary presence/absence DataFrame
    matrix = pd.DataFrame(
        0,
        index=sorted(all_genes),
        columns=sorted(gene_sets.keys()),
        dtype=int,
    )
    for name, genes in gene_sets.items():
        matrix.loc[list(genes), name] = 1

    logger.info(f"Created binary matrix with {matrix.shape[0]} rows × {matrix.shape[1]} columns")
    return matrix

def main():
    p = argparse.ArgumentParser(description="Build gene×contrast binary matrix")
    p.add_argument("root_dir",
                   help="Root folder under which analysis results (contrasts.csv & deg_*.csv) live")
    p.add_argument("--output", "-o",
                   default="integration_matrix.csv",
                   help="Path to write matrix CSV")
    p.add_argument("--lfc", "-l", type=float, default=1.0,
                   help="Log fold change threshold (absolute value, default=1.0)")
    p.add_argument("--padj", "-p", type=float, default=0.05,
                   help="Adjusted p-value threshold (default=0.05)")
    p.add_argument("--verbose", "-v", action="store_true",
                   help="Enable verbose logging")
    args = p.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    logger.info(f"Starting analysis on directory: {args.root_dir}")
    logger.info(f"Using thresholds: |logFC| >= {args.lfc}, adj.P.Val <= {args.padj}")

    matrix = build_binary_matrix(args.root_dir, args.lfc, args.padj)
    matrix.to_csv(args.output)
    logger.info(f"Wrote integration matrix to: {args.output}")
    print(f"Wrote integration matrix: {args.output} "
          f"({matrix.shape[0]} genes × {matrix.shape[1]} contrasts)")

if __name__ == "__main__":
    main()

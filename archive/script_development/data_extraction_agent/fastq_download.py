#!/usr/bin/env python3
"""
sra_download.py
Download SRR runs with prefetch → fasterq-dump (robust version).
"""
import argparse, subprocess, pathlib, sys, pandas as pd, shutil, os, re

# -------------------- helpers ------------------------------------------------
def shell(cmd: str):
    """Run a shell command with printing & error propagation."""
    print("▶", cmd, file=sys.stderr, flush=True)
    subprocess.run(cmd, shell=True, check=True)

def sra_path_for(srr: str, root: pathlib.Path) -> pathlib.Path:
    """Return expected .sra file path (prefetch folder layout)."""
    p = root / srr / f"{srr}.sra"
    if p.exists():
        return p
    # fall-back: look anywhere under root (handles legacy flat layouts)
    hits = list(root.rglob(f"{srr}.sra"))
    return hits[0] if hits else p   # may be non-existent

def fastq_done(srr: str, fastq_dir: pathlib.Path):
    """True if compressed FASTQs already exist and are >0 B."""
    for r in (1, 2):      # paired
        f = fastq_dir / f"{srr}_{r}.fastq.gz"
        if not f.exists() or f.stat().st_size == 0:
            return False
    return True

# -------------------- CLI ----------------------------------------------------
ap = argparse.ArgumentParser(
    description="Download SRR runs with prefetch/fasterq-dump")
g  = ap.add_mutually_exclusive_group(required=True)
g.add_argument("--srr-file", help="file: one SRR per line")
g.add_argument("--csv",
               help="CSV that contains an SRR column (default column name SRR)")
ap.add_argument("--csv-column", default="SRR")
ap.add_argument("--outdir", required=True)
ap.add_argument("--threads", type=int, default=6)
ap.add_argument("--spots", type=int, help="-X for fasterq-dump (down-sample)")
ap.add_argument("--temp", default="/tmp")
args = ap.parse_args()

# -------------------- collect SRRs ------------------------------------------
if args.srr_file:
    srrs = [l.strip() for l in open(args.srr_file) if l.strip()]
else:
    df = pd.read_csv(args.csv)
    col = args.csv_column
    if col not in df.columns:
        sys.exit(f"Column {col} not found in {args.csv}")
    srrs = df[col].dropna().astype(str).tolist()

if not srrs:
    sys.exit("No SRR accessions provided.")
print(f"{len(srrs):,} SRR accessions queued.")

# -------------------- directory setup ---------------------------------------
outdir       = pathlib.Path(args.outdir).resolve()
prefetch_dir = outdir / "sra"
fastq_dir    = outdir / "fastq"
for d in (prefetch_dir, fastq_dir):
    d.mkdir(parents=True, exist_ok=True)

# -------------------- 1. prefetch (.sra) ------------------------------------
need_fetch = [s for s in srrs if not sra_path_for(s, prefetch_dir).exists()]
if need_fetch:
    tmp_lst = outdir / "srr_to_fetch.txt"
    tmp_lst.write_text("\n".join(need_fetch) + "\n")
    shell(f"prefetch --option-file {tmp_lst} "
          f"-O {prefetch_dir} "
          "-t http "             # force HTTPS transport
          "--progress")
else:
    print("✓ All .sra files already present – skipping prefetch.")

# -------------------- 2. fasterq-dump → FASTQ.gz -----------------------------
for srr in srrs:
    if fastq_done(srr, fastq_dir):
        print(f"✓ {srr}: FASTQ already compressed – skipping conversion.")
        continue

    sra_file = sra_path_for(srr, prefetch_dir)
    if not sra_file.exists():
        print(f"⚠ {srr}: .sra still missing, skip.", file=sys.stderr)
        continue

    cmd = ["fasterq-dump", str(sra_file),
           "--threads", str(args.threads),
           "--split-files",
           "-O", str(fastq_dir),
           "--temp", args.temp]
    if args.spots:
        cmd += ["-X", str(args.spots)]
    shell(" ".join(cmd))

    # gzip / pigz
    pigz = shutil.which("pigz")
    for fq in fastq_dir.glob(f"{srr}_*.fastq"):
        shell(f"{pigz or 'gzip'} -{'p ' + str(args.threads) if pigz else 'f'} {fq}")

print("\n✅  Finished. FASTQs stored in", fastq_dir)

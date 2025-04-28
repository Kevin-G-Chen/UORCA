#!/usr/bin/env python3
"""
sra_download.py
Download SRR runs with prefetch ➜ fasterq-dump.

USAGE
  # (a) from a plain text file with one SRR per line
  python sra_download.py --srr-file srr_ids.txt --outdir /data/fastq

  # (b) directly from the long CSV we built earlier
  python sra_download.py --csv meta_long.csv  --outdir /data/fastq \
                         --csv-column SRR

OPTIONS
  --threads N        (default 6)
  --spots  N         download only first N spots (passes -X to fasterq-dump)
  --temp  DIR        scratch space (default /tmp)
"""
import argparse, subprocess, pathlib, sys, csv, re, pandas as pd, shutil, os

# ---------- shell helper ----------------------------------------------------
def run(cmd, **kw):
    print("▶", cmd, file=sys.stderr); sys.stderr.flush()
    return subprocess.run(cmd, shell=True, check=True, **kw)

# ---------- CLI -------------------------------------------------------------
ap = argparse.ArgumentParser()
g  = ap.add_mutually_exclusive_group(required=True)
g.add_argument("--srr-file", help="file with one SRR per line")
g.add_argument("--csv", help="CSV that contains an SRR column")
ap.add_argument("--csv-column", default="SRR",
                help="column name if --csv is given (default SRR)")
ap.add_argument("--outdir", required=True, help="where FASTQ will go")
ap.add_argument("--threads", type=int, default=6)
ap.add_argument("--spots", type=int, default=None,
                help="limit reads per run (passes -X to fasterq-dump)")
ap.add_argument("--temp", default="/tmp")
args = ap.parse_args()

# ---------- gather SRR list --------------------------------------------------
if args.srr_file:
    srrs = [l.strip() for l in open(args.srr_file) if l.strip()]
else:
    df   = pd.read_csv(args.csv)
    if args.csv_column not in df.columns:
        sys.exit(f"Column {args.csv_column} not found in {args.csv}")
    srrs = df[args.csv_column].dropna().astype(str).tolist()

if not srrs:
    sys.exit("No SRR accessions found.")

print(f"{len(srrs)} SRR runs queued.")

# ---------- ensure dirs ------------------------------------------------------
outdir = pathlib.Path(args.outdir).resolve()
fastq_dir = outdir / "fastq"
prefetch_dir = outdir / "sra"
for d in (outdir, fastq_dir, prefetch_dir):
    d.mkdir(parents=True, exist_ok=True)

# ---------- 1. prefetch all .sra files --------------------------------------
srr_list_file = outdir / "srr_list.txt"
srr_list_file.write_text("\n".join(srrs) + "\n")

run(f"prefetch --option-file {srr_list_file} "
    f"-O {prefetch_dir} "
    "--progress"                                  # show per-run progress
)

# ---------- 2. convert -> FASTQ ------------------------------------------------
for srr in srrs:
    sra_path = prefetch_dir / f"{srr}.sra"
    if not sra_path.exists():
        print(f"⚠  {sra_path} missing, skipping", file=sys.stderr)
        continue
    cmd = ["fasterq-dump", str(sra_path),
           "--threads", str(args.threads),
           "--split-files",
           "-O", str(fastq_dir),
           "--temp", args.temp]
    if args.spots:
        cmd += ["-X", str(args.spots)]
    run(" ".join(cmd))

    # compress on the fly (pigz if available)
    for fq in fastq_dir.glob(f"{srr}_*.fastq"):
        run(f"pigz -p {args.threads} {fq}" if shutil.which("pigz")
            else f"gzip -f {fq}")

print("\n✅  All done. FASTQ lives in", fastq_dir)

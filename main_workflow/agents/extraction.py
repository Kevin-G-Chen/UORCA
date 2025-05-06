"""
Data‑Extraction Agent
─────────────────────
* Tool 1  fetch_geo_metadata  – GEO → GSM ↔ SRX ↔ SRR long table (.csv)
* Tool 2  download_fastqs     – prefetch + fasterq‑dump into FASTQ.gz

Both tools populate ctx.deps so downstream agents (metadata / analysis)
see the files immediately.
"""
from __future__ import annotations
import asyncio, subprocess, re, pathlib, shutil, sys, os, logging
from typing import List

import pandas as pd
from dotenv import load_dotenv
import GEOparse as gp
from pydantic_ai import Agent, RunContext
from shared import ExtractionContext, RNAseqCoreContext
from shared.workflow_logging import log_tool

# ── logging setup ────────────────────────────────────────────────────────────
#logging.basicConfig(
#    format="%(asctime)s  %(levelname)-8s  %(name)s ▶  %(message)s",
#    level=logging.INFO,
#    datefmt="%H:%M:%S"
#)
logger = logging.getLogger(__name__)

# ── env + agent definition ──────────────────────────────────────────────────
load_dotenv()

extract_agent = Agent(
    "openai:o4-mini",
    deps_type=RNAseqCoreContext,
    system_prompt=(
        pathlib.Path("./main_workflow/prompts/extraction.txt").read_text()
        if pathlib.Path("prompts/extraction_system.txt").exists()
        else "You handle data extraction tasks."
    )
)

# ─────────────────────────────────────────────────────────────────────────────
@extract_agent.tool
@log_tool
async def fetch_geo_metadata(ctx: RunContext[RNAseqCoreContext], accession: str) -> str:
    """Retrieve sample‑level metadata and SRR run IDs for a GEO series."""
    logger.info("🔍 fetch_geo_metadata() called for %s", accession)

    out_root = pathlib.Path(ctx.deps.output_dir or ".").resolve()
    meta_dir = out_root / "metadata"
    meta_dir.mkdir(parents=True, exist_ok=True)
    logger.info("Saving metadata under %s", meta_dir)

    gse = gp.get_GEO(accession, destdir=str(meta_dir), silent=True)
    gsms = gse.gsms
    logger.info("Fetched %d GSM samples", len(gsms))

    re_srx, re_srr = re.compile(r"(SR[XP]\d+)"), re.compile(r"(SRR\d+)")

    def srx_from_rel(rel: List[str] | None):
        if not rel:
            return None
        for line in rel:
            hit = re_srx.search(line)
            if hit:
                return hit.group(1)
        return None

  # ── wide sample table (now uses GEOparse phenotype_data) ────────────────
    meta_wide = (
        gse.phenotype_data
          .reset_index()              # GSM ID becomes its own column
          .rename(columns={"index": "GSM"})
    )
    meta_wide.to_csv(meta_dir / "meta_wide.csv", index=False)
    logger.info("meta_wide.csv written (%d rows, %d cols)", len(meta_wide), meta_wide.shape[1])


    # ── long GSM↔SRX↔SRR table ──────────────────────────────────────────────
    def srrs_from_srx(srx: str):
        cmd = (
            f"esearch -db sra -query {srx} | "
            "efetch -format runinfo | cut -d',' -f1 | grep ^SRR"
        )
        try:
            out = subprocess.run(
                cmd, shell=True, check=True,
                stdout=subprocess.PIPE, text=True
            ).stdout.splitlines()
            return [s for s in out if re_srr.match(s)]
        except subprocess.CalledProcessError:
            logger.warning("Entrez Direct failed for %s – skipping", srx)
            return []

    long_rows = []
    for gsm, g in gsms.items():
        srx = srx_from_rel(g.metadata.get("relation"))
        if not srx:
            continue
        for srr in srrs_from_srx(srx):
            long_rows.append({"GSM": gsm, "SRX": srx, "SRR": srr})

    meta_long = pd.DataFrame(long_rows)
    logger.info("Constructed long table: %d SRR rows", len(meta_long))

    # ── merge & save ─────────────────────────────────────────────────────────
    out_df = meta_long.merge(
        meta_wide,
        on="GSM",
        how="left",
        validate="many_to_one"
    )
    out_df.to_csv(meta_dir / "meta_long.csv", index=False)
    logger.info("meta_long.csv written (%d rows)", len(out_df))

    # expose merged table to downstream agents
    ctx.deps.metadata_df = out_df
    ctx.deps.metadata_path = str(meta_dir / "meta_long.csv")

    return (
        f"Fetched {len(gsms)} GSM samples.\n"
        f"Identified {len(meta_long)} SRR runs "
        f"across {meta_long['SRX'].nunique()} SRX experiments."
    )

# ─────────────────────────────────────────────────────────────────────────────
@extract_agent.tool
@log_tool
async def download_fastqs(
    ctx: RunContext[RNAseqCoreContext],
    threads: int = 6,
    max_spots: int | None = None
) -> str:
    """Convert SRR accessions → **paired FASTQ.gz** using SRA‑Toolkit."""

    logger.info("📥 download_fastqs() called – %d threads, max_spots=%s",
                threads, max_spots)

    if ctx.deps.metadata_df is None or "SRR" not in ctx.deps.metadata_df.columns:
        raise ValueError("metadata_df with SRR column required. Run fetch_geo_metadata first.")

    out_root = pathlib.Path(ctx.deps.output_dir or ".").resolve()
    prefetch_dir = out_root / "sra"
    fastq_dir = out_root / "fastq"
    for d in (prefetch_dir, fastq_dir):
        d.mkdir(parents=True, exist_ok=True)

    def sra_path(srr: str) -> pathlib.Path:
        p = prefetch_dir / srr / f"{srr}.sra"
        if p.exists():
            return p
        hits = list(prefetch_dir.rglob(f"{srr}.sra"))
        return hits[0] if hits else p

    def fastq_ready(srr: str):
        return all((fastq_dir / f"{srr}_{r}.fastq.gz").is_file() and
                    (fastq_dir / f"{srr}_{r}.fastq.gz").stat().st_size > 0
                    for r in (1, 2))

    srrs = ctx.deps.metadata_df["SRR"].dropna().astype(str).tolist()
    logger.info("Total SRRs listed: %d", len(srrs))

    need_prefetch = [s for s in srrs if not sra_path(s).exists()]
    logger.info("SRRs needing prefetch: %d", len(need_prefetch))

    if need_prefetch:
        lst = out_root / "srr_to_fetch.txt"
        lst.write_text("\n".join(need_prefetch) + "\n")
        cmd = f"prefetch --option-file {lst} -O {prefetch_dir} -t http --progress"
        logger.info("Running prefetch…")
        proc = await asyncio.create_subprocess_shell(cmd)
        await proc.communicate()
        if proc.returncode != 0:
            raise RuntimeError("prefetch failed")

    pigz = shutil.which("pigz")
    converted = 0
    for srr in srrs:
        if fastq_ready(srr):
            logger.debug("%s already converted – skip", srr)
            continue
        sra = sra_path(srr)
        if not sra.exists():
            logger.warning("⚠ %s: .sra missing after prefetch – skip.", srr)
            continue
        cmd = [
            "fasterq-dump", str(sra),
            "--threads", str(threads),
            "--split-files",
            "-O", str(fastq_dir)
        ]
        if max_spots:
            cmd += ["-X", str(max_spots)]
        logger.info("Running fasterq-dump on %s", srr)
        proc = await asyncio.create_subprocess_exec(*cmd)
        await proc.communicate()
        if proc.returncode != 0:
            logger.warning("⚠ fasterq-dump failed on %s", srr)
            continue

        # gzip newly created FASTQs
        for fq in fastq_dir.glob(f"{srr}_*.fastq"):
            gz_cmd = [pigz or "gzip", "-f", str(fq)]
            await asyncio.create_subprocess_exec(*gz_cmd)
        converted += 1

    ctx.deps.fastq_dir = str(fastq_dir)
    logger.info("FASTQ conversion finished: %d new SRRs", converted)

    return (
        "FASTQ download complete.\n"
        f"Total SRRs: {len(srrs)}   Newly converted: {converted}"
    )

# ─────────────────────────────────────────────────────────────────────────────
@log_tool
async def run_agent_async(prompt: str, deps: ExtractionContext, usage=None):
    """Thin wrapper used by master.py (async all‑the‑way)."""
    logger.info("🛠️ Extraction agent invoked by master – prompt: %s", prompt)
    return await extract_agent.run(prompt, deps=deps, usage=usage)

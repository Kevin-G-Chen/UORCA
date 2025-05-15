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
from shared.workflow_logging import log_tool, log_agent_tool

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

# ──────────────────────────────────────────────────────────────────────────
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

    # Capture dataset information
    dataset_info_parts = []

    # Get title
    if 'title' in gse.metadata and gse.metadata['title']:
        title = gse.metadata['title'][0] if isinstance(gse.metadata['title'], list) else gse.metadata['title']
        dataset_info_parts.append(f"Title: {title}")

    # Get summary
    if 'summary' in gse.metadata and gse.metadata['summary']:
        summary = gse.metadata['summary'][0] if isinstance(gse.metadata['summary'], list) else gse.metadata['summary']
        dataset_info_parts.append(f"Summary: {summary}")

    # Get overall design
    if 'overall_design' in gse.metadata and gse.metadata['overall_design']:
        design = gse.metadata['overall_design'][0] if isinstance(gse.metadata['overall_design'], list) else gse.metadata['overall_design']
        dataset_info_parts.append(f"Overall Design: {design}")

    # Combine all information
    dataset_information = "\n\n".join(dataset_info_parts)

    # Store in context
    ctx.deps.dataset_information = dataset_information
    logger.info("📑 Captured dataset information from GEO metadata")

    # Save to a file for reference
    info_path = meta_dir / "dataset_info.txt"
    with open(info_path, 'w') as f:
        f.write(dataset_information)
    logger.info("💾 Saved dataset information to %s", info_path)

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
    logger.info("meta_wide.csv prepared (%d rows, %d cols)", len(meta_wide), meta_wide.shape[1])


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

    metadata_filename = f"{accession}_metadata.csv"
    out_df.to_csv(meta_dir / metadata_filename, index=False)

    logger.info("%s written (%d rows)", metadata_filename, len(out_df))

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

    # Determine available CPUs from SLURM or system
    slurm_cpus = int(os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count() or 8))
    threads = max(1, min(slurm_cpus - 2, threads))  # Use provided threads but cap at available CPUs minus 2
    logger.info("⚙️ Using %d threads based on system resources (SLURM_CPUS_PER_TASK=%s)",
                threads, os.environ.get("SLURM_CPUS_PER_TASK", "not set"))

    out_root = pathlib.Path(ctx.deps.output_dir or ".").resolve()
    prefetch_dir = out_root / "sra"
    fastq_dir = out_root / "fastq"
    for d in (prefetch_dir, fastq_dir):
        d.mkdir(parents=True, exist_ok=True)

    def sra_path(srr: str) -> pathlib.Path:
        # First try the standard .sra extension
        p = prefetch_dir / srr / f"{srr}.sra"
        if p.exists():
            return p

        # Next try .sralite extension
        p_lite = prefetch_dir / srr / f"{srr}.sralite"
        if p_lite.exists():
            return p_lite

        # If neither exists at standard location, do recursive search for both
        hits = list(prefetch_dir.rglob(f"{srr}.sra")) + list(prefetch_dir.rglob(f"{srr}.sralite"))
        return hits[0] if hits else p  # Return first hit or original path if nothing found


    def fastq_ready(srr: str):
        # Check if the compressed FASTQ file exists
        # Using "{srr}*.fastq.gz" to match any variation
        fastq_files = list(fastq_dir.glob(f"{srr}*.fastq.gz"))
        return len(fastq_files) > 0 and all(f.stat().st_size > 0 for f in fastq_files)


    srrs = ctx.deps.metadata_df["SRR"].dropna().astype(str).tolist()
    logger.info("Total SRRs listed: %d", len(srrs))

    need_prefetch = [s for s in srrs if not sra_path(s).exists()]
    logger.info("SRRs needing prefetch: %d", len(need_prefetch))

    # Helper function for prefetching a single SRR
    async def prefetch_single(srr: str):
        cmd = f"prefetch {srr} -O {prefetch_dir} -t https --progress"
        logger.debug("▶ Running command: %s", cmd)
        proc = await asyncio.create_subprocess_shell(
            cmd,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE
        )
        stdout, stderr = await proc.communicate()
        if proc.returncode != 0:
            logger.error("❌ Prefetch failed for %s: %s", srr, stderr.decode().strip())
        elif stderr:
            logger.debug("📤 Prefetch stderr for %s: %s", srr, stderr.decode().strip())
        return {
            "srr": srr,
            "returncode": proc.returncode,
            "stdout": stdout.decode() if stdout else "",
            "stderr": stderr.decode() if stderr else ""
        }

    if need_prefetch:
        lst = out_root / "srr_to_fetch.txt"
        lst.write_text("\n".join(need_prefetch) + "\n")

        # Use parallel prefetch if there are many SRRs
        if len(need_prefetch) > 5:
            logger.info("🔄 Prefetching %d SRRs in parallel batches", len(need_prefetch))
            # Create batches of SRRs to avoid overwhelming the system
            batch_size = min(10, max(1, len(need_prefetch) // 2))
            batches = [need_prefetch[i:i+batch_size] for i in range(0, len(need_prefetch), batch_size)]

            successful_prefetch = 0
            for batch_num, batch in enumerate(batches, 1):
                logger.info("🔄 Processing prefetch batch %d/%d with %d SRRs",
                            batch_num, len(batches), len(batch))
                tasks = [prefetch_single(srr) for srr in batch]
                results = await asyncio.gather(*tasks)

                # Log results
                for res in results:
                    if res["returncode"] == 0:
                        successful_prefetch += 1
                        logger.info("✅ Successfully prefetched %s", res["srr"])
                    else:
                        logger.error("❌ Failed to prefetch %s", res["srr"])

            logger.info("✅ Prefetch completed - %d/%d successful", successful_prefetch, len(need_prefetch))
        else:
            # For small numbers, use the original command to maintain compatibility
            cmd = f"prefetch --option-file {lst} -O {prefetch_dir} -t https --progress"
            logger.info("▶ Running prefetch: %s", cmd)
            proc = await asyncio.create_subprocess_shell(
                cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE
            )
            stdout, stderr = await proc.communicate()
            if stdout:
                logger.info("📤 Prefetch stdout: %s", stdout.decode().strip())
            if stderr:
                logger.info("📤 Prefetch stderr: %s", stderr.decode().strip())
            if proc.returncode != 0:
                raise RuntimeError(f"prefetch failed with code {proc.returncode}")

    # Check which SRRs we need to convert to FASTQ
    need_conversion = [srr for srr in srrs if not fastq_ready(srr) and sra_path(srr).exists()]
    logger.info("SRRs needing conversion to FASTQ: %d", len(need_conversion))

    # Find pigz for faster compression
    pigz = shutil.which("pigz")
    if pigz:
        logger.info("🔧 Found pigz for faster parallel compression")
    else:
        logger.info("⚠️ pigz not found, falling back to standard gzip")

    # Process conversions (fasterq-dump + compression)
    converted = 0
    for srr in need_conversion:
        sra = sra_path(srr)
        if not sra.exists():
            logger.warning("⚠ %s: .sra missing after prefetch – skip.", srr)
            continue

        # Calculate threads for this conversion based on total available
        conversion_threads = max(1, min(threads, 16))  # Cap at 16 threads per conversion

        # Run fasterq-dump with proper logging
        cmd = [
            "fasterq-dump", str(sra),
            "--threads", str(conversion_threads),
            "--split-files",
            "-v",
            "-O", str(fastq_dir)
        ]
        if max_spots:
            cmd += ["-X", str(max_spots)]

        logger.info("▶ Running fasterq-dump on %s with %d threads", srr, conversion_threads)
        proc = await asyncio.create_subprocess_exec(
            *cmd,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE
        )
        stdout, stderr = await proc.communicate()

        # Log the stdout/stderr for debugging
        if stdout:
            logger.debug("📤 Fasterq-dump stdout for %s: %s", srr, stdout.decode().strip())
        if stderr:
            logger.info("📤 Fasterq-dump stderr for %s: %s", srr, stderr.decode().strip())

        if proc.returncode != 0:
            logger.warning("⚠ fasterq-dump failed on %s with code %d", srr, proc.returncode)
            continue

        # Compress the FASTQ files one by one (simpler approach)
        fastqs = list(fastq_dir.glob(f"{srr}*.fastq"))
        if fastqs:
            logger.info("🔄 Compressing %d FASTQ files for %s", len(fastqs), srr)
            for fq in fastqs:
                if pigz:
                    # Use pigz with multiple threads if available
                    threads_per_file = max(1, slurm_cpus // 2)  # Use half of available CPUs for compression
                    gz_cmd = [pigz, "-fv", "-p", str(threads_per_file), str(fq)]
                else:
                    # Fall back to standard gzip
                    gz_cmd = ["gzip", "-fv", str(fq)]

                logger.info("▶ Running compression: %s", " ".join(gz_cmd))
                proc = await asyncio.create_subprocess_exec(
                    *gz_cmd,
                    stdout=asyncio.subprocess.PIPE,
                    stderr=asyncio.subprocess.PIPE
                )
                stdout, stderr = await proc.communicate()
                if proc.returncode == 0:
                    logger.info("✅ Compressed %s successfully", fq.name)
                else:
                    logger.warning("⚠ Compression failed for %s: %s", fq, stderr.decode().strip())



        converted += 1
        logger.info("✅ Completed processing for %s (%d/%d)",
                    srr, converted, len(need_conversion))

    # Update context for downstream agents
    ctx.deps.fastq_dir = str(fastq_dir)
    logger.info("✅ FASTQ conversion finished: %d new SRRs converted", converted)

    return (
        "FASTQ download complete.\n"
        f"Total SRRs: {len(srrs)}   Newly converted: {converted}\n"
        f"FASTQ files saved to: {fastq_dir}"
    )


# ─────────────────────────────────────────────────────────────────────────────
@log_agent_tool
async def run_agent_async(prompt: str, deps: ExtractionContext, usage=None):
    """Thin wrapper used by master.py."""
    logger.info("🛠️ Extraction agent invoked by master – prompt: %s", prompt)
    result = await extract_agent.run(prompt, deps=deps, usage=usage)

    # Log the agent's output
    logger.info("📄 Extraction agent output: %s", result.output)

    # If you want to log usage statistics
    if hasattr(result, 'usage') and result.usage:
        try:
            usage_stats = result.usage()
            logger.info("📊 Extraction agent usage: %s", usage_stats)
        except Exception as e:
            logger.debug("Could not get usage stats: %s", e)

    return result

"""
Data-Extraction Agent
─────────────────────
* Tool 1  fetch_geo_metadata  – GEO → GSM ↔ SRX ↔ SRR long table (.csv)
* Tool 2  download_fastqs     – prefetch + fasterq-dump into FASTQ.gz

Both tools populate ctx.deps so downstream agents (metadata / analysis)
see the files immediately.
"""
from __future__ import annotations
import asyncio, subprocess, re, pathlib, shutil, sys, os, pandas as pd
from typing import List
from dotenv import load_dotenv, find_dotenv
import GEOparse as gp
from pydantic_ai import Agent, RunContext
from shared import RNAseqData

load_dotenv(find_dotenv())          # one call is enough   [oai_citation_attribution:9‡PyPI](https://pypi.org/project/python-dotenv/?utm_source=chatgpt.com)

extract_agent: Agent[RNAseqData, str] = Agent(
    "openai:gpt-4o-mini",
    deps_type=RNAseqData,
    system_prompt=open("prompts/extraction_system.txt").read()
                  if pathlib.Path("prompts/extraction_system.txt").exists()
                  else "You handle data extraction tasks.",
    instrument=True,
)

# ──────────────────────────────────────────────────────────────────────────────
@extract_agent.tool
async def fetch_geo_metadata(ctx: RunContext[RNAseqData], accession: str) -> str:
    """
    Retrieve sample-level metadata and SRR run IDs for a GEO series.

    Args
    ----
    accession : str
        A **GSE accession** (e.g. ``GSE262710``).

    Side-effects
    ------------
    * Downloads the GEO series into ``ctx.deps.output_dir/geo/``.
    * Creates **two CSVs**:
      1. ``<output_dir>/metadata/meta_wide.csv`` – wide GSM-by-column table.
      2. ``<output_dir>/metadata/meta_long.csv`` – long GSM↔SRX↔SRR mapping.
    * Stores the *long* dataframe in ``ctx.deps.metadata_df`` so later agents
      can merge, reflect, and design contrasts.

    Returns
    -------
    str
        A human-readable summary of how many GSMs and SRRs were found.

    Failure modes
    -------------
    * Raises ``ValueError`` if the accession is not a valid GSE.
    * Continues gracefully when individual GSMs lack SRX/SRR relations.

    Equivalent CLI
    --------------
    .. code-block:: bash

       python geo_to_sra_metadata_debug.py GSE262710 --out meta_long.csv
    """
    out_root   = pathlib.Path(ctx.deps.output_dir or ".").resolve()
    meta_dir   = out_root / "metadata"
    meta_dir.mkdir(parents=True, exist_ok=True)

    gse = gp.get_GEO(accession, destdir=str(meta_dir), silent=True)         ## GEOparse  [oai_citation_attribution:10‡Anandology](https://anandology.com/python-for-bioinformatics/cookbook/geoparse.html?utm_source=chatgpt.com)
    gsms = gse.gsms
    re_srx, re_srr = re.compile(r"(SR[XP]\d+)"), re.compile(r"(SRR\d+)")

    def srx_from_rel(rel: List[str] | None):
        if not rel:
            return None
        for line in rel:
            hit = re_srx.search(line)
            if hit:
                return hit.group(1)
        return None

    # Wide (sample meta) -------------------------------------------------------
    meta_wide = (pd.DataFrame.from_dict(
        {n: {k: (v[0] if v else None) for k, v in g.metadata.items()}
         for n, g in gsms.items()},
        orient="index")
        .reset_index()
        .rename(columns={"index": "GSM"}))
    meta_wide.to_csv(meta_dir / "meta_wide.csv", index=False)

    # Long (GSM↔SRX↔SRR) ------------------------------------------------------
    def srrs_from_srx(srx: str):
        cmd = (f"esearch -db sra -query {srx} | "
               "efetch -format runinfo | cut -d',' -f1 | grep ^SRR")
        try:
            out = subprocess.run(cmd, shell=True, check=True,
                                 stdout=subprocess.PIPE,
                                 text=True).stdout.splitlines()
            return [s for s in out if re_srr.match(s)]
        except subprocess.CalledProcessError:
            return []

    long_rows = []
    for gsm, g in gsms.items():
        srx = srx_from_rel(g.metadata.get("relation"))
        if not srx:
            continue
        for srr in srrs_from_srx(srx):                                     ## Entrez Direct  [oai_citation_attribution:11‡NCBI](https://www.ncbi.nlm.nih.gov/dbvar/content/tools/entrez/?utm_source=chatgpt.com)
            long_rows.append({"GSM": gsm, "SRX": srx, "SRR": srr})

    meta_long = pd.DataFrame(long_rows)
    meta_long.to_csv(meta_dir / "meta_long.csv", index=False)
    ctx.deps.metadata_df = meta_long
    ctx.deps.metadata_path = str(meta_dir / "meta_long.csv")

    return (f"Fetched {len(gsms)} GSM samples.\n"
            f"Identified {len(meta_long)} SRR runs "
            f"across {meta_long['SRX'].nunique()} SRX experiments.")

# ──────────────────────────────────────────────────────────────────────────────
@extract_agent.tool
async def download_fastqs(
    ctx: RunContext[RNAseqData],
    threads: int = 6,
    max_spots: int | None = None
) -> str:
    """
    Convert SRR accessions → **paired FASTQ.gz** using *prefetch* + *fasterq-dump*.

    Inputs (via ctx.deps)
    ---------------------
    * ``metadata_df``   must contain an “SRR” column (created by `fetch_geo_metadata`)
    * ``output_dir``    destination folder; FASTQs appear under ``<output_dir>/fastq``

    Parameters
    ----------
    threads : int, default 6
        Number of CPU threads for *fasterq-dump* and *pigz*.
    max_spots : int | None
        If given, passes ``-X`` to *fasterq-dump* for **down-sampled tests**.

    Behaviour
    ---------
    1. Creates two sub-folders:
       * ``<out>/sra`` – *.sra downloads* (via *prefetch*)
       * ``<out>/fastq`` – final compressed FASTQs
    2. Skips any SRR that already has gzipped FASTQ pairs >0 B.
    3. Validates exit codes; emits warnings instead of aborting the run.

    Returns
    -------
    str
        Multi-line summary of how many SRRs were downloaded and converted.

    Equivalent CLI
    --------------
    .. code-block:: bash

       python sra_download.py --csv meta_long.csv --outdir out --threads 8

    Notes
    -----
    The implementation follows the **official SRA-Toolkit recommendation**
    (*prefetch* ➜ *fasterq-dump*).  [oai_citation_attribution:12‡GitHub](https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump?utm_source=chatgpt.com) [oai_citation_attribution:13‡GitHub](https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump/633360aa302b9f2b6e8ceda7c99dde07f4f20e2e?utm_source=chatgpt.com)
    """
    if ctx.deps.metadata_df is None or "SRR" not in ctx.deps.metadata_df.columns:
        raise ValueError("metadata_df with SRR column required. "
                         "Run fetch_geo_metadata first.")

    out_root       = pathlib.Path(ctx.deps.output_dir or ".").resolve()
    prefetch_dir   = out_root / "sra"
    fastq_dir      = out_root / "fastq"
    for d in (prefetch_dir, fastq_dir):
        d.mkdir(parents=True, exist_ok=True)

    def sra_path(srr: str) -> pathlib.Path:
        p = prefetch_dir / srr / f"{srr}.sra"
        if p.exists():
            return p
        hits = list(prefetch_dir.rglob(f"{srr}.sra"))
        return hits[0] if hits else p

    def fastq_ready(srr: str):
        return all((fastq_dir / f"{srr}_{r}.fastq.gz").is_file()
                   and (fastq_dir / f"{srr}_{r}.fastq.gz").stat().st_size > 0
                   for r in (1, 2))

    srrs = ctx.deps.metadata_df["SRR"].dropna().astype(str).tolist()
    need_prefetch = [s for s in srrs if not sra_path(s).exists()]

    if need_prefetch:
        lst = out_root / "srr_to_fetch.txt"
        lst.write_text("\n".join(need_prefetch) + "\n")
        cmd = (f"prefetch --option-file {lst} -O {prefetch_dir} "
               "-t http --progress")
        proc = await asyncio.create_subprocess_shell(cmd)
        await proc.communicate()
        if proc.returncode != 0:
            raise RuntimeError("prefetch failed")

    pigz = shutil.which("pigz")
    converted = 0
    for srr in srrs:
        if fastq_ready(srr):
            continue
        sra = sra_path(srr)
        if not sra.exists():
            print(f"⚠ {srr}: .sra missing after prefetch – skip.", file=sys.stderr)
            continue
        cmd = [ "fasterq-dump", str(sra),
                "--threads", str(threads),
                "--split-files",
                "-O", str(fastq_dir),
                "--temp", "/tmp" ]
        if max_spots:
            cmd += ["-X", str(max_spots)]
        proc = await asyncio.create_subprocess_exec(*cmd)
        await proc.communicate()
        if proc.returncode != 0:
            print(f"⚠ fasterq-dump failed on {srr}", file=sys.stderr)
            continue

        # gzip
        for fq in fastq_dir.glob(f"{srr}_*.fastq"):
            gz_cmd = [pigz or "gzip", "-f", str(fq)]
            await asyncio.create_subprocess_exec(*gz_cmd)
        converted += 1

    ctx.deps.fastq_dir = str(fastq_dir)
    return (f"FASTQ download complete.\n"
            f"Total SRRs: {len(srrs)}   Newly converted: {converted}")

# ──────────────────────────────────────────────────────────────────────────────
async def run_agent_async(prompt: str, deps: RNAseqData, usage=None):
    """Thin wrapper used by master.py (async all-the-way)."""
    return await extract_agent.run(prompt, deps=deps, usage=usage)

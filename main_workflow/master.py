from pydantic_ai import Agent, RunContext
from dotenv import load_dotenv, find_dotenv
import os, pathlib

from shared import RNAseqData
from agents import extraction, analysis, reporting

# ── load env once ───────────────────────────
load_dotenv(find_dotenv())                      # walks up dirs  [oai_citation_attribution:5‡Stack Overflow](https://stackoverflow.com/questions/41125023/how-to-use-python-dotenvs-find-dotenv-method?utm_source=chatgpt.com)
assert os.getenv("OPENAI_API_KEY"), "API key missing"

# ── master agent definition ─────────────────
master = Agent(
    "openai:o4-mini",
    system_prompt=(
        "You are the pipeline orchestrator.\n"
        "Available tools:\n"
        "• extract(acc) – download/locate data\n"
        "• meta()       – process metadata + contrasts\n"
        "• analyse()    – quant + DEG\n"
        "• report()     – final HTML/PDF\n\n"
        "Plan which tools to call – skip ones that aren't needed."
    ),
    instrument=True,
)

# ── tool wrappers (all async, all await) ────
@master.tool
async def extract(ctx: RunContext[RNAseqData], accession: str) -> str:
    """
    High-level extraction wrapper.
    The LLM will decide whether to call only metadata fetching, only FASTQ
    downloading, or both.
    """
    r = await extraction.run_agent_async(
        f"Prepare raw data for {accession}. "
        "If FASTQs exist, just fetch GEO metadata; otherwise download files.",
        deps=ctx.deps, usage=ctx.usage
    )
    return r.output

@master.tool
async def analyse(ctx: RunContext[RNAseqData]) -> str:
    r = await analysis.run_agent_async(
        "Run Kallisto and DEG (skip GSEA)", deps=ctx.deps, usage=ctx.usage
    )
    return r.output

@master.tool
async def report(ctx: RunContext[RNAseqData]) -> str:
    png_dir   = ctx.deps.output_dir + "/figures"
    rst_dir   = ctx.deps.output_dir + "/rst"
    stamp     = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    sphinx_dir = f"{ctx.deps.output_dir}/{stamp}"
    log_file   = f"{sphinx_dir}/sphinx_build.log"
    rep_deps = ReportContext(
          png_dir=png_dir,
          rst_folder=rst_dir,
          sphinx_output_folder=sphinx_dir,
          log_path=log_file
      )
    r = await reporting.run_agent_async(
        "Generate publication-ready report", deps=ctx.deps, usage=ctx.usage
    )
    return r.output

# ── CLI entrypoint ──────────────────────────
if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--accession",
        required=True)
    ap.add_argument("--output_dir",
        default=None,
        help="Root directory for all outputs")
    ap.add_argument("--kallisto_index_dir",
        default="./resources/kallisto_indices",
        help="Folder that contains one or more *.idx files")
    ap.add_argument("--tx2gene",
        default=None,
        help="Path to a transcript‑to‑gene mapping file (optional)")
    args = ap.parse_args()

    output_dir = args.output_dir
    if not output_dir:
        output_dir = f"./analysis/{args.accession}"

    deps = RNAseqData(
        fastq_dir=f"{output_dir}/fastq",
        metadata_path=f"{output_dir}/metadata/meta_long.csv",
        kallisto_index_dir=args.kallisto_index_dir,
        tx2gene_path=args.tx2gene,
        output_dir=output_dir,
    )
    pathlib.Path(deps.output_dir).mkdir(parents=True, exist_ok=True)

    run = master.run_sync(f"Analyse {args.accession} by first extracting the data, performing an analysis on it, then finally generating a report", deps=deps)
    print(run.output)
    print("Token usage:", run.usage())

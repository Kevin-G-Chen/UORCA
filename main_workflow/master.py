from pydantic_ai import Agent, RunContext
from dotenv import load_dotenv, find_dotenv
import os, pathlib

from shared import RNAseqData
from agents import extraction, metadata, analysis, reporting

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
async def meta(ctx: RunContext[RNAseqData]) -> str:
    r = await metadata.run_agent_async(
        "Clean metadata, build contrasts, maybe reflect.",
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
    r = await reporting.run_agent_async(
        "Generate publication-ready report", deps=ctx.deps, usage=ctx.usage
    )
    return r.output

# ── CLI entrypoint ──────────────────────────
if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--accession", required=True)
    args = ap.parse_args()

    deps = RNAseqData(
        fastq_dir=f"./data/{args.accession}/fastq",
        metadata_path=f"./data/{args.accession}/{args.accession}_metadata.csv",
        kallisto_index_dir="./resources/kallisto_indices",
        output_dir=f"./analysis/{args.accession}",
    )
    pathlib.Path(deps.output_dir).mkdir(parents=True, exist_ok=True)

    run = master.run_sync(f"Analyse {args.accession}", deps=deps)
    print(run.output)
    print("Token usage:", run.usage())

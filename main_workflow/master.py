from pydantic_ai import Agent, RunContext
from dotenv import load_dotenv, find_dotenv
import os, pathlib
import datetime
from shared import ExtractionContext, AnalysisContext, ReportingContext, RNAseqCoreContext
from agents import extraction, analysis, reporting
import argparse

# ── load env once ───────────────────────────
load_dotenv(find_dotenv())
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
async def extract(ctx: RunContext[ExtractionContext], accession: str) -> str:
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
async def analyse(ctx: RunContext[AnalysisContext]) -> str:
    r = await analysis.run_agent_async(
        "Run Kallisto and DEG (skip GSEA)", deps=ctx.deps, usage=ctx.usage
    )
    return r.output

@master.tool
async def report(ctx: RunContext[ReportingContext]) -> str:
    r = await reporting.run_agent_async(
        "Generate publication-ready report", deps=ctx.deps, usage=ctx.usage
    )
    return r.output

# ── CLI entrypoint ──────────────────────────
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--accession", required=True)
    ap.add_argument("--output_dir", default=None)
    ap.add_argument("--kallisto_index_dir", default="./resources/kallisto_indices")
    ap.add_argument("--tx2gene", default=None)
    args = ap.parse_args()

    output_dir = args.output_dir or f"./analysis/{args.accession}"

    # Prepare a generic context (Core), enriched as the pipeline runs.
    ctx = RNAseqCoreContext(
        accession=args.accession,
        output_dir=output_dir,
        fastq_dir=f"{output_dir}/fastq",
        metadata_path=f"{output_dir}/metadata/meta_long.csv",
        kallisto_index_dir=args.kallisto_index_dir,
        organism="human",
        tx2gene_path=args.tx2gene
    )

    # START THE AGENTIC PIPELINE: Let agent decide what tools to call!
    initial_prompt = (
        f"Analyse {args.accession} by first extracting the data, "
        f"performing an analysis on it, then finally generating a report. "
        f"Skip any steps if the data already exists and is up to date. "
        f"Document each tool invocation and output."
    )
    run = master.run_sync(
        initial_prompt,
        deps=ctx
    )

    print(run.output)
    try:
        print("Token usage:", run.usage())
    except Exception:
        pass

if __name__ == "__main__":
    main()

from shared.workflow_logging import setup_logging

import logging
from pydantic_ai import Agent, RunContext
from dotenv import load_dotenv, find_dotenv
import os, pathlib
import datetime
from shared import ExtractionContext, AnalysisContext, ReportingContext, RNAseqCoreContext
import argparse

# ── load env once ───────────────────────────
load_dotenv(find_dotenv())
assert os.getenv("OPENAI_API_KEY"), "API key missing"

# ── master agent definition ─────────────────

# read in system prompt

master_prompt_path = "./main_workflow/prompts/master.txt"
try:
    master_prompt = pathlib.Path(master_prompt_path).read_text()
except Exception as e:
    logger.warning("Could not read analysis system prompt: %s – using fallback", e)
    master_prompt = """
    You are a bioinformatics expert who oversees the execution of a bioinformatic analysis. You will not need to perform any of the analysis yourself, but instead have an expert team of specialised agents who will perform the analysis for you.
    """

master = Agent(
    "openai:o4-mini",
    system_prompt=master_prompt,
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
        f"""
        Prepare raw data for {accession} using the following steps:

        1. First, fetch metadata from the GEO database using the fetch_geo_metadata tool:
           - Download the GEO Series SOFT file
           - Parse all GSM records
           - Derive the chain GSM → SRX → SRR

        2. Then, download FASTQ files using the download_fastqs tool:
           - Use the prefetch → fasterq-dump pipeline
           - Ensure proper thread allocation for SLURM jobs
           - Compress output files using pigz if available

        Be attentive to resource usage, logging, and error handling during this process.
        """,
        deps=ctx.deps, usage=ctx.usage
    )
    return r.output

@master.tool
async def analyse(ctx: RunContext[AnalysisContext]) -> str:
    r = await analysis.run_agent_async(
        f"""
        Perform RNA-seq analysis on data from organism: {ctx.deps.organism}

        Follow these sequential steps:

        1. Identify the appropriate Kallisto index file for {ctx.deps.organism} in the resource directory
           - Ensure you select the correct species-specific index (.idx file)

        2. Run Kallisto quantification:
           - Use the FASTQ files from the extraction step
           - Optimize threads for parallelization in the SLURM environment
           - Generate abundance files

        3. Process metadata:
           - Use the metadata agent to identify biologically relevant columns
           - Merge columns if necessary into a single grouping variable
           - Extract unique values and construct appropriate contrast matrices

        4. Perform differential expression analysis using edgeR/limma:
           - Use the appropriate tx2gene file for {ctx.deps.organism}
           - Follow best practices for normalization, filtering, and statistical analysis
           - Generate visualizations (MDS plots, heatmaps, volcano plots)

        Provide detailed explanations for each step and handle errors gracefully.
        Skip GSEA analysis for this run.""",
        deps=ctx.deps,
        usage=ctx.usage
    )
    return r.output

@master.tool
async def report(ctx: RunContext[ReportingContext]) -> str:
    if not isinstance(ctx.deps, ReportingContext):
        ctx.deps = ReportingContext(**ctx.deps.dict())  # copies all existing attrs
    if not ctx.deps.png_dir:
        # Look for plot directories in standard locations created by analysis agent
        plot_dir = f"{ctx.deps.output_dir}/plots"
        ctx.deps.png_dir = plot_dir

    if not ctx.deps.rst_folder:
        ctx.deps.rst_folder = f"{ctx.deps.output_dir}/report/rst"

    if not ctx.deps.sphinx_output_folder:
        ctx.deps.sphinx_output_folder = f"{ctx.deps.output_dir}/report/sphinx"

    if not ctx.deps.log_path:
        ctx.deps.log_path = f"{ctx.deps.output_dir}/report/sphinx_build.log"

    os.makedirs(ctx.deps.rst_folder, exist_ok=True)
    os.makedirs(os.path.dirname(ctx.deps.log_path), exist_ok=True)

    r = await reporting.run_agent_async(
        "Generate a report", deps=ctx.deps, usage=ctx.usage
    )
    return r.output

# ── CLI entrypoint ──────────────────────────
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--accession", required=True)
    ap.add_argument("--output_dir", default=None)
    ap.add_argument("--resource_dir", default=None)
    ap.add_argument("--organism", default="human")
    args = ap.parse_args()

    # -------- create & remember the chosen run folder -----------------------
    output_dir = f"{args.output_dir}/{args.accession}" or f"./analysis/{args.accession}"
    os.makedirs(output_dir, exist_ok=True)

    # -------- configure logging *inside* that folder ------------------------
    log_dir = pathlib.Path(output_dir) / "logs"
    log_path = setup_logging(log_dir)
    logger = logging.getLogger(__name__)  # Add this line to get a logger
    logger.info("🚀 Starting UORCA master agent - logging to %s", log_path)

    import importlib, sys

    logger.info("📚 Importing agent modules")
    globals()["extraction"] = importlib.import_module("agents.extraction")
    globals()["analysis"]   = importlib.import_module("agents.analysis")
    globals()["reporting"]  = importlib.import_module("agents.reporting")
    logger.info("✅ All agent modules imported successfully")

    # -------- build the initial CoreContext and run the orchestrator --------
    ctx = RNAseqCoreContext(
        accession=args.accession,
        output_dir=output_dir,
        organism=args.organism,
        resource_dir=args.resource_dir,
    )

    logger.info("🧩 Built initial context with accession: %s", args.accession)

    initial_prompt = (
        f"Analyse {args.accession}, an RNAseq dataset for looking at the {args.organism} species, by first extracting the data, performing an analysis on it, then finally generating a report. Skip any steps if the data already exists and is up to date. Document each tool invocation and output."
    )

    logger.info("🤖 Running master agent with prompt: %s", initial_prompt)
    run = master.run_sync(initial_prompt, deps=ctx)

    logger.info("✅ Master agent completed execution")
    logger.info("📝 Final output: %s", run.output)

    try:
        usage_stats = run.usage()
        logger.info("📊 Token usage: %s", usage_stats)
    except Exception as e:
        logger.warning("⚠️ Could not get token usage: %s", e)


if __name__ == "__main__":
    main()

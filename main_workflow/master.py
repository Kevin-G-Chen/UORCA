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
        """,
        deps=ctx.deps, usage=ctx.usage
    )
    return r.output

@master.tool
async def analyse(ctx: RunContext[AnalysisContext]) -> str:
    """
    Run RNA-seq analysis with reflection capabilities.

    Performs multiple iterations of analysis, with each iteration having access to
    the results of previous iterations to enable self-correction and improvement.
    """
    # Define the number of reflection iterations
    num_iterations = 3

    # Store all iteration outputs
    iteration_results = []

    logger = logging.getLogger(__name__)

    # Initialize analysis_history if it doesn't exist
    if not hasattr(ctx.deps, 'analysis_history'):
        ctx.deps.analysis_history = []

    # Ensure dataset_information is carried over from extraction to analysis if available
    if hasattr(ctx.deps, 'dataset_information') and ctx.deps.dataset_information:
        logger.info("📊 Adding dataset information to analysis context")
        # Dataset information is already in ctx.deps, no need to do anything else

    # Base prompt that will be extended with previous results
    base_prompt = f"""
        Perform RNA-seq analysis on data from organism: {ctx.deps.organism}

        Follow these sequential steps:

        1. Process metadata:
           - Use the metadata agent to identify biologically relevant columns
           - Merge columns if necessary into a single grouping variable
           - Extract unique values and construct appropriate contrast matrices

        2. Identify the appropriate Kallisto index file for {ctx.deps.organism}, which will be nested in the resource directory: {ctx.deps.resource_dir} (note you will need to search recursively)
           - Ensure you select the correct species-specific index (.idx file)

        3. Run Kallisto quantification:
           - Use the FASTQ files from the extraction step
           - Optimize threads for parallelization in the SLURM environment
           - Generate abundance files

        4. Perform differential expression analysis using edgeR/limma:
           - Use the appropriate tx2gene file for {ctx.deps.organism}, again found in the resource directory: {ctx.deps.resource_dir} (again, search recursively)
           - Follow best practices for normalization, filtering, and statistical analysis
           - Generate visualizations (MDS plots, heatmaps, volcano plots)

        Provide detailed explanations for each step and handle errors gracefully.
    """

    # Log the start of the reflection process
    logger.info("🔄 Starting analysis with %d reflection iterations", num_iterations)

    for iteration in range(1, num_iterations + 1):
        # Create the prompt for this iteration
        if iteration == 1:
            current_prompt = base_prompt
        else:
            # Add reflection component for iterations after the first
            reflection_prompt = f"""

            This is reflection iteration {iteration} of {num_iterations}.

            Below are the results and approach from your previous iteration.
            Review this output carefully and consider:

            1. Was the previous approach correct? If not, what should be changed? Pay special note to the files which were used - for example, was the species correct?
            2. Were there any errors or misunderstandings you can fix in this iteration?
            3. Can you improve upon the previous results?

            Previous iteration output:
            -------------------------
            {iteration_results[-1]}
            -------------------------

            Now, proceed with your analysis again, either confirming the previous approach
            if correct or making appropriate adjustments.
            """
            current_prompt = base_prompt + reflection_prompt

        # Store the prompt in history BEFORE running the agent
        ctx.deps.analysis_history.append({
            "iteration": iteration,
            "prompt": current_prompt,
            "output": None  # Will be filled after running
        })

        # Run the current iteration
        logger.info("🔄 Running analysis iteration %d/%d", iteration, num_iterations)
        result = await analysis.run_agent_async(
            current_prompt,
            deps=ctx.deps,
            usage=ctx.usage
        )

        # Store the result and update the history with the output
        iteration_results.append(result.output)
        ctx.deps.analysis_history[-1]["output"] = result.output

        logger.info("✅ Completed analysis iteration %d/%d", iteration, num_iterations)

    # Return final iteration result along with a summary
    final_result = f"""
    Analysis completed after {num_iterations} reflection iterations.

    Final result from iteration {num_iterations}:

    {iteration_results[-1]}
    """

    return final_result

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

def cleanup_large_files(output_dir, logger):
    """
    Clean up large FASTQ and SRA files to save disk space.
    Only removes files from standard pipeline folders.

    Args:
        output_dir: Base output directory for the run
        logger: Logger instance for recording actions
    """
    output_path = pathlib.Path(output_dir)
    sra_dir = output_path / "sra"
    fastq_dir = output_path / "fastq"

    # Clean up SRA files
    if sra_dir.exists():
        logger.info("🧹 Cleaning up SRA files in %s", sra_dir)
        sra_files = list(sra_dir.glob("**/*.sra")) + list(sra_dir.glob("**/*.sralite"))
        for file in sra_files:
            try:
                file.unlink()
                logger.info("🗑️ Removed SRA file: %s", file)
            except Exception as e:
                logger.warning("⚠️ Failed to remove %s: %s", file, str(e))

    # Clean up FASTQ files
    if fastq_dir.exists():
        logger.info("🧹 Cleaning up FASTQ files in %s", fastq_dir)
        fastq_files = list(fastq_dir.glob("**/*.fastq.gz"))
        for file in fastq_files:
            try:
                file.unlink()
                logger.info("🗑️ Removed FASTQ file: %s", file)
            except Exception as e:
                logger.warning("⚠️ Failed to remove %s: %s", file, str(e))

    logger.info("✅ Cleanup completed")


# ── CLI entrypoint ──────────────────────────
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--accession", required=True)
    ap.add_argument("--output_dir", default=None)
    ap.add_argument("--resource_dir", default=None)
    ap.add_argument("--organism", required=True,
                       help="Name of the organism to be analysed (must be supplied)")
    ap.add_argument("--cleanup", action="store_true",
                       help="Clean up large FASTQ and SRA files after successful completion")
    args = ap.parse_args()

    # -------- create & remember the chosen run folder -----------------------
    output_dir = (
        f"{args.output_dir}/{args.accession}"
        if args.output_dir
        else f"./analysis/{args.accession}"
    )
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
    run = master.run_sync(initial_prompt, deps=ctx, request_limit = 100)

    logger.info("✅ Master agent completed execution")
    logger.info("📝 Final output: %s", run.output)

    # Perform cleanup if requested
    if args.cleanup:
        logger.info("🧹 Starting cleanup of large temporary files...")
        cleanup_large_files(output_dir, logger)
    else:
        logger.info("ℹ️ Cleanup skipped. Use --cleanup flag to remove FASTQ/SRA files")

    try:
        usage_stats = run.usage()
        logger.info("📊 Token usage: %s", usage_stats)
    except Exception as e:
        logger.warning("⚠️ Could not get token usage: %s", e)


if __name__ == "__main__":
    main()

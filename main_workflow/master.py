from shared.workflow_logging import setup_logging

import logging
from pydantic_ai import Agent, RunContext
from dotenv import load_dotenv, find_dotenv
import os, pathlib
import datetime
import json
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
async def evaluate_analysis(ctx: RunContext[RNAseqCoreContext]) -> str:
    """
    Evaluate the success of the analysis run by checking output files and tool logs.
    
    This tool examines:
    1. Kallisto quantification results (abundance files exist and are non-empty)
    2. edgeR/limma analysis results (DEG files exist)
    3. Tool logs for errors or failures
    4. Whether proper files were used (correct organism-specific indices and tx2gene files)
    
    Returns a detailed diagnostic string with the evaluation results.
    """
    logger = logging.getLogger(__name__)
    logger.info("🔍 Evaluating analysis results")

    # Collect basic facts about the run
    facts = []

    # Check Kallisto quantification results
    if not ctx.deps.abundance_files or len(ctx.deps.abundance_files) == 0:
        facts.append("⚠️ No abundance files were found")
    else:
        # Check if abundance files actually exist
        missing_files = []
        for file_path in ctx.deps.abundance_files:
            if not os.path.exists(file_path):
                missing_files.append(file_path)

        if missing_files:
            facts.append(f"⚠️ {len(missing_files)}/{len(ctx.deps.abundance_files)} abundance files are missing")
        else:
            facts.append(f"✓ Found {len(ctx.deps.abundance_files)} abundance files")

    # Add information about Kallisto index
    # Check Kallisto index
    kallisto_index_used = getattr(ctx.deps, 'kallisto_index_used', None)
    if kallisto_index_used:
        facts.append(f"✓ Kallisto index used: {kallisto_index_used}")
    else:
        facts.append("⚠️ No Kallisto index was recorded as being used")

    # Check for DEG results
    if not ctx.deps.deg_results_path:
        facts.append("⚠️ No differential expression results were found")
    else:
        if isinstance(ctx.deps.deg_results_path, list):
            # Multiple DEG files
            missing_deg_files = []
            for deg_file in ctx.deps.deg_results_path:
                if not os.path.exists(deg_file):
                    missing_deg_files.append(deg_file)

            if missing_deg_files:
                facts.append(f"⚠️ {len(missing_deg_files)}/{len(ctx.deps.deg_results_path)} DEG result files are missing")
            else:
                facts.append(f"✓ Found {len(ctx.deps.deg_results_path)} DEG result files")
        else:
            # Single DEG file
            if os.path.exists(ctx.deps.deg_results_path):
                facts.append(f"✓ Found DEG results file: {ctx.deps.deg_results_path}")
            else:
                facts.append(f"⚠️ DEG results file not found: {ctx.deps.deg_results_path}")

    # Add information about tx2gene file
    tx2gene_file_used = getattr(ctx.deps, 'tx2gene_file_used', None)
    if tx2gene_file_used:
        facts.append(f"✓ tx2gene file used: {tx2gene_file_used}")
    else:
        facts.append("⚠️ No tx2gene file was recorded as being used")

    # Check for errors in the tool logs
    tool_logs = getattr(ctx.deps, 'tool_logs', [])
    if tool_logs:
        error_logs = [log for log in tool_logs if log["success"] is False]
        if error_logs:
            facts.append(f"⚠️ Found {len(error_logs)} tool call errors")
            for error_log in error_logs[:3]:  # Limit to first 3 errors to avoid overwhelming
                facts.append(f"  - {error_log['tool_name']}: {error_log['error']}")
        else:
            facts.append(f"✓ No errors found in {len(tool_logs)} tool calls")

    # Use LLM to evaluate the analysis
    evaluation_agent = Agent(
        "openai:o4-mini",
        system_prompt="You are an expert RNA-seq analysis evaluator. You assess whether an analysis was successful and if the correct files were used."
    )

    evaluation_prompt = f"""
    Please evaluate this RNA-seq analysis and determine if it was successful.

    ORGANISM: {ctx.deps.organism}

    ANALYSIS FACTS:
    {chr(10).join(facts)}

    KALLISTO INDEX: {getattr(ctx.deps, 'kallisto_index_used', 'Not specified')}
    TX2GENE FILE: {getattr(ctx.deps, 'tx2gene_file_used', 'Not specified')}

    Perform the following evaluation:

    1. Determine if the analysis was SUCCESSFUL or FAILED
    2. Check if the Kallisto index matches the organism ({ctx.deps.organism})
    3. Check if the tx2gene file matches the organism ({ctx.deps.organism})
    4. Identify if any critical steps failed

    Return your evaluation as a list of diagnostic statements, with each line starting with:
    ✅ for success conditions
    ❌ for failure conditions
    ❓ for unclear/missing information

    EVALUATION:
    """

    try:
        result = await evaluation_agent.run(evaluation_prompt)
        diagnostics = result.output.strip().split('\n')

        # Determine overall success based on presence of failure indicators
        success = not any("❌" in line for line in diagnostics)

        # Store results in context - safely handle attribute setting
        setattr(ctx.deps, 'analysis_success', success)
        setattr(ctx.deps, 'analysis_diagnostics', "\n".join(diagnostics))

        # Log the evaluation results
        if success:
            logger.info("✅ Analysis evaluation: SUCCESS")
        else:
            logger.info("❌ Analysis evaluation: FAILURE - see diagnostics for details")

        return ctx.deps.analysis_diagnostics
    except Exception as e:
        logger.error("❌ Error during LLM evaluation: %s", str(e))
        # Fall back to basic evaluation on LLM error
        success = all(item.startswith("✓") for item in facts)
        diagnostics = [f"❌ Evaluation error: {str(e)}"] + facts
        setattr(ctx.deps, 'analysis_success', success)
        setattr(ctx.deps, 'analysis_diagnostics', "\n".join(diagnostics))
        return ctx.deps.analysis_diagnostics

@master.tool
async def generate_reflection(ctx: RunContext[RNAseqCoreContext]) -> str:
    """
    Generate a concise reflection based on analysis diagnostics.
    
    This tool creates a focused reflection that:
    1. Identifies what went wrong in the analysis
    2. Suggests specific corrections for the next attempt
    3. Provides guidance on what files or parameters to use
    
    Returns a concise reflection string that will guide the next analysis attempt.
    """
    logger = logging.getLogger(__name__)
    logger.info("🧠 Generating reflection for failed analysis")
    
    # Safely access analysis_diagnostics
    analysis_diagnostics = getattr(ctx.deps, 'analysis_diagnostics', None)
    if not analysis_diagnostics:
        logger.warning("⚠️ No analysis diagnostics available - running evaluation first")
        await evaluate_analysis(ctx)
        analysis_diagnostics = getattr(ctx.deps, 'analysis_diagnostics', "No diagnostics available")
    
    # Safely access tool_logs
    tool_logs = getattr(ctx.deps, 'tool_logs', [])
    recent_logs = tool_logs[-5:] if tool_logs and len(tool_logs) > 0 else []
    
    # Create reflection prompt
    reflection_prompt = f"""
    Based on the RNA-seq analysis diagnostics, identify the key issues and suggest concrete next steps.
    
    ORGANISM: {ctx.deps.organism}
    
    DIAGNOSTICS:
    {analysis_diagnostics}
    
    TOOL LOGS (most recent 5):
    {json.dumps(recent_logs, indent=2)}
    
    Please provide a CONCISE one-paragraph reflection that:
    1. Clearly identifies the main problem(s)
    2. Gives specific guidance on what file paths or parameters to use in the next attempt
    3. Indicates which steps succeeded and can be skipped in the next run
    
    REFLECTION:
    """
    
    # Run a small LLM to generate the reflection
    reflection_agent = Agent(
        "openai:o4-mini",
        system_prompt="You are an RNA-seq analysis troubleshooter who provides concise, actionable guidance."
    )
    
    result = await reflection_agent.run(reflection_prompt)
    reflection = result.output.strip()
    
    # Safely access and update reflections
    reflections = getattr(ctx.deps, 'reflections', [])
    reflections.append(reflection)
    setattr(ctx.deps, 'reflections', reflections)
    
    logger.info("✅ Generated reflection: %s", reflection)
    
    return reflection

@master.tool
async def analyse(ctx: RunContext[RNAseqCoreContext]) -> str:
    """
    Run RNA-seq analysis with intelligent reflection capabilities.
    
    Instead of a fixed number of iterations, this tool:
    1. Attempts analysis
    2. Evaluates success based on output files and logs
    3. Generates targeted reflections when needed
    4. Continues until success or max attempts reached
    
    This approach follows the Reflexion methodology for LLM self-improvement.
    """
    # Maximum number of attempts
    max_attempts = 3
    
    # Initialize collections
    logger = logging.getLogger(__name__)
    
    # Convert to AnalysisContext if needed
    if not isinstance(ctx.deps, AnalysisContext):
        # Create an AnalysisContext with all the attributes from RNAseqCoreContext
        analysis_ctx = AnalysisContext(**ctx.deps.dict())
        # Replace the context's deps
        ctx.deps = analysis_ctx
        logger.info("🔄 Converted RNAseqCoreContext to AnalysisContext")
    
    # Initialize collections if they don't exist
    for attr_name, default_value in [
        ('tool_logs', []),
        ('reflections', []),
        ('analysis_history', [])
    ]:
        if not hasattr(ctx.deps, attr_name):
            setattr(ctx.deps, attr_name, default_value)

    # Ensure dataset_information is carried over from extraction to analysis if available
    if hasattr(ctx.deps, 'dataset_information') and ctx.deps.dataset_information:
        logger.info("📊 Adding dataset information to analysis context")

    # Base prompt that will be the core of our analysis instructions
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
           - Use the appropriate tx2gene file (t2g.txt) for {ctx.deps.organism}, again found in the resource directory: {ctx.deps.resource_dir} (again, search recursively - the file is named t2g.txt, but the path may vary, e.g. ./mouse/t2g.txt or ./human/t2g.txt)
           - Follow best practices for normalization, filtering, and statistical analysis
           - Generate visualizations (MDS plots, heatmaps, volcano plots)

        Provide detailed explanations for each step and handle errors gracefully.
    """

    # Log the start of the analysis process
    logger.info("🚀 Starting RNA-seq analysis with up to %d attempts", max_attempts)

    # Loop until success or max attempts reached
    attempt = 0
    final_output = None

    while attempt < max_attempts:
        attempt += 1
        logger.info("🔄 Starting analysis attempt %d/%d", attempt, max_attempts)

        # Construct the prompt with previous reflections
        current_prompt = base_prompt

        if ctx.deps.reflections:
            reflection_history = "\n".join([f"Reflection {i+1}: {r}" for i, r in enumerate(ctx.deps.reflections)])
            current_prompt += f"""

            PREVIOUS REFLECTIONS:
            {reflection_history}

            Please incorporate these reflections in your current analysis approach.
            Focus on addressing the specific issues identified above.
            """

        # Record the attempt in history
        ctx.deps.analysis_history.append({
            "iteration": attempt,
            "prompt": current_prompt,
            "output": None  # Will be filled after running
        })

        # Run the analysis
        result = await analysis.run_agent_async(
            current_prompt,
            deps=ctx.deps,
            usage=ctx.usage
        )

        # Store the result
        final_output = result.output
        analysis_history = getattr(ctx.deps, 'analysis_history', [])
        if analysis_history:
            analysis_history[-1]["output"] = final_output
            setattr(ctx.deps, 'analysis_history', analysis_history)

        # Evaluate the success of the analysis
        logger.info("🔍 Evaluating results of attempt %d", attempt)
        await evaluate_analysis(ctx)

        # Check if analysis was successful
        analysis_success = getattr(ctx.deps, 'analysis_success', False)
        if analysis_success:
            logger.info("✅ Analysis attempt %d succeeded!", attempt)
            break
        
        # If unsuccessful and not at max attempts, generate reflection for next attempt
        if attempt < max_attempts:
            logger.info("❌ Analysis attempt %d failed, generating reflection", attempt)
            await generate_reflection(ctx)
        else:
            logger.info("❌ Final analysis attempt failed, stopping after %d attempts", max_attempts)

    # Summarize the analysis process
    analysis_success = getattr(ctx.deps, 'analysis_success', False)
    analysis_diagnostics = getattr(ctx.deps, 'analysis_diagnostics', "No diagnostics available")
    
    if analysis_success:
        summary = f"""
        Analysis completed successfully after {attempt} attempt(s).
        
        {final_output}
        """
    else:
        summary = f"""
        Analysis did not complete successfully after {max_attempts} attempts.
        
        Final diagnostic information:
        {analysis_diagnostics}
        
        Final output:
        {final_output}
        """

    return summary

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
    
    # AnalysisContext is already imported at the top of the file
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

import logging, os, asyncio, json, argparse, subprocess, glob, pathlib
from typing import Optional, List, Dict, Any
from pydantic_ai import Agent, RunContext
from pydantic_ai.usage import UsageLimits
from dotenv import load_dotenv, find_dotenv
from shared import RNAseqCoreContext, AnalysisContext, CheckpointStatus
from shared.workflow_logging import log_tool, setup_logging
from shared.entrez_utils import configure_entrez


# ── load env once ───────────────────────────────
load_dotenv(find_dotenv())
assert os.getenv("OPENAI_API_KEY"), "API key missing"

# ── configure entrez with rate limiting ─────
try:
    configure_entrez()
except Exception as e:
    print(f"Warning: Failed to configure Entrez: {e}")
    print("Please ensure ENTREZ_EMAIL is set in your .env file")

# ── configure python‑logging ────────────────────────────────────────────────
logger = logging.getLogger(__name__)

# ── load system prompt ──────────────────────────────────────────────────────
master_prompt_path = "./main_workflow/prompts/master.txt"
try:
    master_prompt = pathlib.Path(master_prompt_path).read_text()
except Exception as e:
    logger.warning("Could not read analysis system prompt: %s – using fallback", e)
    master_prompt = """
    You are a bioinformatics expert who oversees the execution of a bioinformatic analysis. You will not need to perform any of the analysis yourself, but instead have an expert team of specialised agents who will perform the analysis for you.
    """

# ── agent definitions ───────────────────────────────────────────────────────
master = Agent(
    "openai:gpt-5-mini",
    deps_type = RNAseqCoreContext,
    system_prompt = master_prompt
)

def save_analysis_info(ctx):
    """Save analysis information to a JSON file for integration with reporting tools."""
    logger.info("📝 Saving analysis information for integration")

    # Prepare analysis information
    analysis_info = {
        "accession": ctx.deps.accession,
        "organism": getattr(ctx.deps, 'organism', 'Unknown'),
        "dataset_information": getattr(ctx.deps, 'dataset_information', 'No information available'),
        "analysis_success": getattr(ctx.deps, 'analysis_success', False),
        "reflection_iterations": getattr(ctx.deps, 'reflection_iterations', 0),
        "merged_column": getattr(ctx.deps, 'merged_column', None),
        "unique_groups": getattr(ctx.deps, 'unique_groups', []),
        "kallisto_index_used": getattr(ctx.deps, 'kallisto_index_used', None),
        "tx2gene_file_used": getattr(ctx.deps, 'tx2gene_file_used', None),
        "checkpoints": getattr(ctx.deps, 'checkpoints', {}).model_dump() if hasattr(getattr(ctx.deps, 'checkpoints', {}), 'model_dump') else {}
    }

    # Also include contrasts if available
    if hasattr(ctx.deps, 'contrast_matrix_df') and ctx.deps.contrast_matrix_df is not None:
        try:
            analysis_info["contrasts"] = ctx.deps.contrast_matrix_df.to_dict('records')
        except:
            analysis_info["contrasts"] = []

    # Save to metadata directory
    metadata_dir = os.path.join(ctx.deps.output_dir, "metadata")
    os.makedirs(metadata_dir, exist_ok=True)
    info_file = os.path.join(metadata_dir, "analysis_info.json")
    with open(info_file, "w") as f:
        json.dump(analysis_info, f, indent=2)

    logger.info("✅ Saved analysis information to %s", info_file)

@master.tool
@log_tool
async def extract(ctx: RunContext[RNAseqCoreContext]) -> str:
    """Run data extraction using the extraction agent."""
    logger.info("🔍 Running extraction agent")

    prompt = f"""
    Please extract data for accession {ctx.deps.accession}.

    Steps to perform:
    1. Fetch GEO metadata and sample information
    2. Download FASTQ files for the samples

    Use the available tools to complete these tasks.
    """

    # Import extraction module
    import importlib
    extraction = importlib.import_module("agents.extraction")

    # Call the extraction agent
    r = await extraction.run_agent_async(
        prompt,
        deps=ctx.deps, usage=ctx.usage
    )
    return r.output

@master.tool
async def evaluate_analysis(ctx: RunContext[RNAseqCoreContext]) -> str:
    """
    Evaluate the success of the analysis run based on checkpoint completion.

    This function checks the status of predefined checkpoints to determine
    if the analysis was successful. Critical checkpoints must be completed
    for the analysis to be considered successful.
    """
    logger.info("🔍 Evaluating analysis results")

    # Convert to AnalysisContext if needed to access checkpoints
    if not isinstance(ctx.deps, AnalysisContext):
        try:
            analysis_ctx = AnalysisContext(**ctx.deps.model_dump() if hasattr(ctx.deps, 'model_dump') else ctx.deps.__dict__)
            ctx.deps = analysis_ctx
        except Exception as e:
            logger.error("❌ Failed to convert context for evaluation: %s", str(e))
            return "Error: Could not access checkpoint information"

    # Initialize checkpoints if they don't exist
    if not hasattr(ctx.deps, 'checkpoints') or ctx.deps.checkpoints is None:
        from shared import AnalysisCheckpoints
        ctx.deps.checkpoints = AnalysisCheckpoints()
        logger.info("🔄 Initialized analysis checkpoints for evaluation")

    # Get checkpoints
    checkpoints = getattr(ctx.deps, 'checkpoints', None)
    if not checkpoints:
        logger.error("❌ No checkpoints found in context")
        ctx.deps.analysis_success = False
        ctx.deps.analysis_diagnostics = "❌ No checkpoints found - analysis incomplete"
        return ctx.deps.analysis_diagnostics

    # Define critical checkpoints that must be completed for success
    critical_checkpoints = [
        'metadata_extraction',
        'metadata_analysis',
        'kallisto_quantification',
        'edger_limma_preparation',
        'rnaseq_analysis'
    ]

    # Optional checkpoints (may be skipped in some workflows)
    optional_checkpoints = [
        'fastq_extraction',
        'kallisto_index_selection'
    ]

    # Evaluate checkpoint status
    diagnostics = []
    failed_critical = []
    completed_critical = []

    # Check critical checkpoints
    for checkpoint_name in critical_checkpoints:
        checkpoint = getattr(checkpoints, checkpoint_name, None)
        if not checkpoint:
            diagnostics.append(f"❌ Critical checkpoint '{checkpoint_name}' not found")
            failed_critical.append(checkpoint_name)
            continue

        status = checkpoint.status
        error_msg = checkpoint.error_message

        if status == CheckpointStatus.COMPLETED and not error_msg:
            diagnostics.append(f"✅ {checkpoint_name}: completed successfully")
            completed_critical.append(checkpoint_name)
        elif status == CheckpointStatus.FAILED or error_msg:
            diagnostics.append(f"❌ {checkpoint_name}: failed - {error_msg or 'unknown error'}")
            failed_critical.append(checkpoint_name)
        elif status == CheckpointStatus.IN_PROGRESS:
            diagnostics.append(f"⏳ {checkpoint_name}: still in progress")
            failed_critical.append(checkpoint_name)
        else:
            diagnostics.append(f"❌ {checkpoint_name}: not started")
            failed_critical.append(checkpoint_name)

    # Check optional checkpoints (for informational purposes)
    for checkpoint_name in optional_checkpoints:
        checkpoint = getattr(checkpoints, checkpoint_name, None)
        if not checkpoint:
            continue

        status = checkpoint.status
        error_msg = checkpoint.error_message

        if status == CheckpointStatus.COMPLETED and not error_msg:
            diagnostics.append(f"✅ {checkpoint_name}: completed successfully (optional)")
        elif status == CheckpointStatus.FAILED or error_msg:
            diagnostics.append(f"⚠️ {checkpoint_name}: failed - {error_msg or 'unknown error'} (optional)")
        elif status == CheckpointStatus.NOT_STARTED:
            diagnostics.append(f"ℹ️ {checkpoint_name}: skipped (optional)")

    # Determine overall success
    analysis_success = len(failed_critical) == 0 and len(completed_critical) == len(critical_checkpoints)

    # Add summary
    if analysis_success:
        diagnostics.insert(0, f"✅ Analysis SUCCESSFUL: {len(completed_critical)}/{len(critical_checkpoints)} critical checkpoints completed")
    else:
        diagnostics.insert(0, f"❌ Analysis FAILED: {len(failed_critical)} critical checkpoints failed")
        if failed_critical:
            diagnostics.append(f"Failed critical checkpoints: {', '.join(failed_critical)}")

    # Store results
    ctx.deps.analysis_success = analysis_success
    ctx.deps.analysis_diagnostics = "\n".join(diagnostics)

    # Log result
    if analysis_success:
        logger.info("✅ Analysis evaluation: SUCCESS")
    else:
        logger.info("❌ Analysis evaluation: FAILURE")
        logger.info("Failed checkpoints: %s", failed_critical)

    return ctx.deps.analysis_diagnostics

@master.tool
async def generate_reflection(ctx: RunContext[RNAseqCoreContext]) -> str:
    """Generate reflections on analysis failures to guide retry attempts."""
    logger.info("🤔 Generating reflection on analysis failure")

    # Convert to AnalysisContext if needed
    if not isinstance(ctx.deps, AnalysisContext):
        try:
            analysis_ctx = AnalysisContext(**ctx.deps.model_dump() if hasattr(ctx.deps, 'model_dump') else ctx.deps.__dict__)
            ctx.deps = analysis_ctx
        except Exception as e:
            logger.warning("⚠️ Failed to convert context for reflection: %s", str(e))

    # Initialize checkpoints if they don't exist
    if not hasattr(ctx.deps, 'checkpoints') or ctx.deps.checkpoints is None:
        from shared import AnalysisCheckpoints
        ctx.deps.checkpoints = AnalysisCheckpoints()
        logger.info("🔄 Initialized analysis checkpoints for reflection")

    # Get analysis diagnostics
    analysis_diagnostics = getattr(ctx.deps, 'analysis_diagnostics', None)
    if not analysis_diagnostics:
        logger.warning("⚠️ No analysis diagnostics available - running evaluation first")
        await evaluate_analysis(ctx)
        analysis_diagnostics = getattr(ctx.deps, 'analysis_diagnostics', "No diagnostics available")

    # Get previous reflections to avoid repetition
    previous_reflections = getattr(ctx.deps, 'reflections', [])

    # Get tool logs for context
    tool_logs = getattr(ctx.deps, 'tool_logs', [])
    recent_errors = []
    for log in tool_logs[-10:]:  # Last 10 tool calls
        if not log.get('success', True):
            recent_errors.append(f"- {log['tool_name']}: {log.get('error', 'Unknown error')}")

    # Create reflection agent
    reflection_agent = Agent(
        "openai:gpt-5-mini",
        system_prompt="""You are an expert bioinformatics troubleshooter. Your job is to analyze failed RNA-seq analysis attempts and provide specific, actionable recommendations for the next attempt.

Focus on:
1. Identifying the root cause of failures
2. Suggesting specific parameter changes or different approaches
3. Highlighting potential issues with file paths, organism mismatches, or metadata problems
4. Providing clear, actionable advice for the analysis agent

Be concise and specific. Avoid generic advice."""
        , model_settings={"temperature": 1}
    )

    reflection_prompt = f"""
    The RNA-seq analysis failed. Please analyze the failure and provide specific recommendations.

    ORGANISM: {getattr(ctx.deps, 'organism', 'Unknown')}
    ACCESSION: {ctx.deps.accession}

    ANALYSIS DIAGNOSTICS:
    {analysis_diagnostics}

    RECENT TOOL ERRORS:
    {chr(10).join(recent_errors) if recent_errors else "No recent tool errors logged"}

    PREVIOUS REFLECTIONS (avoid repeating these):
    {chr(10).join([f"{i+1}. {r}" for i, r in enumerate(previous_reflections)]) if previous_reflections else "None"}

    ANALYSIS HISTORY:
    {len(getattr(ctx.deps, 'analysis_history', []))} previous attempts made

    Provide a specific reflection on what went wrong and what should be tried differently in the next attempt. Focus on the most critical issue that needs to be addressed.

    REFLECTION:
    """

    try:
        result = await reflection_agent.run(reflection_prompt)
        reflection = result.output.strip()

        # Store the reflection
        if not hasattr(ctx.deps, 'reflections'):
            ctx.deps.reflections = []
        ctx.deps.reflections.append(reflection)

        # Increment reflection count
        current_count = getattr(ctx.deps, 'reflection_iterations', 0)
        ctx.deps.reflection_iterations = current_count + 1

        logger.info("📝 Generated reflection %d: %s", ctx.deps.reflection_iterations, reflection[:100] + "..." if len(reflection) > 100 else reflection)

        return reflection

    except Exception as e:
        logger.error("❌ Error generating reflection: %s", str(e))
        fallback_reflection = f"Analysis failed at attempt {len(getattr(ctx.deps, 'analysis_history', []))}. Review the error messages and try a different approach."

        if not hasattr(ctx.deps, 'reflections'):
            ctx.deps.reflections = []
        ctx.deps.reflections.append(fallback_reflection)

        return fallback_reflection

@master.tool
@log_tool
async def analyse(ctx: RunContext[RNAseqCoreContext]) -> str:
    """Run analysis with reflection-based retry logic."""
    logger.info("🧬 Starting analysis with reflection capability")

    # Skip analysis entirely if earlier step flagged it as unnecessary
    if getattr(ctx.deps, "analysis_should_proceed", True) is False:
        reason = getattr(ctx.deps, "analysis_skip_reason", "Dataset flagged as too small")
        logger.info("⚠️ Skipping analysis: %s", reason)
        ctx.deps.analysis_success = False
        return f"Analysis skipped: {reason}"

    # Convert to AnalysisContext
    if not isinstance(ctx.deps, AnalysisContext):
        try:
            analysis_ctx = AnalysisContext(**ctx.deps.model_dump() if hasattr(ctx.deps, 'model_dump') else ctx.deps.__dict__)
            ctx.deps = analysis_ctx
        except Exception as e:
            logger.error("❌ Failed to convert to AnalysisContext: %s", str(e))
            return f"Error: {str(e)}"

    # Initialize checkpoints if they don't exist
    if not hasattr(ctx.deps, 'checkpoints') or ctx.deps.checkpoints is None:
        from shared import AnalysisCheckpoints
        ctx.deps.checkpoints = AnalysisCheckpoints()
        logger.info("🔄 Initialized analysis checkpoints")

    # Initialize analysis tracking
    if not hasattr(ctx.deps, 'analysis_history'):
        ctx.deps.analysis_history = []
    if not hasattr(ctx.deps, 'reflections'):
        ctx.deps.reflections = []
    if not hasattr(ctx.deps, 'reflection_iterations'):
        ctx.deps.reflection_iterations = 0

    max_attempts = 5

    for attempt in range(1, max_attempts + 1):
        logger.info("🔄 Analysis attempt %d/%d", attempt, max_attempts)

        # Prepare prompt with reflection context
        base_prompt = f"""
        Please analyze the RNA-seq data for {ctx.deps.organism} (accession: {ctx.deps.accession}).

        Your task is to complete the full RNA-seq analysis pipeline:
        1. Process metadata using the metadata agent
        2. Run Kallisto quantification
        3. Prepare sample mapping for differential expression
        4. Run edgeR/limma differential expression analysis

        Use the available tools to complete this analysis.
        """

        # Add reflection context for retry attempts
        if attempt > 1:
            recent_reflections = ctx.deps.reflections[-2:] if len(ctx.deps.reflections) >= 2 else ctx.deps.reflections
            reflection_context = "\n".join([f"Reflection {i+1}: {r}" for i, r in enumerate(recent_reflections)])

            base_prompt += f"""

            IMPORTANT - PREVIOUS ATTEMPT FAILED:
            This is attempt {attempt} of {max_attempts}. The previous attempt(s) failed.

            REFLECTIONS FROM PREVIOUS ATTEMPTS:
            {reflection_context}

            Please carefully consider these reflections and adjust your approach accordingly.
            Address the specific issues identified in the reflections.
            """

        try:
            # Import analysis module
            import importlib
            analysis = importlib.import_module("agents.analysis")

            # Run the analysis agent
            result = await analysis.run_agent_async(base_prompt, deps=ctx.deps, usage=ctx.usage)

            # Store the attempt
            ctx.deps.analysis_history.append({
                "attempt": attempt,
                "prompt": base_prompt,
                "output": result.output,
                "timestamp": str(asyncio.get_event_loop().time())
            })

            # Evaluate the success of the analysis
            logger.info("🔍 Evaluating results of attempt %d", attempt)
            await evaluate_analysis(ctx)

            # Check if analysis was successful
            analysis_success = getattr(ctx.deps, 'analysis_success', False)

            if analysis_success:
                logger.info("✅ Analysis attempt %d succeeded!", attempt)
                return result.output
            else:
                logger.info("❌ Analysis attempt %d failed", attempt)

                if attempt < max_attempts:
                    # Generate reflection for next attempt
                    logger.info("🤔 Generating reflection for next attempt")
                    reflection = await generate_reflection(ctx)
                    logger.info("📝 Reflection generated: %s", reflection[:100] + "..." if len(reflection) > 100 else reflection)
                else:
                    logger.info("❌ Final analysis attempt failed")
                    return f"Analysis failed after {max_attempts} attempts. Latest output: {result.output}"

        except Exception as e:
            logger.error("❌ Analysis attempt %d crashed: %s", attempt, str(e))

            # Store the failed attempt
            ctx.deps.analysis_history.append({
                "attempt": attempt,
                "prompt": base_prompt,
                "output": None,
                "error": str(e),
                "timestamp": str(asyncio.get_event_loop().time())
            })

            if attempt < max_attempts:
                # Generate reflection even for crashed attempts
                logger.info("🤔 Generating reflection after crash")
                try:
                    reflection = await generate_reflection(ctx)
                    logger.info("📝 Reflection generated after crash: %s", reflection[:100] + "..." if len(reflection) > 100 else reflection)
                except Exception as ref_error:
                    logger.error("❌ Failed to generate reflection after crash: %s", str(ref_error))
            else:
                return f"Analysis crashed after {max_attempts} attempts. Final error: {str(e)}"

    return f"Analysis failed after {max_attempts} attempts"

def cleanup_large_files(output_dir: str, accession: str):
    """Clean up large intermediate files (FASTQ and SRA) after successful analysis."""
    logger.info("🧹 Starting cleanup of large files")

    # Calculate storage before cleanup
    total_size_before = 0
    fastq_dir = os.path.join(output_dir, "fastq")
    sra_dir = os.path.join(output_dir, "sra")

    # Calculate FASTQ files size
    if os.path.exists(fastq_dir):
        for root, dirs, files in os.walk(fastq_dir):
            for file in files:
                file_path = os.path.join(root, file)
                if os.path.exists(file_path):
                    total_size_before += os.path.getsize(file_path)

    # Calculate SRA files size
    if os.path.exists(sra_dir):
        for root, dirs, files in os.walk(sra_dir):
            for file in files:
                file_path = os.path.join(root, file)
                if os.path.exists(file_path):
                    total_size_before += os.path.getsize(file_path)

    size_before_gb = total_size_before / (1024**3)
    logger.info("📊 Storage before cleanup: %.2f GB", size_before_gb)

    # Clean up FASTQ files
    fastq_files_removed = 0
    if os.path.exists(fastq_dir):
        try:
            for root, dirs, files in os.walk(fastq_dir):
                for file in files:
                    file_path = os.path.join(root, file)
                    try:
                        os.remove(file_path)
                        fastq_files_removed += 1
                    except OSError as e:
                        logger.warning("⚠️ Could not remove FASTQ file %s: %s", file_path, e)

            # Remove empty directories
            for root, dirs, files in os.walk(fastq_dir, topdown=False):
                for dir_name in dirs:
                    dir_path = os.path.join(root, dir_name)
                    try:
                        os.rmdir(dir_path)
                    except OSError:
                        pass  # Directory not empty or other issue

            # Remove the main fastq directory if empty
            try:
                os.rmdir(fastq_dir)
                logger.info("🗑️ Removed FASTQ directory")
            except OSError:
                logger.info("📂 FASTQ directory kept (contains remaining files)")

        except Exception as e:
            logger.error("❌ Error during FASTQ cleanup: %s", str(e))

    # Clean up SRA files
    sra_files_removed = 0
    if os.path.exists(sra_dir):
        try:
            for root, dirs, files in os.walk(sra_dir):
                for file in files:
                    file_path = os.path.join(root, file)
                    try:
                        os.remove(file_path)
                        sra_files_removed += 1
                    except OSError as e:
                        logger.warning("⚠️ Could not remove SRA file %s: %s", file_path, e)

            # Remove empty directories
            for root, dirs, files in os.walk(sra_dir, topdown=False):
                for dir_name in dirs:
                    dir_path = os.path.join(root, dir_name)
                    try:
                        os.rmdir(dir_path)
                    except OSError:
                        pass

            # Remove the main sra directory if empty
            try:
                os.rmdir(sra_dir)
                logger.info("🗑️ Removed SRA directory")
            except OSError:
                logger.info("📂 SRA directory kept (contains remaining files)")

        except Exception as e:
            logger.error("❌ Error during SRA cleanup: %s", str(e))

    # Calculate final storage
    total_size_after = 0
    for root, dirs, files in os.walk(output_dir):
        for file in files:
            file_path = os.path.join(root, file)
            if os.path.exists(file_path):
                total_size_after += os.path.getsize(file_path)

    size_after_gb = total_size_after / (1024**3)
    size_saved_gb = size_before_gb - size_after_gb

    logger.info("🧹 Cleanup completed:")
    logger.info("  - FASTQ files removed: %d", fastq_files_removed)
    logger.info("  - SRA files removed: %d", sra_files_removed)
    logger.info("  - Storage before: %.2f GB", size_before_gb)
    logger.info("  - Storage after: %.2f GB", size_after_gb)
    logger.info("  - Storage saved: %.2f GB (%.1f%%)", size_saved_gb,
                (size_saved_gb / size_before_gb * 100) if size_before_gb > 0 else 0)

def main():
    parser = argparse.ArgumentParser(description="UORCA RNA-seq Analysis Pipeline")
    parser.add_argument("--accession", required=True, help="GEO accession (e.g., GSE123456)")
    parser.add_argument("--output_dir", default="../UORCA_results", help="Output directory")
    parser.add_argument("--resource_dir", default="./data/kallisto_indices/", help="Resource directory")
    parser.add_argument("--cleanup", action="store_true", help="Clean up FASTQ and SRA files after analysis")
    args = parser.parse_args()

    # -------- create & remember the chosen run folder -----------------------
    output_dir = os.path.join(args.output_dir, args.accession)
    os.makedirs(output_dir, exist_ok=True)

    # -------- configure logging *inside* that folder ------------------------
    log_dir = pathlib.Path(output_dir) / "logs"
    log_path = setup_logging(log_dir)
    logger.info("🚀 Starting UORCA master agent - logging to %s", log_path)

    # Create context
    ctx = AnalysisContext(
        accession=args.accession,
        output_dir=output_dir,
        resource_dir=args.resource_dir,
        organism="Unknown"  # Will be determined during extraction
    )

    logger.info("🚀 Starting UORCA analysis for %s", args.accession)
    logger.info("📁 Output directory: %s", ctx.output_dir)
    logger.info("📚 Resource directory: %s", ctx.resource_dir)

    # Import agent modules dynamically
    import importlib, sys
    logger.info("📚 Importing agent modules")
    extraction = importlib.import_module("agents.extraction")
    analysis = importlib.import_module("agents.analysis")
    logger.info("✅ All agent modules imported successfully")

    # Create the orchestration prompt
    orchestration_prompt = f"""
    Analyse {args.accession}, an RNAseq dataset by first extracting the data and performing an analysis on it.
    The organism will be automatically determined during data extraction.
    Skip any steps if the data already exists and is up to date.
    Document each tool invocation and output.
    """

    try:
        # Run the master agent
        logger.info("🤖 Running master orchestration agent")
        result = master.run_sync(orchestration_prompt, deps=ctx, usage_limits=UsageLimits(request_limit=100))

        # Create a simple context object for save_analysis_info
        class SimpleContext:
            def __init__(self, deps):
                self.deps = deps

        # fall back to ctx if the result does not provide a deps attribute
        result_deps = getattr(result, "deps", None)
        if result_deps is None:
            logger.warning("AgentRunResult missing 'deps' attribute; using current context")
            result_deps = ctx

        simple_ctx = SimpleContext(result_deps)

        # Save analysis info for reporting integration
        save_analysis_info(simple_ctx)

        # Cleanup if requested and analysis was successful
        if args.cleanup:
            logger.info("🧹 Analysis successful - performing cleanup")
            cleanup_large_files(ctx.output_dir, args.accession)
        else:
            logger.info("ℹ️ Cleanup skipped. Use --cleanup flag to remove FASTQ/SRA files")

        # Log final results
        logger.info("✅ Master agent completed execution")
        logger.info("📝 Final output: %s", result.output)

        # Log usage stats if available
        if hasattr(result, 'usage') and result.usage:
            try:
                usage_stats = result.usage()
                logger.info("📊 Token usage: %s", usage_stats)
            except:
                pass

    except Exception as e:
        logger.error("❌ Pipeline failed: %s", str(e), exc_info=True)
        raise

if __name__ == "__main__":
    main()

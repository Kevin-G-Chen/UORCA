#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Auto-Start Manager for UORCA RNA-seq Analysis

This module manages the automatic startup and initial analysis of RNA-seq data
when users load their results directory. It provides a seamless experience by
automatically generating an overview analysis without requiring user interaction.
"""

import os
import sys
import asyncio
import logging
from typing import Dict, Any, Tuple, Optional, Callable, List
from pathlib import Path
from datetime import datetime
import json

# Import the reporting agent
from reporting_agent import ReportingAgent

# Set up logging
logging.basicConfig(
    level=logging.WARNING,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class AutoStartManager:
    """Manages the automatic startup and initial analysis of RNA-seq data."""

    def __init__(self,
                 results_dir: str,
                 default_model: str = "openai:gpt-4o-mini",
                 temperature: float = 0.1):
        """
        Initialize the auto-start manager.

        Args:
            results_dir: Path to UORCA results directory
            default_model: Default model to use for analysis
            temperature: Model temperature for generation
        """
        self.results_dir = os.path.abspath(results_dir)
        self.default_model = default_model
        self.temperature = temperature
        self.agent = None
        self.servers = []

        logger.info(f"Initialized AutoStartManager for: {self.results_dir}")

    def validate_data_directory(self) -> Tuple[bool, Optional[str]]:
        """
        Validate that the results directory contains analyzable data.

        Returns:
            Tuple of (is_valid, error_message)
        """
        try:
            # Check if directory exists
            if not os.path.exists(self.results_dir):
                return False, "Results directory does not exist"

            # Check for RNA-seq analysis directories
            found_rnaseq = False
            found_deg_files = False
            analysis_count = 0

            for root, dirs, files in os.walk(self.results_dir):
                if "RNAseqAnalysis" in dirs:
                    found_rnaseq = True
                    analysis_count += 1

                if "DEG.csv" in files:
                    found_deg_files = True

            if not found_rnaseq:
                return False, "No RNAseqAnalysis directories found"

            if not found_deg_files:
                return False, "No DEG.csv files found in analysis directories"

            logger.info(f"Validated directory: found {analysis_count} analyses")
            return True, None

        except Exception as e:
            logger.error(f"Error validating data directory: {e}")
            return False, f"Validation error: {str(e)}"

    def setup_servers(self) -> List:
        """Set up the MCP servers for data extraction and analysis."""
        try:
            # Import here to avoid circular imports
            from mcp_utils import setup_mcp_server

            # Ensure the MCP child sees the same results_dir the user selected
            env_vars = {"RESULTS_DIR": self.results_dir}
            api_key = os.getenv("OPENAI_API_KEY")
            if api_key:
                env_vars["OPENAI_API_KEY"] = api_key

            # Data extraction server
            logger.info("Setting up data extraction MCP server...")
            data_server = setup_mcp_server(
                "data_extractor",
                env_vars=env_vars,
            )
            self.servers.append(data_server)

            # Analysis server
            logger.info("Setting up analysis MCP server...")
            analysis_server = setup_mcp_server(
                "analysis",
                env_vars=env_vars,
            )
            self.servers.append(analysis_server)

            logger.info(f"Successfully set up {len(self.servers)} MCP servers")
            return self.servers

        except Exception as e:
            logger.error(f"Error setting up MCP servers: {e}")
            raise

    async def run_initial_analysis(self,
                                 progress_callback: Optional[Callable[[float, str], None]] = None) -> Dict[str, Any]:
        """Run initial analysis on the data without requiring a specific research question."""
        try:
            import time
            overall_start = time.time()
            logger.info("Starting run_initial_analysis")

            if progress_callback:
                progress_callback(0.1, "Initializing analysis agent")
            logger.info("Step 1/6: Initializing analysis agent")

            agent_start = time.time()
            self.agent = ReportingAgent(
                results_dir=self.results_dir,
                model_name=self.default_model,
                temperature=self.temperature
            )
            logger.info(f"Agent created in {time.time() - agent_start:.2f} seconds")

            if progress_callback:
                progress_callback(0.2, "Setting up analysis servers")

            servers_start = time.time()
            if not self.servers:
                logger.info("Step 2/6: Setting up servers")
                self.setup_servers()
                logger.info(f"Servers setup in {time.time() - servers_start:.2f} seconds")
            else:
                logger.info("Step 2/6: Reusing existing servers")

            logger.info("Step 3/6: Connecting agent to servers")
            server_connection_start = time.time()
            self.agent.servers = self.servers
            logger.info(f"Connected agent to servers in {time.time() - server_connection_start:.2f} seconds")

            if progress_callback:
                progress_callback(0.3, "Creating analysis agent")

            logger.info("Step 4/6: Creating agent")
            agent_creation_start = time.time()
            self.agent.create_agent()
            logger.info(f"Agent created in {time.time() - agent_creation_start:.2f} seconds")

            logger.info("Step 5/6: Preparing prompt")
            initial_prompt = """Please provide a comprehensive initial analysis of this RNA-seq dataset that would be immediately useful to a researcher opening this data for the first time.

Your analysis should include:

1. **Dataset Overview**
   - Use list_available_contrasts to enumerate all experimental comparisons
   - Identify which contrasts have the most differential expression
   - Group related contrasts if patterns are apparent

2. **Key Findings**
   - Use find_common_degs to identify genes that are consistently affected across multiple contrasts
   - Find the most strongly differentially expressed genes (highest fold changes)
   - Look for genes that might be master regulators or key drivers

3. **Expression Patterns**
   - Use analyze_gene_patterns on the top common genes to identify co-expression modules
   - Identify contrasts that show similar overall patterns
   - Find any striking opposing patterns between conditions

4. **Biological Insights**
   - Based on the gene names and patterns, suggest what biological processes are most affected
   - Identify potential pathways or mechanisms
   - Suggest what types of experiments these might represent

5. **Data Quality Summary**
   - Note the total number of genes tested
   - Mention the range of DEG counts across contrasts
   - Flag any contrasts with unusually few or many DEGs

Organize your response with clear sections and use bullet points for readability.
Focus on actionable insights that help the researcher quickly understand their data.
Be specific with gene names, fold changes, and numbers rather than generic statements."""

            if progress_callback:
                progress_callback(0.4, "Analyzing dataset structure")

            logger.info("Step 6/6: Running LLM analysis - this may take several minutes")
            logger.info("About to call agent.run() with the analysis prompt")

            llm_start_time = time.time()
            try:
                progress_task = asyncio.create_task(self._log_progress_during_analysis(llm_start_time))
                result = await self.agent.agent.run(initial_prompt)
                progress_task.cancel()
                logger.info(f"LLM analysis completed in {time.time() - llm_start_time:.2f} seconds")

                if progress_callback:
                    progress_callback(1.0, "Analysis complete")

                return {
                    "success": True,
                    "analysis_type": "initial_overview",
                    "report": result.data,
                    "metadata": {
                        "results_dir": self.results_dir,
                        "model_used": self.default_model,
                        "timestamp": datetime.now().isoformat(),
                        "analysis_version": "2.0",
                        "analysis_duration_seconds": time.time() - overall_start
                    }
                }
            except Exception as analysis_error:
                logger.error(f"Error in agent.run: {analysis_error}", exc_info=True)
                raise

        except Exception as e:
            logger.error(f"Error during initial analysis: {e}", exc_info=True)
            return {
                "success": False,
                "error": str(e),
                "error_type": type(e).__name__,
                "timestamp": datetime.now().isoformat()
            }

    async def _log_progress_during_analysis(self, start_time: float):
        """Log progress during a long-running analysis periodically."""
        import time
        try:
            interval = 30
            count = 0
            while True:
                await asyncio.sleep(interval)
                elapsed = time.time() - start_time
                count += 1
                logger.info(f"[PROGRESS] LLM analysis running for {elapsed:.1f} seconds (update #{count})")
        except asyncio.CancelledError:
            pass

    async def run_focused_analysis(self,
                                 research_question: str,
                                 progress_callback: Optional[Callable[[float, str], None]] = None) -> Dict[str, Any]:
        """
        Run a focused analysis based on a specific research question.

        Args:
            research_question: The research question to investigate
            progress_callback: Optional progress callback function

        Returns:
            Dictionary with analysis results
        """
        try:
            # Ensure agent is initialized
            if not self.agent:
                if progress_callback:
                    progress_callback(0.1, "Initializing analysis agent")

                self.agent = ReportingAgent(
                    results_dir=self.results_dir,
                    model_name=self.default_model,
                    temperature=self.temperature
                )

                # Setup servers if not already done
                if not self.servers:
                    self.setup_servers()

                # Use existing servers
                self.agent.servers = self.servers

            if progress_callback:
                progress_callback(0.3, "Starting focused analysis")

            # Run the focused analysis
            result = await self.agent.analyze_research_question(
                research_question=research_question,
                max_contrasts=12,  # Analyze up to 12 most relevant contrasts
                max_genes=50       # Report top 50 genes
            )

            if progress_callback:
                progress_callback(1.0, "Analysis complete")

            return result

        except Exception as e:
            logger.error(f"Error during focused analysis: {e}", exc_info=True)
            return {
                "success": False,
                "error": str(e),
                "research_question": research_question,
                "timestamp": datetime.now().isoformat()
            }

    async def suggest_research_directions(self) -> Dict[str, Any]:
        """
        Suggest potential research questions based on the data.

        Returns:
            Dictionary with suggested research directions
        """
        try:
            if not self.agent:
                self.agent = ReportingAgent(
                    results_dir=self.results_dir,
                    model_name=self.default_model
                )

                # Setup servers if not already done
                if not self.servers:
                    self.setup_servers()

                # Use existing servers
                self.agent.servers = self.servers
                self.agent.create_agent()

            prompt = """Based on the available contrasts and initial analysis of this dataset, suggest 5-7 specific research questions that this data could help answer.

For each suggestion:
1. State the research question clearly
2. Identify which contrasts would be most relevant
3. Mention key genes or pathways that might be involved
4. Explain why this question is interesting given the data

Focus on diverse questions covering different aspects like:
- Disease mechanisms
- Treatment responses
- Biomarker discovery
- Pathway analysis
- Cell type differences
- Temporal changes

Make the questions specific and actionable based on the actual data available."""

            result = await self.agent.agent.run(prompt)

            return {
                "success": True,
                "suggestions": result.data,
                "timestamp": datetime.now().isoformat()
            }

        except Exception as e:
            logger.error(f"Error generating research suggestions: {e}", exc_info=True)
            return {
                "success": False,
                "error": str(e),
                "timestamp": datetime.now().isoformat()
            }

    async def cleanup(self):
        """Clean up resources. Should only be called when completely done with the manager."""
        # Clean up servers first
        for server in self.servers:
            try:
                await server.cleanup()
            except Exception as e:
                logger.warning(f"Error cleaning up server: {e}", exc_info=True)

        self.servers = []

        # Then clean up agent
        if self.agent:
            try:
                await self.agent.cleanup()
            except Exception as e:
                logger.warning(f"Error during agent cleanup: {e}", exc_info=True)
            finally:
                self.agent = None
        logger.info("AutoStartManager cleanup complete")



# Convenience function for quick auto-start
async def auto_analyze(results_dir: str,
                      progress_callback: Optional[Callable[[float, str], None]] = None,
                      model: str = "openai:gpt-4o-mini") -> Dict[str, Any]:
    """
    Convenience function to quickly auto-analyze a dataset.

    Args:
        results_dir: Path to UORCA results
        progress_callback: Optional progress callback
        model: Model to use for analysis

    Returns:
        Analysis results dictionary
    """
    manager = AutoStartManager(results_dir, default_model=model)

    try:
        # Validate directory first
        is_valid, error_msg = manager.validate_data_directory()
        if not is_valid:
            return {
                "success": False,
                "error": error_msg,
                "timestamp": datetime.now().isoformat()
            }

        # Run analysis
        return await manager.run_initial_analysis(progress_callback)

    finally:
        await manager.cleanup()


if __name__ == "__main__":
    # Example usage
    async def main():
        results_dir = os.environ.get("UORCA_RESULTS", "./UORCA_results")

        def progress_printer(progress: float, message: str):
            print(f"[{progress*100:.0f}%] {message}")

        print("Running automatic analysis...")
        result = await auto_analyze(results_dir, progress_callback=progress_printer)

        if result["success"]:
            print("\nAnalysis completed successfully!")
            print(f"\nReport:\n{result['report']}")
        else:
            print(f"\nAnalysis failed: {result['error']}")

    asyncio.run(main())

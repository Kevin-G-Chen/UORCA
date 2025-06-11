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
import time

# Import the reporting agent
from reporting_agent import ReportingAgent
from mcp_utils import verify_server_responsive, setup_mcp_server, is_server_process_alive, restart_server

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
        self.server_env_vars = None  # Store env vars for server restart
        self.last_health_check = 0  # Timestamp of last health check

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
            # Clean up any existing servers first
            if self.servers:
                logger.info("Cleaning up existing servers before setup")
                try:
                    asyncio.run(self._cleanup_servers_sync())
                except Exception as cleanup_error:
                    logger.warning(f"Error cleaning up existing servers: {cleanup_error}")
                self.servers = []

            # Ensure the MCP child sees the same results_dir the user selected
            env_vars = {"RESULTS_DIR": self.results_dir}
            api_key = os.getenv("OPENAI_API_KEY")
            if api_key:
                env_vars["OPENAI_API_KEY"] = api_key
            
            # Store env vars for potential server restart
            self.server_env_vars = env_vars

            # Data extraction server
            logger.info("Setting up data extraction MCP server...")
            data_server = setup_mcp_server(
                "data_extractor",
                env_vars=env_vars,
                max_retries=3,
                retry_delay=1.5
            )
            
            # Verify server is responsive
            if not self._verify_server_health(data_server, "data_extractor"):
                logger.warning("Data extraction server failed initial health check, attempting restart...")
                data_server = restart_server("data_extractor", env_vars)
                if not data_server or not self._verify_server_health(data_server, "data_extractor"):
                    raise RuntimeError("Data extraction server failed health check after restart")
            
            self.servers.append(data_server)

            # Analysis server
            logger.info("Setting up analysis MCP server...")
            analysis_server = setup_mcp_server(
                "analysis",
                env_vars=env_vars,
                max_retries=3,
                retry_delay=1.5
            )
            
            # Verify server is responsive
            if not self._verify_server_health(analysis_server, "analysis"):
                logger.warning("Analysis server failed initial health check, attempting restart...")
                analysis_server = restart_server("analysis", env_vars)
                if not analysis_server or not self._verify_server_health(analysis_server, "analysis"):
                    raise RuntimeError("Analysis server failed health check after restart")
            
            self.servers.append(analysis_server)

            logger.info(f"Successfully set up {len(self.servers)} MCP servers")
            self.last_health_check = time.time()
            return self.servers

        except Exception as e:
            logger.error(f"Error setting up MCP servers: {e}")
            # Clean up any partially created servers
            try:
                asyncio.run(self._cleanup_servers_sync())
            except Exception:
                pass
            self.servers = []
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

    def _verify_server_health(self, server, server_name: str, timeout: float = 3.0) -> bool:
        """Verify that a server is healthy and responsive."""
        try:
            logger.info(f"Verifying health of {server_name} server...")
            
            # Give server a moment to fully initialize
            time.sleep(0.5)
            
            # Quick process check first
            if not is_server_process_alive(server):
                logger.error(f"{server_name} server process terminated")
                return False
            
            # Use the verify function from mcp_utils
            is_healthy = verify_server_responsive(server, timeout)
            
            if is_healthy:
                logger.info(f"{server_name} server health check passed")
            else:
                logger.error(f"{server_name} server health check failed")
                
            return is_healthy
            
        except Exception as e:
            logger.error(f"Error verifying {server_name} server health: {e}")
            return False

    def verify_all_servers_healthy(self, quick_check: bool = False) -> bool:
        """Verify that all servers are healthy and responsive."""
        if not self.servers:
            logger.warning("No servers to verify")
            return False
            
        server_names = ["data_extractor", "analysis"]
        
        # Skip health check if we did one recently (for quick_check mode)
        if quick_check and hasattr(self, 'last_health_check'):
            time_since_check = time.time() - self.last_health_check
            if time_since_check < 30:  # Skip if checked within last 30 seconds
                logger.debug("Skipping health check, done recently")
                return all(is_server_process_alive(server) for server in self.servers)
        
        healthy_count = 0
        for i, server in enumerate(self.servers):
            server_name = server_names[i] if i < len(server_names) else f"server_{i}"
            if self._verify_server_health(server, server_name):
                healthy_count += 1
            else:
                logger.error(f"Server {server_name} failed health check")
                
        all_healthy = healthy_count == len(self.servers)
        if all_healthy:
            logger.info("All servers passed health checks")
            self.last_health_check = time.time()
        else:
            logger.warning(f"Only {healthy_count}/{len(self.servers)} servers are healthy")
            
        return all_healthy

    async def _cleanup_servers_sync(self):
        """Internal method to clean up servers synchronously."""
        cleanup_tasks = []
        for i, server in enumerate(self.servers):
            try:
                logger.info(f"Cleaning up server {i+1}/{len(self.servers)}")
                cleanup_tasks.append(server.cleanup())
            except Exception as e:
                logger.warning(f"Error initiating cleanup for server {i}: {e}")
        
        # Wait for all cleanup tasks to complete
        if cleanup_tasks:
            try:
                await asyncio.gather(*cleanup_tasks, return_exceptions=True)
            except Exception as e:
                logger.warning(f"Error during server cleanup gather: {e}")

    async def cleanup(self):
        """Clean up resources. Should only be called when completely done with the manager."""
        cleanup_successful = True
        
        # Clean up servers first
        logger.info(f"Starting cleanup of {len(self.servers)} servers")
        
        try:
            await self._cleanup_servers_sync()
            logger.info("Server cleanup completed")
        except Exception as e:
            cleanup_successful = False
            logger.error(f"Error during server cleanup: {e}", exc_info=True)

        self.servers = []

        # Then clean up agent
        if self.agent:
            try:
                logger.info("Cleaning up agent")
                await self.agent.cleanup()
                logger.info("Agent cleanup successful")
            except Exception as e:
                cleanup_successful = False
                logger.warning(f"Error during agent cleanup: {e}", exc_info=True)
            finally:
                self.agent = None
        
        logger.info(f"AutoStartManager cleanup complete (successful: {cleanup_successful})")
        return cleanup_successful

    def restart_servers(self) -> bool:
        """Restart all MCP servers. Returns True if successful."""
        try:
            logger.info("Restarting MCP servers...")
            
            # Clean up existing servers
            try:
                asyncio.run(self._cleanup_servers_sync())
            except Exception as e:
                logger.warning(f"Error during server cleanup before restart: {e}")
            
            self.servers = []
            
            # Setup new servers
            self.setup_servers()
            
            # Verify all servers are healthy
            if self.verify_all_servers_healthy():
                logger.info("Server restart successful")
                return True
            else:
                logger.error("Server restart failed health verification")
                return False
                
        except Exception as e:
            logger.error(f"Error restarting servers: {e}", exc_info=True)
            return False

    def check_and_recover_servers(self) -> bool:
        """Check server health and attempt recovery if needed."""
        try:
            if not self.servers:
                logger.info("No servers to check, setting up new servers")
                self.setup_servers()
                return len(self.servers) > 0
            
            # Quick health check
            if self.verify_all_servers_healthy(quick_check=True):
                return True
            
            logger.warning("Some servers are unhealthy, attempting recovery...")
            
            # Try individual server restart first
            server_names = ["data_extractor", "analysis"]
            recovered = True
            
            for i, server in enumerate(self.servers[:]):  # Use slice to avoid modification during iteration
                server_name = server_names[i] if i < len(server_names) else f"server_{i}"
                
                if not is_server_process_alive(server):
                    logger.warning(f"Server {server_name} is dead, restarting...")
                    try:
                        # Try to restart individual server
                        new_server = restart_server(server_name, self.server_env_vars)
                        if new_server and self._verify_server_health(new_server, server_name):
                            self.servers[i] = new_server
                            logger.info(f"Successfully restarted {server_name}")
                        else:
                            logger.error(f"Failed to restart {server_name}")
                            recovered = False
                    except Exception as e:
                        logger.error(f"Error restarting {server_name}: {e}")
                        recovered = False
            
            if not recovered:
                logger.warning("Individual server restart failed, trying full restart...")
                return self.restart_servers()
            
            return recovered
            
        except Exception as e:
            logger.error(f"Error during server recovery: {e}", exc_info=True)
            return False



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

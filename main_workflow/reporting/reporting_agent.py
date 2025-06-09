#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Reporting Agent Controller for RNA-seq Analysis

This module provides the agent controller that integrates MCP servers with LLMs
to generate intelligent reports and analyses of RNA-seq data. It coordinates
between data extraction, analysis, and natural language generation.
"""

import os
import sys
import json
import logging
import asyncio
from typing import List, Dict, Any, Optional, Tuple
from pathlib import Path
from datetime import datetime

# Import MCP server utilities from parent directory
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mcp_utils import setup_mcp_server

# Import agent framework
from pydantic_ai import Agent
from pydantic_ai.models import Model

# Set up logging
logging.basicConfig(
    level=logging.WARNING,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class ReportingAgent:
    """Controller for the RNA-seq reporting agent."""

    def __init__(self,
                 results_dir: str,
                 model_name: str = "openai:gpt-4o-mini",
                 temperature: float = 0.1):
        """
        Initialize the reporting agent.

        Args:
            results_dir: Path to UORCA results directory
            model_name: Model identifier (e.g., "openai:gpt-4o-mini", "ollama:llama3.3:latest")
            temperature: Model temperature for response generation
        """
        self.results_dir = os.path.abspath(results_dir)
        self.model_name = model_name
        self.temperature = temperature
        self.servers = []
        self.agent = None

        logger.info(f"Initialized ReportingAgent with results_dir: {self.results_dir}")
        logger.info(f"Using model: {self.model_name}")

    def setup_servers(self) -> List:
        """Set up the MCP servers for data extraction and analysis."""
        try:
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

    def _get_model_instance(self) -> Model:
        """Get the appropriate model instance based on model name."""
        if self.model_name.startswith("openai:"):
            # Using OpenAI model directly
            return self.model_name
        elif self.model_name.startswith("ollama:"):
            # Local Ollama model
            from pydantic_ai.models.openai import OpenAIModel
            model_name = self.model_name.replace("ollama:", "")
            return OpenAIModel(
                model=model_name,
                base_url='http://localhost:11434/v1',
                api_key="ollama"  # Ollama doesn't require API key
            )
        else:
            # Default to OpenAI
            return f"openai:{self.model_name}"

    def create_agent(self, research_question: str = None) -> Agent:
        """
        Create the agent with the specified research question.

        Args:
            research_question: Optional specific research question to focus on

        Returns:
            Configured Agent instance
        """
        # Set up servers if not already done
        if not self.servers:
            self.setup_servers()

        # Get model instance
        model_instance = self._get_model_instance()

        # Create the agent
        self.agent = Agent(
            model=model_instance,
            model_settings={"temperature": self.temperature},
            mcp_servers=self.servers,
            system=self._get_system_prompt(research_question)
        )

        logger.info("Created agent with MCP servers and system prompt")
        return self.agent

    def _get_system_prompt(self, research_question: str = None) -> str:
        """Generate the system prompt for the agent."""
        base_prompt = """You are an expert computational biologist specializing in RNA-seq data analysis and interpretation.

Your role is to help researchers understand their differential expression results by:
1. Identifying the most relevant contrasts for their research questions
2. Finding key differentially expressed genes and their patterns
3. Providing biological interpretation and mechanistic insights
4. Suggesting potential pathways and therapeutic implications

You have access to MCP tools for data extraction and analysis:
- list_available_contrasts: See all contrasts in the dataset
- get_deg_counts: Get differential expression statistics for contrasts
- find_common_degs: Find genes differentially expressed across multiple contrasts
- get_gene_info: Get detailed expression data for specific genes
- rank_contrasts_by_relevance: Rank contrasts by relevance to research questions
- analyze_gene_patterns: Analyze expression patterns across genes
- calculate_gene_statistics: Get comprehensive statistics for genes

Always provide data-driven insights backed by specific numbers and patterns from the analysis."""

        if research_question:
            base_prompt += f"\n\nCurrent research focus: '{research_question}'"
            base_prompt += "\n\nPrioritize findings and interpretations that directly address this research question."

        return base_prompt

    async def analyze_dataset_overview(self) -> Dict[str, Any]:
        """
        Generate an overview analysis of the entire dataset.

        Returns:
            Dictionary containing overview analysis results
        """
        if not self.agent:
            self.create_agent()

        prompt = """Please provide a comprehensive overview of this RNA-seq dataset:

1. First, list all available contrasts using the list_available_contrasts tool
2. Identify the contrasts with the most DEGs
3. Find genes that are consistently differentially expressed across multiple contrasts
4. Provide a summary of the key biological themes apparent in the data

Focus on providing actionable insights that would help a researcher understand what experiments were performed and what the main findings are."""

        try:
            result = await self.agent.run(prompt)
            return {
                "success": True,
                "analysis": result.data,
                "timestamp": datetime.now().isoformat()
            }
        except Exception as e:
            logger.error(f"Error in analyze_dataset_overview: {e}", exc_info=True)
            return {
                "success": False,
                "error": str(e),
                "timestamp": datetime.now().isoformat()
            }

    async def analyze_research_question(self,
                                      research_question: str,
                                      max_contrasts: int = 10,
                                      max_genes: int = 50) -> Dict[str, Any]:
        """
        Analyze the dataset with a specific research question in mind.

        Args:
            research_question: The biological question to investigate
            max_contrasts: Maximum number of contrasts to analyze
            max_genes: Maximum number of genes to report

        Returns:
            Dictionary containing analysis results
        """
        # Create agent with research question focus
        self.create_agent(research_question)

        prompt = f"""Analyze this RNA-seq dataset to address the research question: "{research_question}"

Please follow this structured approach:

1. Use rank_contrasts_by_relevance to identify the {max_contrasts} most relevant contrasts for this research question

2. For the top relevant contrasts:
   - Get DEG counts using get_deg_counts
   - Find common DEGs across these contrasts using find_common_degs
   - Focus on genes that appear in multiple relevant contrasts

3. For the top {max_genes} most important genes:
   - Use get_gene_info to get their expression patterns
   - Use analyze_gene_patterns to understand their relationships
   - Calculate statistics using calculate_gene_statistics

4. Synthesize findings into:
   - A direct answer to the research question
   - Key genes and their roles
   - Biological mechanisms and pathways involved
   - Potential therapeutic implications or future research directions

Be specific with numbers, gene names, and fold changes. Avoid generic statements."""

        try:
            result = await self.agent.run(prompt)
            return {
                "success": True,
                "research_question": research_question,
                "analysis": result.data,
                "parameters": {
                    "max_contrasts": max_contrasts,
                    "max_genes": max_genes
                },
                "timestamp": datetime.now().isoformat()
            }
        except Exception as e:
            logger.error(f"Error in analyze_research_question: {e}", exc_info=True)
            return {
                "success": False,
                "error": str(e),
                "research_question": research_question,
                "timestamp": datetime.now().isoformat()
            }

    async def find_biomarkers(self,
                            condition_keywords: List[str],
                            min_frequency: float = 0.8,
                            consistency_threshold: float = 0.9) -> Dict[str, Any]:
        """
        Find potential biomarker genes for specified conditions.

        Args:
            condition_keywords: Keywords to identify relevant contrasts
            min_frequency: Minimum frequency of differential expression
            consistency_threshold: Minimum consistency in direction

        Returns:
            Dictionary containing biomarker analysis
        """
        if not self.agent:
            self.create_agent()

        keywords_str = ", ".join(condition_keywords)
        prompt = f"""Identify potential biomarker genes for conditions related to: {keywords_str}

1. First, find all contrasts that might be relevant to these keywords
2. Use identify_consistent_degs with min_frequency={min_frequency} and consistency_threshold={consistency_threshold}
3. For the most consistent genes:
   - Verify they are relevant to the conditions of interest
   - Check their expression patterns across all relevant contrasts
   - Assess their potential as biomarkers

Provide a ranked list of potential biomarkers with:
- Consistency scores
- Expression patterns
- Biological relevance
- Literature support (if known)
- Potential clinical applications"""

        try:
            result = await self.agent.run(prompt)
            return {
                "success": True,
                "condition_keywords": condition_keywords,
                "biomarker_analysis": result.data,
                "parameters": {
                    "min_frequency": min_frequency,
                    "consistency_threshold": consistency_threshold
                },
                "timestamp": datetime.now().isoformat()
            }
        except Exception as e:
            logger.error(f"Error in find_biomarkers: {e}", exc_info=True)
            return {
                "success": False,
                "error": str(e),
                "condition_keywords": condition_keywords,
                "timestamp": datetime.now().isoformat()
            }

    async def generate_pathway_analysis(self,
                                      genes: List[str],
                                      research_context: str = None) -> Dict[str, Any]:
        """
        Generate pathway and network analysis for a set of genes.

        Args:
            genes: List of genes to analyze
            research_context: Optional context for the analysis

        Returns:
            Dictionary containing pathway analysis
        """
        if not self.agent:
            self.create_agent(research_context)

        genes_str = ", ".join(genes[:20])  # Limit to first 20 for prompt
        if len(genes) > 20:
            genes_str += f" and {len(genes) - 20} more"

        prompt = f"""Analyze potential pathways and networks for these genes: {genes_str}

1. Use analyze_gene_patterns to understand how these genes behave together
2. Use find_expression_clusters to identify co-expression modules
3. For each major pattern or cluster:
   - Describe the expression pattern
   - Suggest biological pathways that might be involved
   - Identify potential regulatory relationships

If research context is provided: {research_context or 'General analysis'}

Provide insights on:
- Major pathways likely affected
- Potential upstream regulators
- Downstream effects
- Therapeutic targeting opportunities"""

        try:
            result = await self.agent.run(prompt)
            return {
                "success": True,
                "genes_analyzed": len(genes),
                "pathway_analysis": result.data,
                "research_context": research_context,
                "timestamp": datetime.now().isoformat()
            }
        except Exception as e:
            logger.error(f"Error in generate_pathway_analysis: {e}", exc_info=True)
            return {
                "success": False,
                "error": str(e),
                "genes_analyzed": len(genes),
                "timestamp": datetime.now().isoformat()
            }

    async def compare_conditions(self,
                               condition1_keywords: List[str],
                               condition2_keywords: List[str]) -> Dict[str, Any]:
        """
        Compare gene expression between two conditions.

        Args:
            condition1_keywords: Keywords for first condition
            condition2_keywords: Keywords for second condition

        Returns:
            Dictionary containing comparison analysis
        """
        if not self.agent:
            self.create_agent()

        prompt = f"""Compare gene expression between two conditions:
Condition 1: {', '.join(condition1_keywords)}
Condition 2: {', '.join(condition2_keywords)}

1. Identify contrasts relevant to each condition
2. Find genes that are:
   - Uniquely changed in condition 1
   - Uniquely changed in condition 2
   - Changed in both but in opposite directions
   - Changed in both in the same direction

3. For the most significant differences:
   - Analyze the biological implications
   - Suggest what cellular processes differ between conditions
   - Identify potential therapeutic targets

Provide a clear comparison with specific genes and pathways that distinguish these conditions."""

        try:
            result = await self.agent.run(prompt)
            return {
                "success": True,
                "condition1": condition1_keywords,
                "condition2": condition2_keywords,
                "comparison_analysis": result.data,
                "timestamp": datetime.now().isoformat()
            }
        except Exception as e:
            logger.error(f"Error in compare_conditions: {e}", exc_info=True)
            return {
                "success": False,
                "error": str(e),
                "timestamp": datetime.now().isoformat()
            }

    async def cleanup(self):
        """Clean up resources and close MCP servers."""
        for server in self.servers:
            try:
                # MCPServerStdio uses cleanup() to terminate the subprocess
                await server.cleanup()
            except Exception as e:
                logger.warning(f"Error cleaning up server: {e}", exc_info=True)

        self.servers = []
        self.agent = None
        logger.info("Cleaned up reporting agent resources")


# Convenience functions for common analyses
async def quick_analysis(results_dir: str,
                        research_question: str,
                        model: str = "openai:gpt-4o-mini") -> Dict[str, Any]:
    """
    Perform a quick analysis of RNA-seq data for a research question.

    Args:
        results_dir: Path to UORCA results
        research_question: The research question to investigate
        model: Model to use for analysis

    Returns:
        Analysis results dictionary
    """
    agent = ReportingAgent(results_dir, model_name=model)
    try:
        return await agent.analyze_research_question(research_question)
    finally:
        try:
            agent.terminate()
        except Exception as e:
            logger.warning(f"Error during quick_analysis cleanup: {e}", exc_info=True)


async def dataset_summary(results_dir: str,
                         model: str = "openai:gpt-4o-mini") -> Dict[str, Any]:
    """
    Generate a summary overview of an RNA-seq dataset.

    Args:
        results_dir: Path to UORCA results
        model: Model to use for analysis

    Returns:
        Dataset overview dictionary
    """
    agent = ReportingAgent(results_dir, model_name=model)
    try:
        return await agent.analyze_dataset_overview()
    finally:
        try:
            agent.terminate()
        except Exception as e:
            logger.warning(f"Error during dataset_summary cleanup: {e}", exc_info=True)


if __name__ == "__main__":
    # Example usage
    import asyncio

    async def main():
        results_dir = os.environ.get("UORCA_RESULTS", "./UORCA_results")

        # Example 1: Dataset overview
        print("Generating dataset overview...")
        overview = await dataset_summary(results_dir)
        print(json.dumps(overview, indent=2))

        # Example 2: Research question analysis
        print("\nAnalyzing specific research question...")
        analysis = await quick_analysis(
            results_dir,
            "What genes are involved in the inflammatory response?"
        )
        print(json.dumps(analysis, indent=2))

    asyncio.run(main())

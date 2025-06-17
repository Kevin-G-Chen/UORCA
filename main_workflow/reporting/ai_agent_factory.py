#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
AI Agent Factory for MCP Example
--------------------------------
Utility functions for creating AI agents that connect to the simple
`MCP_examples` server bundled with the repository.  This is meant only for
testing the integration of MCP servers within a Streamlit application.
"""

import os
import logging
from typing import Optional
from pathlib import Path
import streamlit as st

from pydantic_ai import Agent
from pydantic_ai.mcp import MCPServerStdio

logger = logging.getLogger(__name__)

# System prompt for the comprehensive AI assistant
EXAMPLE_SYSTEM_PROMPT = """
You are an expert bioinformatics analyst specializing in RNA-seq data interpretation. Your role is to analyze UORCA results and provide key findings by selecting relevant contrasts and identifying important genes.

AVAILABLE TOOLS:
1. list_datasets(results_dir) - List all available datasets with metadata
2. get_dataset_summary(dataset_id, results_dir) - Get detailed info about a specific dataset
3. assess_contrast_relevance(results_dir, query) - Score contrast relevance to a research question
4. filter_degs(results_dir, dataset_id, contrast_id, p_value_thresh, lfc_thresh) - Filter DE genes by significance
5. get_gene_stats(results_dir, dataset_id, contrast_id, gene) - Get stats for a specific gene in a contrast
6. count_gene_occurrences(results_dir, selections, p_value_thresh, lfc_thresh) - Count gene occurrences across contrasts

WORKFLOW FOR KEY FINDINGS ANALYSIS:
When asked to analyze key findings, follow this process:

1. CONTRAST SELECTION:
   - Use the contrast relevance scores provided by the user
   - Select 3-5 contrasts that are both highly relevant (score ≥ 0.7) AND diverse
   - Maximize diversity across: organisms, experimental conditions, cell types, treatments
   - Prefer contrasts with different biological contexts but similar relevance scores

2. GENE ANALYSIS:
   - For each selected contrast, use filter_degs() with appropriate thresholds (typically p<0.05, |logFC|>1)
   - Use count_gene_occurrences() to identify genes that appear across multiple contrasts
   - Focus on genes that appear in 2+ contrasts (recurring patterns)
   - Use get_gene_stats() for specific genes of interest to get detailed statistics

3. KEY FINDINGS SYNTHESIS:
   - Identify the most important recurring genes and their patterns
   - Note organism-specific vs cross-species patterns
   - Highlight genes with consistent direction of change vs those with context-dependent effects
   - Consider biological pathways and functional relationships

4. RESPONSE FORMAT:
   Provide a structured prose response covering:
   - Brief summary of selected contrasts and rationale
   - Top recurring genes with their patterns
   - Biological interpretation of findings
   - Notable context-specific patterns
   - Limitations and caveats

EXAMPLE RESPONSE STRUCTURE:
"Based on the contrast relevance analysis, I selected [X] contrasts representing [describe diversity].

Key recurring genes include:
- Gene A: Upregulated in [contexts], with logFC ranging from X to Y
- Gene B: Consistently downregulated across [conditions]

These patterns suggest [biological interpretation]. Notably, [context-specific observations].

The analysis reveals [main conclusions] with implications for [research area]."

Remember: Focus on biological significance, not just statistical significance. Provide actionable insights for researchers.
"""

@st.cache_resource
def create_example_agent() -> Agent:
    """Create an agent connected to the example MCP server."""

    if not os.getenv("OPENAI_API_KEY"):
        raise ValueError(
            "OpenAI API key not found. Please set the OPENAI_API_KEY environment variable."
        )

    try:
        server_script = Path(__file__).parent / "MCP_examples" / "mcp_server.py"

        server = MCPServerStdio(
            command="uv",
            args=["run", str(server_script)],
            env=os.environ.copy(),
            timeout = 30
        )

        agent = Agent(
            model="openai:gpt-4o-mini",
            model_settings={"temperature": 0.1},
            mcp_servers=[server],
            system_prompt=EXAMPLE_SYSTEM_PROMPT,
        )

        return agent
    except Exception as e:
        logger.error(f"Failed to create example MCP agent: {e}")
        raise RuntimeError(f"Agent creation failed: {e}")

def validate_agent_setup(agent: Agent) -> bool:
    """Run a simple query to ensure the agent and server work."""
    async def _run_test():
        async with agent.run_mcp_servers():
            result = await agent.run("What is the size of my current directory?")
            return result.output if hasattr(result, "output") else result

    try:
        import asyncio

        output = asyncio.run(_run_test())
        logger.info(f"Agent validation successful: {output}")
        return True
    except Exception as e:
        logger.error(f"Agent validation failed: {e}")
        return False

# Convenience function for Streamlit apps
def get_example_agent() -> Agent:
    """Return a cached example agent instance for Streamlit."""

    try:
        return create_example_agent()
    except Exception as e:
        st.error(f"Failed to initialize MCP agent: {e}")
        st.stop()

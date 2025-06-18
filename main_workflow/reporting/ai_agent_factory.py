#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
AI Agent Factory for UORCA Analysis
-----------------------------------
Factory functions for creating AI agents that connect to the UORCA MCP server
for analyzing RNA-seq results using structured tools.
"""

import os
import logging
from typing import Optional
from pathlib import Path
import streamlit as st

from pydantic_ai import Agent
from pydantic_ai.mcp import MCPServerStdio

logger = logging.getLogger(__name__)

# System prompt for the UORCA analysis agent
UORCA_SYSTEM_PROMPT = """
You are an expert bioinformatics analyst with access to four specialized tools for analyzing RNA-seq differential expression data:

1) get_most_common_genes(lfc_thresh, p_thresh, top_n) - Find genes that are differentially expressed across the most contrasts
2) get_gene_contrast_stats(gene, contrast_id?) - Get detailed statistics for a specific gene across contrasts
3) filter_genes_by_contrast_sets(set_a, set_b, lfc_thresh, p_thresh) - Find genes significant in set A but not set B
4) summarize_contrast(contrast_id, lfc_thresh, p_thresh, max_genes) - Summarize top DEGs and statistics for a contrast

WORKFLOW:
- You will receive a research question and a list of pre-selected contrasts that are relevant to that question.
- Choose reasonable thresholds yourself (typically lfc_thresh around 1.0-2.0, p_thresh around 0.01-0.05).
- Use get_most_common_genes() to identify recurring differential expression patterns across contrasts.
- Use filter_genes_by_contrast_sets() to find genes that are specific to particular experimental contexts.
- Drill into specific genes of interest using get_gene_contrast_stats().
- Optionally use summarize_contrast() to get overviews of individual contrasts.

OUTPUT FORMAT:
Return your analysis in two parts:
1) A structured summary showing the key genes identified for each contrast or gene set
2) A brief 2-3 sentence biological interpretation explaining your rationale and what patterns you discovered

Focus on identifying both shared signatures across multiple contrasts and context-specific gene expression patterns that could provide biological insights.
"""

@st.cache_resource
def create_uorca_agent() -> Agent:
    """Create an agent connected to the UORCA MCP server."""

    if not os.getenv("OPENAI_API_KEY"):
        raise ValueError(
            "OpenAI API key not found. Please set the OPENAI_API_KEY environment variable."
        )

    try:
        server_script = Path(__file__).parent / "mcp_server_core.py"

        server = MCPServerStdio(
            command="uv",
            args=["run", str(server_script), "server"],
            env=os.environ.copy(),
            timeout=30
        )

        agent = Agent(
            model="openai:gpt-4o-mini",
            model_settings={"temperature": 0.1},
            mcp_servers=[server],
            system_prompt=UORCA_SYSTEM_PROMPT,
        )

        return agent
    except Exception as e:
        logger.error(f"Failed to create UORCA MCP agent: {e}")
        raise RuntimeError(f"Agent creation failed: {e}")

def validate_agent_setup(agent: Agent) -> bool:
    """Run a simple query to ensure the agent and server work."""
    async def _run_test():
        async with agent.run_mcp_servers():
            result = await agent.run("Test the get_most_common_genes tool with lfc_thresh=1.0, p_thresh=0.05, top_n=5")
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
def get_uorca_agent() -> Agent:
    """Return a cached UORCA agent instance for Streamlit."""
    try:
        return create_uorca_agent()
    except Exception as e:
        st.error(f"Failed to initialize UORCA MCP agent: {e}")
        st.stop()

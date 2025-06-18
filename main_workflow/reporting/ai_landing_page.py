#!/usr/bin/env python
"""Simplified AI landing page demonstrating MCP integration."""

import asyncio
import logging
from pathlib import Path

import streamlit as st

from ai_agent_factory import get_uorca_agent, validate_agent_setup

logger = logging.getLogger(__name__)


def render_ai_landing_page(results_dir: str) -> None:
    """Render a minimal interface that runs the MCP example."""
    st.title("🤖 Dataset Explorer")
    st.markdown(
        "This tab demonstrates running an MCP server to explore dataset information."
    )

    try:
        agent = get_uorca_agent()
    except Exception as e:
        st.error(f"Failed to initialise agent: {e}")
        return

    # Single button to run the query like in mcp_client.py
    if st.button("Run Dataset Explorer"):
        run_example(agent, results_dir)


def run_example(agent, results_dir: str) -> None:
    """Run the query to get information about the first dataset."""
    async def _query():
        async with agent.run_mcp_servers():
            result = await agent.run(f"Test the get_most_common_genes tool with lfc_thresh=1.0, p_thresh=0.05, top_n=5")
            return result.output if hasattr(result, "output") else result

    with st.spinner("Getting dataset information..."):
        try:
            output = asyncio.run(_query())
            st.write(output)
        except Exception as exc:
            logger.error("Dataset query failed: %s", exc)
            st.error(f"Query failed: {exc}")

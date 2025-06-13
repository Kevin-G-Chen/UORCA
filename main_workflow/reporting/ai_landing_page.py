#!/usr/bin/env python
"""Simplified AI landing page demonstrating MCP integration."""

import asyncio
import logging
from pathlib import Path

import streamlit as st

from ai_agent_factory import get_example_agent, validate_agent_setup

logger = logging.getLogger(__name__)


def render_ai_landing_page(results_dir: str) -> None:
    """Render a minimal interface that runs the MCP example."""
    st.title("🤖 MCP Example")
    st.markdown(
        "This tab demonstrates running an MCP server from a Streamlit app."
    )

    try:
        agent = get_example_agent()
    except Exception as e:
        st.error(f"Failed to initialise agent: {e}")
        return

    if st.button("Run MCP Example"):
        run_example(agent)


def run_example(agent) -> None:
    """Run the example query against the MCP server."""

    async def _query():
        async with agent.run_mcp_servers():
            result = await agent.run("What is the size of my current directory?")
            return result.output if hasattr(result, "output") else result

    with st.spinner("Running example query..."):
        try:
            output = asyncio.run(_query())
            st.write(output)
        except Exception as exc:
            logger.error("Example query failed: %s", exc)
            st.error(f"Query failed: {exc}")


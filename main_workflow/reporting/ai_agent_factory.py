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

# System prompt used for the simple MCP example
EXAMPLE_SYSTEM_PROMPT = (
    "You are an assistant that helps to describe datasets. Only use tools that will help you describe datasets."
)

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

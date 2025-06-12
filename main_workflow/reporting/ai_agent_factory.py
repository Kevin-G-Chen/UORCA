#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
AI Agent Factory for UORCA Explorer
----------------------------------
Creates and configures AI agents for interacting with UORCA RNA-seq analysis results.
"""

import os
import logging
from typing import Optional
import streamlit as st

from pydantic_ai import Agent
from pydantic_ai.models.openai import OpenAIModel

# Import MCP utilities
from mcp_utils import setup_mcp_server

logger = logging.getLogger(__name__)

# System prompt for UORCA data analysis
UORCA_SYSTEM_PROMPT = """
You are an expert bioinformatics assistant specializing in RNA-seq data analysis.
You have access to UORCA (Upstream Open Reading frame Comparative Analysis) results through specialized tools.

Your tools allow you to:
- list_contrasts(): Get all available experimental contrasts
- list_genes(): Get available gene symbols
- get_lfc(contrast, gene): Get log fold change for a specific gene in a contrast
- get_analysis_info(): Get metadata about the analyses

Key guidelines:
1. Always provide clear, scientifically accurate responses
2. When discussing fold changes, explain the biological significance
3. For log fold changes: positive values = upregulated, negative = downregulated
4. If a gene/contrast combination isn't found, suggest similar alternatives when possible
5. Keep responses concise but informative
6. When listing multiple items, present them in an organized, readable format

Context: UORCA focuses on upstream open reading frames (uORFs) which are regulatory elements
that can affect translation efficiency and gene expression.
"""

@st.cache_resource
def create_uorca_agent(results_dir: Optional[str] = None) -> Agent:
    """
    Create and configure an AI agent for UORCA data analysis.

    Args:
        results_dir: Path to UORCA results directory (optional, uses env var if not provided)

    Returns:
        Configured Agent instance

    Raises:
        ValueError: If OpenAI API key is not available
        RuntimeError: If MCP server setup fails
    """
    # Check for OpenAI API key
    if not os.getenv("OPENAI_API_KEY"):
        raise ValueError(
            "OpenAI API key not found. Please set the OPENAI_API_KEY environment variable."
        )

    try:
        # Setup MCP server with results directory
        env_vars = {}
        if results_dir:
            env_vars["UORCA_RESULTS_DIR"] = results_dir

        server = setup_mcp_server("uorca_data", env_vars=env_vars)
        logger.info(f"Successfully set up UORCA data MCP server")

        # Create the agent with OpenAI GPT-4o-mini
        agent = Agent(
            model="openai:gpt-4o-mini",
            model_settings={
                "temperature": 0.1,  # Low temperature for consistent, factual responses
                "max_tokens": 1000,  # Reasonable limit for responses
            },
            mcp_servers=[server],
            system_prompt=UORCA_SYSTEM_PROMPT,
        )

        logger.info("Successfully created UORCA AI agent")
        return agent

    except Exception as e:
        logger.error(f"Failed to create UORCA agent: {e}")
        raise RuntimeError(f"Agent creation failed: {e}")

def validate_agent_setup(agent: Agent) -> bool:
    """
    Validate that the agent is properly configured and can access data.

    Args:
        agent: The agent to validate

    Returns:
        True if agent is working, False otherwise
    """
    try:
        # Simple test query to check if tools are accessible
        import asyncio
        result = asyncio.run(agent.run("How many contrasts are available? Just give me the number."))
        logger.info(f"Agent validation successful: {result}")
        return True
    except Exception as e:
        logger.error(f"Agent validation failed: {e}")
        return False

# Convenience function for Streamlit apps
def get_uorca_agent(results_dir: Optional[str] = None) -> Agent:
    """
    Get a cached UORCA agent instance.

    This is a convenience wrapper around create_uorca_agent that provides
    better error handling for Streamlit applications.

    Args:
        results_dir: Path to UORCA results directory

    Returns:
        Configured Agent instance
    """
    try:
        return create_uorca_agent(results_dir)
    except Exception as e:
        st.error(f"Failed to initialize AI agent: {e}")
        st.info("Please ensure:")
        st.info("• OpenAI API key is set in environment variables")
        st.info("• UORCA results directory contains valid analysis data")
        st.info("• All required dependencies are installed")
        st.stop()

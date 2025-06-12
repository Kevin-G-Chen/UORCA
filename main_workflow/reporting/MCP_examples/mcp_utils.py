#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Utility functions for working with MCP servers in Autonify.
"""

import os
import sys
from pathlib import Path
import logging
from typing import Optional, Dict

from pydantic_ai.mcp import MCPServerStdio

# Import our config
from autonify.config import ServerConfig

logger = logging.getLogger("autonify.mcp_utils")

def setup_mcp_server(server_name: str, env_vars: Optional[Dict[str, str]] = None) -> MCPServerStdio:
    """Set up an MCP server using the stdio transport.

    Args:
        server_name: The name of the server in servers.toml (e.g., "web_search", "docling")
        env_vars: Optional dictionary of environment variables to pass to the server process

    Returns:
        Configured MCPServerStdio instance

    Raises:
        ValueError: If the server configuration cannot be found
    """
    # Get server path from config
    config = ServerConfig()
    server_config = config.get_server_config(server_name)

    if not server_config:
        raise ValueError(f"Server configuration for '{server_name}' not found in servers.toml")

    # Get the path and make it absolute
    server_path = Path(server_config["path"])
    if not server_path.is_absolute():
        # Make path relative to project root
        project_root = Path(__file__).parent.parent.parent
        server_path = project_root / server_path

    # Set up environment variables
    env_dict = dict(os.environ)

    # Add API keys if available
    if api_key := os.getenv('OPENAI_API_KEY'):
        env_dict['OPENAI_API_KEY'] = api_key

    # Add any additional environment variables
    if env_vars:
        env_dict.update(env_vars)

    # Create and return the server
    server = MCPServerStdio(
        sys.executable,
        args=[str(server_path)],  # Convert Path to string
        env=env_dict
    )

    logger.info(f"Configured MCP server '{server_name}': {server_path}")
    return server

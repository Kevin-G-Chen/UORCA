#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
MCP Utilities for UORCA
-----------------------
Utility functions for working with MCP servers in the UORCA project.
"""

import os
import sys
import toml
from pathlib import Path
import logging
from typing import Optional, Dict, Any, List

from pydantic_ai.mcp import MCPServerStdio

logger = logging.getLogger("uorca.mcp_utils")

def find_project_root() -> Path:
    """Find the UORCA project root directory by looking for servers.toml."""
    current = Path(__file__).parent

    # Walk up the directory tree to find servers.toml
    for parent in [current] + list(current.parents):
        if (parent / "servers.toml").exists():
            return parent

    # If not found, assume the current file's grandparent is project root
    return Path(__file__).parent.parent.parent

def load_server_config(server_name: str) -> Dict[str, Any]:
    """
    Load server configuration from servers.toml.

    Args:
        server_name: Name of the server to load config for

    Returns:
        Server configuration dictionary

    Raises:
        ValueError: If server configuration not found
        FileNotFoundError: If servers.toml not found
    """
    project_root = find_project_root()
    servers_file = project_root / "servers.toml"

    if not servers_file.exists():
        raise FileNotFoundError(f"servers.toml not found at {servers_file}")

    try:
        with open(servers_file, 'r') as f:
            config = toml.load(f)
    except Exception as e:
        raise ValueError(f"Failed to parse servers.toml: {e}")

    if "servers" not in config:
        raise ValueError("No 'servers' section found in servers.toml")

    if server_name not in config["servers"]:
        available_servers = list(config["servers"].keys())
        raise ValueError(
            f"Server '{server_name}' not found in servers.toml. "
            f"Available servers: {available_servers}"
        )

    server_config = config["servers"][server_name]

    # Check if server is enabled
    if not server_config.get("enabled", True):
        raise ValueError(f"Server '{server_name}' is disabled in configuration")

    return server_config

def setup_mcp_server(
    server_name: str,
    env_vars: Optional[Dict[str, str]] = None,
    server_args: Optional[List[str]] = None,
) -> MCPServerStdio:
    """
    Set up an MCP server using the stdio transport.

    Args:
        server_name: The name of the server in servers.toml (e.g., "uorca_data")
        env_vars: Optional dictionary of environment variables to pass to the server process
        server_args: Additional command line arguments to pass to the server

    Returns:
        Configured MCPServerStdio instance

    Raises:
        ValueError: If the server configuration cannot be found or is invalid
        FileNotFoundError: If the server script doesn't exist
    """
    # Load server configuration
    server_config = load_server_config(server_name)

    if "path" not in server_config:
        raise ValueError(f"No 'path' specified for server '{server_name}'")

    # Resolve server path
    project_root = find_project_root()
    server_path = Path(server_config["path"])

    if not server_path.is_absolute():
        server_path = project_root / server_path

    # Verify the server script exists
    if not server_path.exists():
        raise FileNotFoundError(f"Server script not found: {server_path}")

    # Set up environment variables
    env_dict = dict(os.environ)

    # Add OpenAI API key if available
    if api_key := os.getenv('OPENAI_API_KEY'):
        env_dict['OPENAI_API_KEY'] = api_key

    # Add any server-specific environment variables from config
    if "env" in server_config:
        env_dict.update(server_config["env"])

    # Add any additional environment variables passed to this function
    if env_vars:
        env_dict.update(env_vars)

    # Create and return the server
    args_list = [str(server_path)]
    if server_args:
        args_list.extend(server_args)

    server = MCPServerStdio(
        sys.executable,
        args=args_list,
        env=env_dict
    )

    logger.info(f"Configured MCP server '{server_name}': {server_path}")
    return server

def validate_mcp_setup() -> bool:
    """
    Validate that MCP setup is correct for UORCA.

    Returns:
        True if setup is valid, False otherwise
    """
    try:
        # Check if servers.toml exists
        project_root = find_project_root()
        servers_file = project_root / "servers.toml"

        if not servers_file.exists():
            logger.error(f"servers.toml not found at {servers_file}")
            return False

        # Try to load the config
        with open(servers_file, 'r') as f:
            config = toml.load(f)

        # Check if uorca_data server is configured
        if "servers" not in config or "uorca_data" not in config["servers"]:
            logger.error("uorca_data server not configured in servers.toml")
            return False

        # Check if server script exists
        server_config = config["servers"]["uorca_data"]
        server_path = Path(server_config["path"])

        if not server_path.is_absolute():
            server_path = project_root / server_path

        if not server_path.exists():
            logger.error(f"Server script not found: {server_path}")
            return False

        logger.info("MCP setup validation passed")
        return True

    except Exception as e:
        logger.error(f"MCP setup validation failed: {e}")
        return False

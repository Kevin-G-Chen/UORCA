#!/usr/bin/env python3
"""Utilities for launching MCP servers used in UORCA tests."""

from __future__ import annotations

import os
import sys
import subprocess
from pathlib import Path
from typing import Dict, Optional

import tomllib


class MCPServer:
    """Simple wrapper around a server subprocess."""

    def __init__(self, process: subprocess.Popen[str]):
        self._process = process

    def cleanup(self) -> None:
        """Terminate the server process if it is still running."""
        if self._process.poll() is None:
            self._process.terminate()
            try:
                self._process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                self._process.kill()


def _load_config(config_path: Optional[str] = None) -> Dict[str, Dict[str, str]]:
    if config_path is None:
        config_path = Path(__file__).parent / "reporting" / "servers.toml"
    with open(config_path, "rb") as fh:
        return tomllib.load(fh)


def setup_mcp_server(
    name: str,
    env_vars: Optional[Dict[str, str]] = None,
    config_path: Optional[str] = None,
) -> MCPServer:
    """Start an MCP server defined in ``servers.toml``.

    Parameters
    ----------
    name:
        Name of the server as defined in the configuration file.
    env_vars:
        Optional environment variables to set for the server process.
    config_path:
        Optional path to an alternative server configuration file.
    """
    config = _load_config(config_path)
    servers = config.get("servers", {})
    if name not in servers:
        raise ValueError(f"Server '{name}' not found in configuration")

    info = servers[name]
    if not info.get("enabled", False):
        raise ValueError(f"Server '{name}' is disabled in configuration")

    script_path = Path(__file__).parent / info["path"]
    if not script_path.exists():
        raise FileNotFoundError(f"MCP server script not found: {script_path}")

    env = os.environ.copy()
    if env_vars:
        env.update(env_vars)

    process = subprocess.Popen(
        [sys.executable, str(script_path)],
        env=env,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    return MCPServer(process)


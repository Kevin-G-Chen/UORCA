from __future__ import annotations

import os
from pathlib import Path
import tomllib
from typing import Mapping, Optional

from pydantic_ai.mcp import MCPServerStdio


_CONFIG_CACHE: Optional[dict] = None


def _load_config() -> dict:
    """Load servers.toml configuration."""
    global _CONFIG_CACHE
    if _CONFIG_CACHE is None:
        cfg_path = Path(__file__).with_name("servers.toml")
        with open(cfg_path, "rb") as f:
            _CONFIG_CACHE = tomllib.load(f).get("servers", {})
    return _CONFIG_CACHE


def setup_mcp_server(name: str, env_vars: Mapping[str, str] | None = None) -> MCPServerStdio:
    """Set up an MCP server defined in ``servers.toml``.

    Parameters
    ----------
    name:
        Name of the server configuration in ``servers.toml``.
    env_vars:
        Environment variables to provide to the subprocess.

    Returns
    -------
    MCPServerStdio
        Handle to the started MCP server which exposes ``cleanup()`` for
        termination.
    """
    servers = _load_config()
    cfg = servers.get(name)
    if not cfg:
        raise ValueError(f"Server '{name}' not found in servers.toml")
    if not cfg.get("enabled", False):
        raise ValueError(f"Server '{name}' is disabled")

    script_path = Path(__file__).parent / cfg["path"]
    if not script_path.exists():
        raise FileNotFoundError(f"Server script not found: {script_path}")

    env = os.environ.copy()
    if env_vars:
        env.update(env_vars)

    # Spawn server using python interpreter
    return MCPServerStdio(
        "python",
        args=[str(script_path), "server"],
        env=env,
    )

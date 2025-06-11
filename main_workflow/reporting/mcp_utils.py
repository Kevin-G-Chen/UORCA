from __future__ import annotations

import os
import sys
from pathlib import Path
import logging
import time
import asyncio
try:
    import tomllib  # Python 3.11+
except ImportError:
    import tomli as tomllib  # Python 3.9-3.10 fallback
from typing import Mapping, Optional

from pydantic_ai.mcp import MCPServerStdio


logger = logging.getLogger(__name__)


_CONFIG_CACHE: Optional[dict] = None


def _load_config() -> dict:
    """Load servers.toml configuration."""
    global _CONFIG_CACHE
    if _CONFIG_CACHE is None:
        cfg_path = Path(__file__).with_name("servers.toml")
        with open(cfg_path, "rb") as f:
            _CONFIG_CACHE = tomllib.load(f).get("servers", {})
    return _CONFIG_CACHE


async def check_server_health(server: MCPServerStdio, timeout: float = 2.0) -> bool:
    """Check if an MCP server is healthy and responsive.
    
    Parameters
    ----------
    server : MCPServerStdio
        The server to check
    timeout : float
        Timeout in seconds for the health check
        
    Returns
    -------
    bool
        True if server is healthy, False otherwise
    """
    try:
        # Check if the underlying process is still running
        if hasattr(server, '_process') and server._process:
            if server._process.poll() is not None:
                logger.warning("MCP server process has terminated")
                return False
        
        # Try to get available tools as a health check
        # This is a simple way to verify the server is responding
        try:
            tools = await asyncio.wait_for(server.list_tools(), timeout=timeout)
            logger.debug(f"Server health check passed: {len(tools)} tools available")
            return True
        except asyncio.TimeoutError:
            logger.warning("MCP server health check timed out")
            return False
        except Exception as e:
            logger.warning(f"MCP server health check failed: {e}")
            return False
            
    except Exception as e:
        logger.error(f"Error during server health check: {e}")
        return False


def setup_mcp_server(name: str, env_vars: Mapping[str, str] | None = None, 
                     max_retries: int = 3, retry_delay: float = 1.0) -> MCPServerStdio:
    """Set up an MCP server defined in ``servers.toml`` with robust error handling.

    Parameters
    ----------
    name:
        Name of the server configuration in ``servers.toml``.
    env_vars:
        Environment variables to provide to the subprocess.
    max_retries:
        Maximum number of startup attempts
    retry_delay:
        Delay between retry attempts in seconds

    Returns
    -------
    MCPServerStdio
        Handle to the started MCP server which exposes ``cleanup()`` for
        termination.
        
    Raises
    ------
    RuntimeError
        If server fails to start after all retry attempts
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

    results_dir = env_vars.get("RESULTS_DIR", "") if env_vars else ""
    args = [str(script_path), "server"]
    if results_dir:
        args.extend(["--results-dir", results_dir])
        logger.info(f"Passing results dir as argument: {results_dir}")

    last_error = None
    
    for attempt in range(max_retries):
        try:
            logger.info(f"Starting MCP server '{name}' (attempt {attempt + 1}/{max_retries})")
            
            # Create the server
            server = MCPServerStdio(
                sys.executable,
                args=args,
                env=env,
            )
            
            # Wait for server to initialize
            time.sleep(retry_delay)
            
            # Verify server is working with a simple health check
            try:
                # Create event loop for health check if needed
                try:
                    loop = asyncio.get_event_loop()
                    if loop.is_running():
                        # If loop is already running, we can't use run_until_complete
                        # Just return the server and hope for the best
                        logger.warning("Event loop already running, skipping health check")
                        return server
                except RuntimeError:
                    loop = asyncio.new_event_loop()
                    asyncio.set_event_loop(loop)
                
                try:
                    # Run health check
                    is_healthy = loop.run_until_complete(
                        asyncio.wait_for(check_server_health(server), timeout=5.0)
                    )
                    
                    if is_healthy:
                        logger.info(f"MCP server '{name}' started successfully")
                        return server
                    else:
                        raise RuntimeError("Server health check failed")
                        
                finally:
                    # Clean up the loop if we created it
                    try:
                        if not loop.is_running():
                            loop.close()
                    except Exception:
                        pass
                        
            except Exception as health_error:
                logger.warning(f"Health check failed for server '{name}': {health_error}")
                # Clean up the failed server
                try:
                    # Use asyncio.run for cleanup since we're not in an async context
                    def cleanup_server():
                        try:
                            loop = asyncio.new_event_loop()
                            asyncio.set_event_loop(loop)
                            loop.run_until_complete(server.cleanup())
                            loop.close()
                        except Exception as cleanup_err:
                            logger.warning(f"Error during server cleanup: {cleanup_err}")
                    cleanup_server()
                except Exception:
                    pass
                raise health_error
                
        except Exception as e:
            last_error = e
            logger.warning(f"Failed to start MCP server '{name}' on attempt {attempt + 1}: {e}")
            
            # Clean up any partial server state
            try:
                if 'server' in locals():
                    # Use asyncio.run for cleanup since we're not in an async context
                    def cleanup_server():
                        try:
                            loop = asyncio.new_event_loop()
                            asyncio.set_event_loop(loop)
                            loop.run_until_complete(server.cleanup())
                            loop.close()
                        except Exception as cleanup_err:
                            logger.warning(f"Error during server cleanup: {cleanup_err}")
                    cleanup_server()
            except Exception:
                pass
                
            # Wait before retrying (except on last attempt)
            if attempt < max_retries - 1:
                logger.info(f"Waiting {retry_delay} seconds before retry...")
                time.sleep(retry_delay)
    
    # If we get here, we've exhausted all retries
    error_msg = f"Failed to start MCP server '{name}' after {max_retries} attempts"
    if last_error:
        error_msg += f". Last error: {last_error}"
    logger.error(error_msg)
    raise RuntimeError(error_msg)


def verify_server_responsive(server: MCPServerStdio, timeout: float = 3.0) -> bool:
    """Verify that an MCP server is responsive.
    
    This is a synchronous wrapper around the health check for use in non-async contexts.
    
    Parameters
    ----------
    server : MCPServerStdio
        The server to verify
    timeout : float
        Timeout for the verification
        
    Returns
    -------
    bool
        True if server is responsive, False otherwise
    """
    try:
        # Check if the underlying process is still running first
        if hasattr(server, '_process') and server._process:
            if server._process.poll() is not None:
                logger.warning("MCP server process has terminated")
                return False
        
        # Try to use existing event loop
        try:
            loop = asyncio.get_event_loop()
            if loop.is_running():
                # In Streamlit context, we can't use run_until_complete
                # Do a simpler check - just verify process is alive
                logger.debug("Event loop is running, doing basic process check")
                return hasattr(server, '_process') and server._process and server._process.poll() is None
        except RuntimeError:
            # No existing loop, create a new one
            pass
        
        # Create new event loop for health check
        try:
            loop = asyncio.new_event_loop()
            asyncio.set_event_loop(loop)
            result = loop.run_until_complete(
                asyncio.wait_for(check_server_health(server), timeout=timeout)
            )
            return result
        finally:
            try:
                loop.close()
                # Reset event loop policy to avoid issues
                asyncio.set_event_loop(None)
            except Exception:
                pass
                
    except Exception as e:
        logger.error(f"Error verifying server responsiveness: {e}")
        return False


def is_server_process_alive(server: MCPServerStdio) -> bool:
    """Quick check if server process is still alive.
    
    Parameters
    ----------
    server : MCPServerStdio
        The server to check
        
    Returns
    -------
    bool
        True if process is alive, False otherwise
    """
    try:
        if hasattr(server, '_process') and server._process:
            return server._process.poll() is None
        return False
    except Exception as e:
        logger.error(f"Error checking server process: {e}")
        return False


def restart_server(name: str, env_vars: Mapping[str, str] | None = None) -> Optional[MCPServerStdio]:
    """Restart an MCP server by creating a new instance.
    
    Parameters
    ----------
    name : str
        Name of the server configuration
    env_vars : Mapping[str, str] | None
        Environment variables for the server
        
    Returns
    -------
    Optional[MCPServerStdio]
        New server instance or None if restart failed
    """
    try:
        logger.info(f"Restarting MCP server '{name}'")
        return setup_mcp_server(name, env_vars, max_retries=2, retry_delay=1.0)
    except Exception as e:
        logger.error(f"Failed to restart server '{name}': {e}")
        return None


def check_mcp_requirements() -> tuple[bool, list[str]]:
    """Check if all required packages for MCP functionality are available.
    
    Returns
    -------
    tuple[bool, list[str]]
        (all_available, missing_packages)
    """
    required_packages = [
        'pydantic_ai',
        'openai',
    ]
    
    missing = []
    
    # Check each required package
    for package in required_packages:
        try:
            __import__(package)
        except ImportError:
            missing.append(package)
    
    # Check for tomllib/tomli
    try:
        import tomllib
    except ImportError:
        try:
            import tomli
        except ImportError:
            missing.append('tomllib or tomli')
    
    return len(missing) == 0, missing


def get_mcp_status_info() -> dict:
    """Get comprehensive status information about MCP setup.
    
    Returns
    -------
    dict
        Status information including packages, API key, and server scripts
    """
    status = {
        'packages_available': True,
        'missing_packages': [],
        'api_key_set': False,
        'servers_toml_exists': False,
        'server_scripts_exist': {},
        'python_version_ok': False
    }
    
    # Check packages
    packages_ok, missing = check_mcp_requirements()
    status['packages_available'] = packages_ok
    status['missing_packages'] = missing
    
    # Check API key
    status['api_key_set'] = bool(os.getenv("OPENAI_API_KEY"))
    
    # Check Python version
    python_version = sys.version_info
    status['python_version_ok'] = python_version >= (3, 9)
    status['python_version'] = f"{python_version.major}.{python_version.minor}.{python_version.micro}"
    
    # Check servers.toml
    servers_toml = Path(__file__).with_name("servers.toml")
    status['servers_toml_exists'] = servers_toml.exists()
    
    # Check server scripts
    base_path = Path(__file__).parent
    server_scripts = {
        'data_extractor': 'mcp_servers/mcp_data_extractor.py',
        'analysis': 'mcp_servers/mcp_analysis.py'
    }
    
    for name, script_path in server_scripts.items():
        full_path = base_path / script_path
        status['server_scripts_exist'][name] = full_path.exists()
    
    return status

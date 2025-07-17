"""
Streamlit App Logging Utilities
==============================

This module provides logging decorators specifically designed for tracing
execution flow in the UORCA Streamlit application. These decorators help
track which functions are being called during app execution.
"""

import logging
import datetime
import pathlib
import inspect
import functools
import asyncio
import os
from typing import Callable, Any

__all__ = [
    "setup_streamlit_logging",
    "log_streamlit_function",
    "log_streamlit_agent",
    "log_streamlit_tab",
]

# Standard log format (not the STREAMLIT format)
STANDARD_FMT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

# Global logger for streamlit app
_streamlit_logger = None
_logging_initialized = False


def setup_streamlit_logging(log_dir: str | os.PathLike = "logs", *, level: int = logging.INFO) -> pathlib.Path:
    """
    Set up Streamlit-specific logging configuration.

    Creates a separate log file for Streamlit app execution tracing.

    Args:
        log_dir: Directory to store log files
        level: Logging level

    Returns:
        Path to the created log file
    """
    global _streamlit_logger, _logging_initialized

    if _streamlit_logger is not None and _logging_initialized:
        # Already initialized - return the expected log file path
        log_dir = pathlib.Path(log_dir)
        return log_dir / "streamlit_app.log"

    log_dir = pathlib.Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / "streamlit_app.log"

    # First, clean up any existing handlers on the root logger that might cause duplicates
    root_logger = logging.getLogger()
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)

    # Configure root logger to use our standard format
    root_handler = logging.FileHandler(log_file, mode='a', encoding="utf-8")
    root_formatter = logging.Formatter(
        fmt=STANDARD_FMT,
        datefmt="%Y-%m-%d %H:%M:%S"
    )
    root_handler.setFormatter(root_formatter)
    root_logger.addHandler(root_handler)
    root_logger.setLevel(level)

    # Create streamlit-specific logger
    _streamlit_logger = logging.getLogger("streamlit_app")
    # Remove existing handlers to avoid duplication
    for h in _streamlit_logger.handlers[:]:
        _streamlit_logger.removeHandler(h)
    _streamlit_logger.setLevel(level)

    # Add the same handler directly to streamlit logger to avoid propagation issues
    streamlit_handler = logging.FileHandler(log_file, mode='a', encoding="utf-8")
    streamlit_handler.setFormatter(root_formatter)
    _streamlit_logger.addHandler(streamlit_handler)

    # Prevent inheriting handlers from root logger to avoid duplicates
    _streamlit_logger.propagate = False

    # Mark as initialized and log only once
    if not _logging_initialized:
        _streamlit_logger.info("Streamlit app logging initialized")
        _logging_initialized = True

    return log_file


def _get_streamlit_logger():
    """Get or create the streamlit logger."""
    global _streamlit_logger, _logging_initialized
    if _streamlit_logger is None or not _logging_initialized:
        setup_streamlit_logging()
    return _streamlit_logger


def _extract_args_for_logging(bound_args, exclude_keys=None):
    """Extract arguments for logging, excluding sensitive or large data."""
    if exclude_keys is None:
        exclude_keys = {'ri', 'ctx', 'results_dir', 'integrator'}

    args_dict = {}
    for k, v in bound_args.arguments.items():
        if k in exclude_keys:
            args_dict[k] = f"<{type(v).__name__}>"
        elif isinstance(v, (list, set, tuple)):
            args_dict[k] = f"<{type(v).__name__}[{len(v)}]>"
        elif isinstance(v, dict):
            args_dict[k] = f"<dict[{len(v)}]>"
        elif isinstance(v, str) and len(v) > 50:
            args_dict[k] = f"{v[:47]}..."
        else:
            args_dict[k] = v

    return args_dict


def log_streamlit_function(func: Callable) -> Callable:
    """
    Decorator for general Streamlit functions to log their execution.

    Use this for tab render functions, helper functions, and UI components.
    """
    logger = _get_streamlit_logger()

    if asyncio.iscoroutinefunction(func):
        @functools.wraps(func)
        async def async_wrapper(*args, **kwargs):
            bound = inspect.signature(func).bind_partial(*args, **kwargs)
            args_dict = _extract_args_for_logging(bound)

            logger.info("%s called", func.__name__)
            try:
                result = await func(*args, **kwargs)
                logger.info("%s completed successfully", func.__name__)
                return result
            except Exception as e:
                logger.error("%s failed: %s", func.__name__, str(e))
                raise
        return async_wrapper

    else:
        @functools.wraps(func)
        def sync_wrapper(*args, **kwargs):
            bound = inspect.signature(func).bind_partial(*args, **kwargs)
            args_dict = _extract_args_for_logging(bound)

            logger.info("%s called", func.__name__)
            try:
                result = func(*args, **kwargs)
                logger.info("%s completed successfully", func.__name__)
                return result
            except Exception as e:
                logger.error("%s failed: %s", func.__name__, str(e))
                raise
        return sync_wrapper


def log_streamlit_agent(func: Callable) -> Callable:
    """
    Decorator for AI agent-related functions in Streamlit.

    Use this for functions that involve AI agent calls, MCP server interactions,
    or other AI-related operations. Provides clear indication of agent activity.
    """
    logger = _get_streamlit_logger()

    if asyncio.iscoroutinefunction(func):
        @functools.wraps(func)
        async def async_wrapper(*args, **kwargs):
            # Don't log detailed args for agent functions (can be very large)
            logger.info("AGENT CALL: %s started", func.__name__)
            try:
                result = await func(*args, **kwargs)
                logger.info("AGENT CALL: %s completed successfully", func.__name__)
                return result
            except Exception as e:
                logger.error("AGENT CALL: %s failed: %s", func.__name__, str(e))
                raise
        return async_wrapper

    else:
        @functools.wraps(func)
        def sync_wrapper(*args, **kwargs):
            # Don't log detailed args for agent functions (can be very large)
            logger.info("AGENT CALL: %s started", func.__name__)
            try:
                result = func(*args, **kwargs)
                logger.info("AGENT CALL: %s completed successfully", func.__name__)
                return result
            except Exception as e:
                logger.error("AGENT CALL: %s failed: %s", func.__name__, str(e))
                raise
        return sync_wrapper


def log_streamlit_tab(tab_name: str):
    """
    Decorator factory for tab render functions.

    Use this to clearly identify which tab is being rendered.

    Example:
        @log_streamlit_tab("Expression Plots")
        def render_expression_plots_tab(...):
            ...
    """
    def decorator(func: Callable) -> Callable:
        logger = _get_streamlit_logger()

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            bound = inspect.signature(func).bind_partial(*args, **kwargs)
            args_dict = _extract_args_for_logging(bound)

            logger.info("TAB RENDER: %s (%s) called", tab_name, func.__name__)
            try:
                result = func(*args, **kwargs)
                logger.info("TAB RENDER: %s (%s) completed", tab_name, func.__name__)
                return result
            except Exception as e:
                logger.error("TAB RENDER: %s (%s) failed: %s", tab_name, func.__name__, str(e))
                raise
        return wrapper
    return decorator


# Convenience functions for manual logging
def log_streamlit_event(message: str, level: int = logging.INFO):
    """Log a general Streamlit app event."""
    logger = _get_streamlit_logger()
    logger.log(level, "EVENT: %s", message)


def log_streamlit_data_load(data_type: str, count: int = None):
    """Log data loading events."""
    logger = _get_streamlit_logger()
    if count is not None:
        logger.info("DATA LOAD: %s (%d items)", data_type, count)
    else:
        logger.info("DATA LOAD: %s", data_type)


def log_streamlit_user_action(action: str, details: str = None):
    """Log user interactions."""
    logger = _get_streamlit_logger()
    if details:
        logger.info("USER ACTION: %s – %s", action, details)
    else:
        logger.info("USER ACTION: %s", action)

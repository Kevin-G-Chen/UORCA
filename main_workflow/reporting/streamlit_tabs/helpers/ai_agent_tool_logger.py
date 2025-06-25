"""
AI Agent Tool Call Logger for UORCA Streamlit App
================================================

File-based logging system for AI agent MCP tool usage in the UORCA Streamlit application.
This replaces the session state approach with temporary JSON log files that persist
across the MCP server process boundary.
"""

import logging
import datetime
import json
import functools
import inspect
import os
from pathlib import Path
from typing import Dict, Any, List, Optional, Union

__all__ = [
    "AIAgentToolLogger",
    "log_ai_agent_tool",
    "get_ai_tool_logger",
    "clear_ai_tool_logs",
    "get_ai_tool_logs_for_display",
    "start_ai_analysis_session",
    "get_current_log_file"
]

# Terminal logging format similar to master agent
logger = logging.getLogger(__name__)

class AIAgentToolLogger:
    """
    File-based tool call logger for AI agent MCP tool usage.

    Stores tool call data in temporary JSON files instead of session state
    to work across process boundaries with MCP servers.
    """

    def __init__(self):
        self.log_file = None
        self.current_analysis_id = None
        self.log_dir = Path("temp_logs")

        # Ensure log directory exists
        self.log_dir.mkdir(exist_ok=True)

    def start_analysis(self, analysis_id: str = None):
        """Start a new analysis session and create a fresh log file."""
        self.current_analysis_id = analysis_id or f"analysis_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}"

        # Use static filename for current analysis
        self.log_file = self.log_dir / "ai_tool_calls_current.json"

        # Delete existing log file if it exists
        if self.log_file.exists():
            self.log_file.unlink()
            logger.info("🧹 Cleared previous AI tool log file")

        # Initialize with empty list
        with open(self.log_file, 'w', encoding='utf-8') as f:
            json.dump([], f, indent=2)

        logger.info("🤖 AI Agent analysis started: %s (log: %s)",
                   self.current_analysis_id, self.log_file)

    def _truncate_output(self, output: Any, max_length: int = 500) -> str:
        """Truncate large outputs to manageable snippets for logging."""
        if output is None:
            return None

        # Convert to string representation
        output_str = str(output)

        if len(output_str) <= max_length:
            return output_str

        # Truncate and add indicator
        return output_str[:max_length] + f"... [truncated, total length: {len(output_str)}]"

    def log_tool_call(self, tool_name: str, parameters: Dict[str, Any],
                     output: Any = None, success: bool = True, error: str = None):
        """
        Log a tool call to the current JSON log file.

        Args:
            tool_name: Name of the MCP tool called
            parameters: Tool parameters (excluding ctx)
            output: Tool output/result
            success: Whether the tool call succeeded
            error: Error message if failed
        """
        if not self.log_file or not self.log_file.exists():
            # If no log file, try to create a default one
            self.start_analysis()

        # Create tool log entry
        tool_log_entry = {
            "tool_name": tool_name,
            "parameters": parameters,
            "timestamp": datetime.datetime.now().isoformat(),
            "success": success,
            "output_snippet": self._truncate_output(output),
            "error": error,
            "analysis_id": self.current_analysis_id
        }

        try:
            # Read current log entries
            with open(self.log_file, 'r', encoding='utf-8') as f:
                log_entries = json.load(f)

            # Add new entry
            log_entries.append(tool_log_entry)

            # Write back to file
            with open(self.log_file, 'w', encoding='utf-8') as f:
                json.dump(log_entries, f, indent=2, ensure_ascii=False)

            # Terminal logging
            if success:
                logger.info("✅ AI Agent tool %s finished with %d chars output",
                           tool_name, len(str(output)) if output else 0)
            else:
                logger.error("❌ AI Agent tool %s failed: %s", tool_name, error)

        except Exception as e:
            logger.error("Failed to log tool call to file %s: %s", self.log_file, e)

    def get_tool_calls_for_display(self) -> List[Dict[str, Any]]:
        """Get all tool calls from the current log file for UI display."""
        if not self.log_file or not self.log_file.exists():
            return []

        try:
            with open(self.log_file, 'r', encoding='utf-8') as f:
                log_entries = json.load(f)

            # Filter to current analysis if specified
            if self.current_analysis_id:
                log_entries = [entry for entry in log_entries
                             if entry.get('analysis_id') == self.current_analysis_id]

            return log_entries

        except Exception as e:
            logger.error("Failed to read tool calls from file %s: %s", self.log_file, e)
            return []

    def get_log_file_path(self) -> Optional[Path]:
        """Get the current log file path."""
        return self.log_file

    def clear_tool_calls(self):
        """Clear all tool calls by creating a fresh log file."""
        if self.log_file and self.log_file.exists():
            with open(self.log_file, 'w', encoding='utf-8') as f:
                json.dump([], f, indent=2)
            logger.info("🧹 AI Agent tool calls cleared from %s", self.log_file)


# Global instance for the Streamlit app
_ai_tool_logger = None

def get_ai_tool_logger() -> AIAgentToolLogger:
    """Get or create the global AI agent tool logger instance."""
    global _ai_tool_logger
    if _ai_tool_logger is None:
        _ai_tool_logger = AIAgentToolLogger()
    return _ai_tool_logger


def log_ai_agent_tool(func):
    """
    Decorator for AI agent MCP tools that provides file-based logging.

    This replaces the session state approach and works across process boundaries.
    Should be applied to the 4 MCP server tools:
    - get_most_common_genes
    - get_gene_contrast_stats
    - filter_genes_by_contrast_sets
    - summarize_contrast
    """
    tool_logger = get_ai_tool_logger()

    @functools.wraps(func)
    async def wrapper(*args, **kwargs):
        # Get tool name and parameters
        tool_name = func.__name__
        bound = inspect.signature(func).bind_partial(*args, **kwargs)

        # Extract parameters (excluding any context-like objects)
        params = {k: v for k, v in bound.arguments.items()
                 if k not in {'ctx', 'self'}}

        # Terminal logging at start
        logger.info("🛠️ AI Agent tool %s called – params=%s", tool_name, params)

        try:
            # Execute the tool
            result = await func(*args, **kwargs)

            # Parse JSON results if they're strings (common for MCP tools)
            parsed_result = result
            if isinstance(result, str):
                try:
                    parsed_result = json.loads(result)
                except json.JSONDecodeError:
                    parsed_result = result

            # Log success
            tool_logger.log_tool_call(
                tool_name=tool_name,
                parameters=params,
                output=parsed_result,
                success=True
            )

            return result

        except Exception as e:
            # Log failure
            error_msg = str(e)
            tool_logger.log_tool_call(
                tool_name=tool_name,
                parameters=params,
                success=False,
                error=error_msg
            )

            # Terminal logging for error
            logger.exception("❌ AI Agent tool %s crashed", tool_name)
            raise

    return wrapper


# Convenience functions for Streamlit components

def clear_ai_tool_logs():
    """Clear AI agent tool logs (convenience function for UI)."""
    tool_logger = get_ai_tool_logger()
    tool_logger.clear_tool_calls()


def get_ai_tool_logs_for_display() -> List[Dict[str, Any]]:
    """Get AI agent tool logs for display (convenience function for UI)."""
    tool_logger = get_ai_tool_logger()
    return tool_logger.get_tool_calls_for_display()


def start_ai_analysis_session(analysis_id: str = None):
    """Start a new AI analysis session (convenience function for UI)."""
    tool_logger = get_ai_tool_logger()
    tool_logger.start_analysis(analysis_id)


def get_current_log_file() -> Optional[Path]:
    """Get the current log file path (convenience function for UI)."""
    tool_logger = get_ai_tool_logger()
    return tool_logger.get_log_file_path()


def read_log_file_contents() -> Optional[str]:
    """Read the raw contents of the current log file."""
    tool_logger = get_ai_tool_logger()
    log_file = tool_logger.get_log_file_path()

    if not log_file or not log_file.exists():
        return None

    try:
        with open(log_file, 'r', encoding='utf-8') as f:
            return f.read()
    except Exception as e:
        logger.error("Failed to read log file contents: %s", e)
        return None

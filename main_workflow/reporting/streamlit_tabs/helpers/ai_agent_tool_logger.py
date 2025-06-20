"""
AI Agent Tool Call Logger for UORCA Streamlit App
================================================

This module provides logging functionality specifically for the AI agent's
MCP tool usage in the UORCA Streamlit application. It mirrors the approach
used in the master agent's workflow_logging.py but is adapted for Streamlit
session state and focuses on the 4 MCP server tools.
"""

import logging
import datetime
import json
import functools
import inspect
import streamlit as st
from typing import Dict, Any, List, Optional

__all__ = [
    "AIAgentToolLogger",
    "log_ai_agent_tool",
    "get_ai_tool_logger",
    "clear_ai_tool_logs",
    "get_ai_tool_logs_for_display"
]

# Terminal logging format similar to master agent
logger = logging.getLogger(__name__)

class AIAgentToolLogger:
    """
    Tool call logger specifically for AI agent MCP tool usage.

    Mirrors the master agent's tool logging approach but stores data
    in Streamlit session state for post-analysis display.
    """

    def __init__(self):
        self.session_key = 'ai_agent_tool_calls'
        self.current_analysis_id = None

    def start_analysis(self, analysis_id: str = None):
        """Start a new analysis session and clear previous tool calls."""
        self.current_analysis_id = analysis_id or f"analysis_{datetime.datetime.now().strftime('%H%M%S')}"

        # Initialize session state if needed
        if self.session_key not in st.session_state:
            st.session_state[self.session_key] = []

        # Clear previous tool calls for new analysis
        st.session_state[self.session_key] = []
        logger.info("🤖 AI Agent analysis started: %s", self.current_analysis_id)

    def log_tool_call(self, tool_name: str, parameters: Dict[str, Any],
                     output: Any = None, success: bool = True, error: str = None):
        """
        Log a tool call with structure similar to master agent's approach.

        Args:
            tool_name: Name of the MCP tool called
            parameters: Tool parameters (excluding ctx)
            output: Tool output/result
            success: Whether the tool call succeeded
            error: Error message if failed
        """
        # Create tool log entry similar to master agent structure
        tool_log_entry = {
            "tool_name": tool_name,
            "parameters": parameters,
            "timestamp": datetime.datetime.now().isoformat(),
            "success": success,
            "output": output,
            "error": error,
            "analysis_id": self.current_analysis_id
        }

        # Initialize session state if needed
        if self.session_key not in st.session_state:
            st.session_state[self.session_key] = []

        # Add to session state
        st.session_state[self.session_key].append(tool_log_entry)

        # Terminal logging similar to master agent
        if success:
            logger.info("✅ AI Agent tool %s finished", tool_name)
        else:
            logger.error("❌ AI Agent tool %s failed: %s", tool_name, error)

    def get_tool_calls_for_display(self) -> List[Dict[str, Any]]:
        """Get all tool calls for the current analysis for UI display."""
        if self.session_key not in st.session_state:
            return []

        # Filter to current analysis if specified
        tool_calls = st.session_state[self.session_key]
        if self.current_analysis_id:
            tool_calls = [tc for tc in tool_calls if tc.get('analysis_id') == self.current_analysis_id]

        return tool_calls

    def clear_tool_calls(self):
        """Clear all tool calls from session state."""
        if self.session_key in st.session_state:
            st.session_state[self.session_key] = []
        logger.info("🧹 AI Agent tool calls cleared")


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
    Decorator for AI agent MCP tools that provides logging similar to master agent's
    log_tool_for_reflection decorator.

    This should be applied to the 4 MCP server tools:
    - get_most_common_genes
    - get_gene_contrast_stats
    - filter_genes_by_contrast_sets
    - summarize_contrast
    """
    tool_logger = get_ai_tool_logger()

    @functools.wraps(func)
    async def wrapper(*args, **kwargs):
        # Get tool name and parameters (similar to master agent approach)
        tool_name = func.__name__
        bound = inspect.signature(func).bind_partial(*args, **kwargs)

        # Extract parameters (excluding any context-like objects)
        params = {k: v for k, v in bound.arguments.items()
                 if k not in {'ctx', 'self'}}

        # Terminal logging at start (similar to master agent)
        logger.info("🛠️ AI Agent tool %s called – params=%s", tool_name, params)

        tool_log_entry = {
            "tool_name": tool_name,
            "parameters": params,
            "timestamp": datetime.datetime.now().isoformat(),
            "success": None,
            "output": None,
            "error": None
        }

        try:
            # Execute the tool
            result = await func(*args, **kwargs)

            # Handle JSON string outputs (common for MCP tools)
            if isinstance(result, str):
                try:
                    # Try to parse JSON for better display
                    parsed_result = json.loads(result)
                    tool_log_entry["output"] = parsed_result
                except json.JSONDecodeError:
                    tool_log_entry["output"] = result
            else:
                tool_log_entry["output"] = result

            tool_log_entry["success"] = True

            # Log success
            tool_logger.log_tool_call(
                tool_name=tool_name,
                parameters=params,
                output=tool_log_entry["output"],
                success=True
            )

            return result

        except Exception as e:
            # Log failure
            error_msg = str(e)
            tool_log_entry["success"] = False
            tool_log_entry["error"] = error_msg

            tool_logger.log_tool_call(
                tool_name=tool_name,
                parameters=params,
                success=False,
                error=error_msg
            )

            # Terminal logging for error (similar to master agent)
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

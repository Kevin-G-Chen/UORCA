"""
Tool Relevance Analyzer for UORCA AI Assistant
==============================================

Analyzes which AI tool calls contributed to the final interpretation using individual
LLM calls for each tool. Creates a filtered log file containing only relevant tools
in the same format as the original tool logs.
"""

import logging
import asyncio
import json
from typing import List, Dict, Any, Optional
from pathlib import Path
import streamlit as st

from pydantic import BaseModel, Field
from pydantic_ai import Agent
from config_loader import get_ai_agent_config

logger = logging.getLogger(__name__)


class ToolRelevanceResult(BaseModel):
    """Simple boolean result for tool relevance."""
    is_relevant: bool = Field(..., description="Whether this tool supports the findings (True/False)")
    brief_reason: str = Field(..., min_length=5, max_length=200, description="Brief explanation for the decision")


def create_tool_relevance_agent() -> Optional[Agent]:
    """Create a simple agent for analyzing individual tool relevance."""
    try:
        import os
        if not os.getenv("OPENAI_API_KEY"):
            logger.warning("OpenAI API key not found - tool relevance analysis skipped")
            return None

        ai_config = get_ai_agent_config()

        system_prompt = """You are evaluating whether a specific AI tool call contributed to a final research interpretation.

Your task is simple: determine if the tool's output was used to support the final interpretation.

RELEVANT tools:
- Provided data directly mentioned in the final interpretation
- Generated results that supported the conclusions
- Identified genes, patterns, or insights referenced in the final text
- Returned useful data that informed the final analysis

NOT RELEVANT tools:
- Were exploratory but yielded no useful information for the final conclusion
- Returned empty results, errors, or irrelevant data
- Provided information that wasn't used in the final interpretation
- Were redundant or superseded by better tool calls

Be strict: only mark as relevant if the tool's output clearly supports the final interpretation."""

        agent = Agent(
            model=ai_config.model,
            model_settings={"temperature": 0.1},  # Very low temperature for consistent decisions
            system_prompt=system_prompt,
            output_type=ToolRelevanceResult,
        )

        return agent

    except Exception as e:
        logger.error(f"Failed to create tool relevance agent: {e}")
        return None


def format_single_tool_for_analysis(tool_call: Dict[str, Any], index: int) -> str:
    """Format a single tool call for relevance analysis."""
    formatted = f"TOOL CALL #{index}:\n"
    formatted += f"Tool Name: {tool_call.get('tool_name', 'Unknown')}\n"
    formatted += f"Success: {tool_call.get('success', True)}\n"

    # Parameters
    params = tool_call.get('parameters', {})
    if params:
        formatted += f"Parameters: {params}\n"

    # Output
    if tool_call.get('success', True):
        output = tool_call.get('full_output') or tool_call.get('output_snippet')
        if output:
            # Truncate very long outputs
            output_str = str(output)
            if len(output_str) > 1500:
                output_str = output_str[:1500] + f"... [truncated]"
            formatted += f"Output: {output_str}\n"
        else:
            formatted += "Output: No output returned\n"
    else:
        error = tool_call.get('error', 'Unknown error')
        formatted += f"Error: {error}\n"

    return formatted


async def analyze_single_tool_relevance(
    agent: Agent,
    tool_call: Dict[str, Any],
    index: int,
    final_interpretation: str
) -> tuple[int, bool, str]:
    """
    Analyze relevance of a single tool call.

    Returns:
        tuple of (index, is_relevant, reason)
    """
    try:
        tool_text = format_single_tool_for_analysis(tool_call, index)

        analysis_input = f"""FINAL INTERPRETATION:
{final_interpretation}

{tool_text}

Question: Does this tool call support or contribute to the final interpretation above?
Answer with True if the tool's output was used in reaching the conclusion, False otherwise."""

        result = await agent.run(analysis_input)

        if hasattr(result, 'data'):
            relevance_result = result.data
            return (index, relevance_result.is_relevant, relevance_result.brief_reason)
        else:
            logger.warning(f"Unexpected result format for tool {index}")
            return (index, False, "Analysis failed")

    except Exception as e:
        logger.error(f"Error analyzing tool {index}: {e}")
        return (index, False, f"Error: {str(e)}")


async def analyze_all_tools_parallel(
    tool_calls: List[Dict[str, Any]],
    final_interpretation: str
) -> Dict[int, tuple[bool, str]]:
    """
    Analyze all tool calls in parallel.

    Returns:
        Dict mapping tool index to (is_relevant, reason)
    """
    if not tool_calls:
        return {}

    agent = create_tool_relevance_agent()
    if not agent:
        logger.warning("Tool relevance agent not available")
        return {}

    logger.info(f"Starting parallel analysis of {len(tool_calls)} tools")

    # Create tasks for all tools
    tasks = []
    for i, tool_call in enumerate(tool_calls):
        task = analyze_single_tool_relevance(agent, tool_call, i, final_interpretation)
        tasks.append(task)

    # Run all analyses in parallel
    try:
        results = await asyncio.gather(*tasks, return_exceptions=True)

        # Process results
        relevance_map = {}
        successful_analyses = 0

        for result in results:
            if isinstance(result, Exception):
                logger.error(f"Tool analysis failed: {result}")
                continue

            index, is_relevant, reason = result
            relevance_map[index] = (is_relevant, reason)
            if is_relevant:
                successful_analyses += 1

        logger.info(f"Tool relevance analysis completed: {successful_analyses}/{len(tool_calls)} tools deemed relevant")
        return relevance_map

    except Exception as e:
        logger.error(f"Error in parallel tool analysis: {e}")
        return {}


def create_relevant_tool_log_file(
    original_tool_calls: List[Dict[str, Any]],
    relevance_map: Dict[int, tuple[bool, str]]
) -> Optional[Path]:
    """
    Create a new log file containing only relevant tool calls.

    Args:
        original_tool_calls: Original tool call logs
        relevance_map: Map of tool index to (is_relevant, reason)

    Returns:
        Path to the relevant tool log file
    """
    try:
        # Use same directory structure as original logger
        from streamlit_tabs.helpers.ai_agent_tool_logger import get_ai_tool_logger

        tool_logger = get_ai_tool_logger()
        log_dir = tool_logger.log_dir

        # Create relevant tools log file
        relevant_log_file = log_dir / "ai_tool_calls_relevant.json"

        # Filter to only relevant tools
        relevant_tool_calls = []
        for i, tool_call in enumerate(original_tool_calls):
            if i in relevance_map and relevance_map[i][0]:  # is_relevant is True
                # Add the relevance reason as metadata but keep original format
                enhanced_call = tool_call.copy()
                enhanced_call['relevance_reason'] = relevance_map[i][1]
                relevant_tool_calls.append(enhanced_call)

        # Write filtered calls to new file
        with open(relevant_log_file, 'w', encoding='utf-8') as f:
            json.dump(relevant_tool_calls, f, indent=2, ensure_ascii=False)

        logger.info(f"Created relevant tool log file with {len(relevant_tool_calls)} tools: {relevant_log_file}")
        return relevant_log_file

    except Exception as e:
        logger.error(f"Failed to create relevant tool log file: {e}")
        return None


def run_tool_relevance_analysis_sync(
    tool_calls: List[Dict[str, Any]],
    final_interpretation: str
) -> tuple[Optional[Path], Dict[int, tuple[bool, str]]]:
    """
    Synchronous wrapper for tool relevance analysis.

    Args:
        tool_calls: Original tool call logs
        final_interpretation: Final AI interpretation text

    Returns:
        tuple of (relevant_log_file_path, relevance_map)
    """
    if not tool_calls:
        logger.info("No tool calls to analyze")
        return None, {}

    try:
        # Run async analysis
        loop = asyncio.new_event_loop()
        asyncio.set_event_loop(loop)

        relevance_map = loop.run_until_complete(
            analyze_all_tools_parallel(tool_calls, final_interpretation)
        )

        # Create filtered log file
        relevant_log_file = create_relevant_tool_log_file(tool_calls, relevance_map)

        return relevant_log_file, relevance_map

    except Exception as e:
        logger.error(f"Error in synchronous tool relevance analysis: {e}")
        return None, {}
    finally:
        try:
            loop.close()
        except:
            pass


def get_relevant_tool_calls() -> List[Dict[str, Any]]:
    """
    Load relevant tool calls from the filtered log file.

    Returns:
        List of relevant tool calls
    """
    try:
        from streamlit_tabs.helpers.ai_agent_tool_logger import get_ai_tool_logger

        tool_logger = get_ai_tool_logger()
        log_dir = tool_logger.log_dir
        relevant_log_file = log_dir / "ai_tool_calls_relevant.json"

        if not relevant_log_file.exists():
            return []

        with open(relevant_log_file, 'r', encoding='utf-8') as f:
            return json.load(f)

    except Exception as e:
        logger.error(f"Failed to load relevant tool calls: {e}")
        return []


def clear_relevant_tool_log():
    """Clear the relevant tool log file."""
    try:
        from streamlit_tabs.helpers.ai_agent_tool_logger import get_ai_tool_logger

        tool_logger = get_ai_tool_logger()
        log_dir = tool_logger.log_dir
        relevant_log_file = log_dir / "ai_tool_calls_relevant.json"

        if relevant_log_file.exists():
            relevant_log_file.unlink()
            logger.info("Cleared relevant tool log file")

    except Exception as e:
        logger.error(f"Failed to clear relevant tool log: {e}")

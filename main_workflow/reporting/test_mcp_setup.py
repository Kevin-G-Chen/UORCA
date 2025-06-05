#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Test script to verify MCP server setup for UORCA RNA-seq analysis.

This script tests:
1. MCP server initialization
2. Tool availability
3. Basic functionality of key tools
4. Agent integration
"""

import os
import sys
import asyncio
import json
from pathlib import Path

# Add parent directory to path for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'scratch', 'MCP_examples'))

from rich.console import Console
from rich.table import Table
from rich import print as rprint

console = Console()


async def test_mcp_servers():
    """Test MCP server functionality."""
    console.rule("Testing MCP Server Setup")
    
    # Test environment setup
    results_dir = os.environ.get("UORCA_RESULTS", "")
    if not results_dir:
        # Try to find a test results directory
        test_dirs = [
            "./UORCA_results",
            "../UORCA_results",
            "../../UORCA_results",
            os.path.join(os.path.dirname(__file__), "..", "..", "UORCA_results")
        ]
        
        for test_dir in test_dirs:
            if os.path.exists(test_dir):
                results_dir = os.path.abspath(test_dir)
                break
        
        if not results_dir:
            console.print("[red]❌ No UORCA_RESULTS directory found. Please set UORCA_RESULTS environment variable.[/red]")
            return False
    
    console.print(f"[green]✓[/green] Using results directory: {results_dir}")
    
    # Test importing required modules
    try:
        from mcp_utils import setup_mcp_server
        console.print("[green]✓[/green] Successfully imported mcp_utils")
    except Exception as e:
        console.print(f"[red]❌ Failed to import mcp_utils: {e}[/red]")
        return False
    
    try:
        from reporting_agent import ReportingAgent
        console.print("[green]✓[/green] Successfully imported ReportingAgent")
    except Exception as e:
        console.print(f"[red]❌ Failed to import ReportingAgent: {e}[/red]")
        return False
    
    # Test MCP server setup
    console.rule("Testing MCP Servers")
    
    servers_to_test = ["data_extractor", "analysis"]
    working_servers = []
    
    for server_name in servers_to_test:
        try:
            console.print(f"\nTesting {server_name} server...")
            server = setup_mcp_server(server_name, env_vars={"RESULTS_DIR": results_dir})
            console.print(f"[green]✓[/green] Successfully created {server_name} server")
            working_servers.append(server)
            
            # Cleanup
            server.cleanup()
            
        except Exception as e:
            console.print(f"[red]❌ Failed to create {server_name} server: {e}[/red]")
            return False
    
    # Test agent creation
    console.rule("Testing Agent Integration")
    
    try:
        agent = ReportingAgent(results_dir=results_dir)
        console.print("[green]✓[/green] Successfully created ReportingAgent")
        
        # Test server setup
        servers = agent.setup_servers()
        console.print(f"[green]✓[/green] Successfully set up {len(servers)} MCP servers")
        
        # Test agent creation
        test_agent = agent.create_agent("Test research question")
        console.print("[green]✓[/green] Successfully created pydantic_ai Agent")
        
        # Cleanup
        agent.cleanup()
        console.print("[green]✓[/green] Successfully cleaned up resources")
        
    except Exception as e:
        console.print(f"[red]❌ Failed agent integration test: {e}[/red]")
        return False
    
    return True


async def test_basic_tools():
    """Test basic tool functionality."""
    console.rule("Testing Basic Tool Functionality")
    
    results_dir = os.environ.get("UORCA_RESULTS", "./UORCA_results")
    if not os.path.exists(results_dir):
        console.print("[yellow]⚠️  Skipping tool tests - no data directory available[/yellow]")
        return True
    
    try:
        from pydantic_ai import Agent
        from mcp_utils import setup_mcp_server
        
        # Setup data extractor server
        data_server = setup_mcp_server("data_extractor", env_vars={"RESULTS_DIR": results_dir})
        
        # Create a simple test agent
        agent = Agent(
            model="openai:gpt-4o-mini",
            model_settings={"temperature": 0.1},
            mcp_servers=[data_server],
            system="You are a test bot. Use the list_available_contrasts tool and return the number of contrasts found."
        )
        
        console.print("\nTesting list_available_contrasts tool...")
        
        # Test the tool
        result = await agent.run("Please list all available contrasts and tell me how many there are.")
        
        console.print(f"[green]✓[/green] Tool execution successful")
        console.print(f"Result preview: {str(result.data)[:200]}...")
        
        # Cleanup
        data_server.cleanup()
        
    except Exception as e:
        console.print(f"[red]❌ Tool test failed: {e}[/red]")
        return False
    
    return True


async def test_auto_start():
    """Test auto-start manager."""
    console.rule("Testing Auto-Start Manager")
    
    results_dir = os.environ.get("UORCA_RESULTS", "./UORCA_results")
    if not os.path.exists(results_dir):
        console.print("[yellow]⚠️  Skipping auto-start test - no data directory available[/yellow]")
        return True
    
    try:
        from auto_start_manager import AutoStartManager
        
        manager = AutoStartManager(results_dir)
        console.print("[green]✓[/green] Successfully created AutoStartManager")
        
        # Test validation
        is_valid, error = manager.validate_data_directory()
        if is_valid:
            console.print("[green]✓[/green] Data directory validation passed")
        else:
            console.print(f"[yellow]⚠️  Data directory validation failed: {error}[/yellow]")
        
        # Cleanup
        manager.cleanup()
        
    except Exception as e:
        console.print(f"[red]❌ Auto-start test failed: {e}[/red]")
        return False
    
    return True


def create_summary_table(test_results):
    """Create a summary table of test results."""
    table = Table(title="MCP Setup Test Summary")
    table.add_column("Test", style="cyan")
    table.add_column("Status", style="green")
    table.add_column("Details")
    
    for test_name, (passed, details) in test_results.items():
        status = "[green]✓ PASSED[/green]" if passed else "[red]✗ FAILED[/red]"
        table.add_row(test_name, status, details)
    
    return table


async def main():
    """Run all tests."""
    console.print("\n[bold blue]UORCA MCP Setup Test Suite[/bold blue]\n")
    
    test_results = {}
    
    # Run tests
    console.print("Running tests...\n")
    
    # Test 1: MCP Servers
    try:
        passed = await test_mcp_servers()
        test_results["MCP Servers"] = (passed, "Server initialization and agent integration")
    except Exception as e:
        test_results["MCP Servers"] = (False, f"Error: {str(e)}")
    
    # Test 2: Basic Tools
    try:
        passed = await test_basic_tools()
        test_results["Basic Tools"] = (passed, "Tool availability and execution")
    except Exception as e:
        test_results["Basic Tools"] = (False, f"Error: {str(e)}")
    
    # Test 3: Auto-Start
    try:
        passed = await test_auto_start()
        test_results["Auto-Start"] = (passed, "Auto-start manager functionality")
    except Exception as e:
        test_results["Auto-Start"] = (False, f"Error: {str(e)}")
    
    # Display summary
    console.print("\n")
    console.print(create_summary_table(test_results))
    
    # Overall result
    all_passed = all(result[0] for result in test_results.values())
    
    console.print("\n")
    if all_passed:
        console.print("[bold green]✅ All tests passed! MCP setup is ready for use.[/bold green]")
    else:
        console.print("[bold red]❌ Some tests failed. Please check the errors above.[/bold red]")
    
    return all_passed


if __name__ == "__main__":
    # Run the test suite
    success = asyncio.run(main())
    sys.exit(0 if success else 1)
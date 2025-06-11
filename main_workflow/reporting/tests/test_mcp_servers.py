#!/usr/bin/env python3
"""
Test script for MCP server functionality in UORCA

This script tests the MCP server setup, health checks, and basic functionality
to help diagnose server management issues.

Usage:
    python test_mcp_servers.py [--results-dir /path/to/results]
"""

import os
import sys
import asyncio
import argparse
import logging
import time
from pathlib import Path

# Add parent directory to path for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    from mcp_utils import setup_mcp_server, verify_server_responsive, check_server_health, is_server_process_alive
    from auto_start_manager import AutoStartManager
    from reporting_agent import ReportingAgent
except ImportError as e:
    print(f"Error importing modules: {e}")
    print("Make sure you're running this from the correct directory")
    sys.exit(1)

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger("mcp_test")

class Colors:
    """ANSI color codes for terminal output"""
    GREEN = '\033[92m'
    RED = '\033[91m'
    YELLOW = '\033[93m'
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    WHITE = '\033[97m'
    BOLD = '\033[1m'
    END = '\033[0m'

def print_test_header(test_name):
    """Print a formatted test header"""
    print(f"\n{Colors.CYAN}{Colors.BOLD}{'='*60}{Colors.END}")
    print(f"{Colors.CYAN}{Colors.BOLD}Testing: {test_name}{Colors.END}")
    print(f"{Colors.CYAN}{Colors.BOLD}{'='*60}{Colors.END}")

def print_success(message):
    """Print a success message"""
    print(f"{Colors.GREEN}✅ {message}{Colors.END}")

def print_error(message):
    """Print an error message"""
    print(f"{Colors.RED}❌ {message}{Colors.END}")

def print_warning(message):
    """Print a warning message"""
    print(f"{Colors.YELLOW}⚠️  {message}{Colors.END}")

def print_info(message):
    """Print an info message"""
    print(f"{Colors.BLUE}ℹ️  {message}{Colors.END}")

def test_environment():
    """Test the environment setup"""
    print_test_header("Environment Setup")
    
    # Check Python version
    python_version = sys.version_info
    if python_version >= (3, 9):
        print_success(f"Python version: {python_version.major}.{python_version.minor}.{python_version.micro}")
    else:
        print_error(f"Python version too old: {python_version.major}.{python_version.minor}.{python_version.micro} (need >= 3.9)")
        return False
    
    # Check required modules
    required_modules = ['pydantic_ai', 'asyncio', 'tomllib']
    try:
        import tomllib
    except ImportError:
        try:
            import tomli as tomllib
            required_modules.append('tomli')
        except ImportError:
            print_error("Neither tomllib nor tomli available")
            return False
    
    for module in required_modules:
        try:
            __import__(module)
            print_success(f"Module {module} available")
        except ImportError:
            print_error(f"Module {module} not available")
            return False
    
    # Check OpenAI API key
    api_key = os.getenv("OPENAI_API_KEY")
    if api_key:
        print_success("OpenAI API key is set")
    else:
        print_warning("OpenAI API key not set (some features may not work)")
    
    # Check servers.toml
    servers_toml = Path(__file__).parent.parent / "servers.toml"
    if servers_toml.exists():
        print_success(f"servers.toml found at {servers_toml}")
    else:
        print_error(f"servers.toml not found at {servers_toml}")
        return False
    
    return True

def test_server_scripts():
    """Test that server scripts exist and are executable"""
    print_test_header("Server Script Availability")
    
    base_path = Path(__file__).parent.parent
    server_scripts = [
        "mcp_servers/mcp_data_extractor.py",
        "mcp_servers/mcp_analysis.py"
    ]
    
    all_good = True
    for script in server_scripts:
        script_path = base_path / script
        if script_path.exists():
            print_success(f"Server script found: {script}")
            # Check if it's executable
            if os.access(script_path, os.X_OK):
                print_success(f"  Script is executable")
            else:
                print_warning(f"  Script is not executable (may still work)")
        else:
            print_error(f"Server script missing: {script}")
            all_good = False
    
    return all_good

def test_single_server_startup(server_name, results_dir=None):
    """Test starting a single MCP server"""
    print_info(f"Testing {server_name} server startup...")
    
    env_vars = {}
    if results_dir:
        env_vars["RESULTS_DIR"] = results_dir
    
    api_key = os.getenv("OPENAI_API_KEY")
    if api_key:
        env_vars["OPENAI_API_KEY"] = api_key
    
    try:
        server = setup_mcp_server(server_name, env_vars, max_retries=2, retry_delay=1.0)
        print_success(f"  {server_name} server started successfully")
        
        # Test process health
        if is_server_process_alive(server):
            print_success(f"  {server_name} server process is alive")
        else:
            print_error(f"  {server_name} server process is not alive")
            return False, None
        
        # Test responsiveness
        if verify_server_responsive(server, timeout=5.0):
            print_success(f"  {server_name} server is responsive")
        else:
            print_error(f"  {server_name} server is not responsive")
            return False, server
        
        return True, server
        
    except Exception as e:
        print_error(f"  Failed to start {server_name} server: {e}")
        return False, None

async def test_server_communication(server, server_name):
    """Test basic communication with an MCP server"""
    print_info(f"Testing {server_name} server communication...")
    
    try:
        # Test listing tools
        tools = await asyncio.wait_for(server.list_tools(), timeout=10.0)
        print_success(f"  {server_name} server has {len(tools)} tools available")
        
        # Print tool names
        tool_names = [tool.name for tool in tools]
        print_info(f"  Available tools: {', '.join(tool_names)}")
        
        return True
        
    except asyncio.TimeoutError:
        print_error(f"  {server_name} server communication timed out")
        return False
    except Exception as e:
        print_error(f"  {server_name} server communication failed: {e}")
        return False

async def test_server_cleanup(server, server_name):
    """Test server cleanup"""
    print_info(f"Testing {server_name} server cleanup...")
    
    try:
        await server.cleanup()
        print_success(f"  {server_name} server cleaned up successfully")
        
        # Wait a moment and check if process is gone
        await asyncio.sleep(1)
        if not is_server_process_alive(server):
            print_success(f"  {server_name} server process terminated")
        else:
            print_warning(f"  {server_name} server process still alive after cleanup")
        
        return True
        
    except Exception as e:
        print_error(f"  {server_name} server cleanup failed: {e}")
        return False

def test_auto_start_manager(results_dir=None):
    """Test the AutoStartManager functionality"""
    print_test_header("AutoStartManager Testing")
    
    if not results_dir:
        print_warning("No results directory provided, skipping AutoStartManager test")
        return False
    
    try:
        manager = AutoStartManager(results_dir)
        
        # Test directory validation
        is_valid, error = manager.validate_data_directory()
        if is_valid:
            print_success("Results directory validation passed")
        else:
            print_error(f"Results directory validation failed: {error}")
            return False
        
        # Test server setup
        manager.setup_servers()
        print_success(f"AutoStartManager set up {len(manager.servers)} servers")
        
        # Test server health
        if manager.verify_all_servers_healthy():
            print_success("All servers are healthy")
        else:
            print_error("Some servers are not healthy")
        
        # Cleanup
        asyncio.run(manager.cleanup())
        print_success("AutoStartManager cleanup completed")
        
        return True
        
    except Exception as e:
        print_error(f"AutoStartManager test failed: {e}")
        return False

def test_reporting_agent(results_dir=None):
    """Test the ReportingAgent functionality"""
    print_test_header("ReportingAgent Testing")
    
    if not results_dir:
        print_warning("No results directory provided, skipping ReportingAgent test")
        return False
    
    try:
        agent = ReportingAgent(results_dir)
        
        # Test server setup
        agent.setup_servers()
        print_success(f"ReportingAgent set up {len(agent.servers)} servers")
        
        # Test server health
        if agent.verify_all_servers_healthy():
            print_success("All servers are healthy")
        else:
            print_error("Some servers are not healthy")
        
        # Test agent creation (without running analysis)
        agent.create_agent("Test biological context")
        print_success("ReportingAgent created successfully")
        
        # Cleanup
        asyncio.run(agent.cleanup())
        print_success("ReportingAgent cleanup completed")
        
        return True
        
    except Exception as e:
        print_error(f"ReportingAgent test failed: {e}")
        return False

async def run_all_tests(results_dir=None):
    """Run all tests"""
    print(f"{Colors.BOLD}{Colors.WHITE}MCP Server Functionality Test Suite{Colors.END}")
    print(f"{Colors.WHITE}{'='*60}{Colors.END}")
    
    test_results = {}
    
    # Test 1: Environment
    test_results['environment'] = test_environment()
    
    # Test 2: Server scripts
    test_results['scripts'] = test_server_scripts()
    
    # Test 3: Individual server startup
    print_test_header("Individual Server Testing")
    
    servers_to_test = ['data_extractor', 'analysis']
    test_results['servers'] = {}
    
    for server_name in servers_to_test:
        success, server = test_single_server_startup(server_name, results_dir)
        test_results['servers'][server_name] = success
        
        if success and server:
            # Test communication
            comm_success = await test_server_communication(server, server_name)
            test_results['servers'][f'{server_name}_comm'] = comm_success
            
            # Test cleanup
            cleanup_success = await test_server_cleanup(server, server_name)
            test_results['servers'][f'{server_name}_cleanup'] = cleanup_success
    
    # Test 4: AutoStartManager
    test_results['auto_start_manager'] = test_auto_start_manager(results_dir)
    
    # Test 5: ReportingAgent
    test_results['reporting_agent'] = test_reporting_agent(results_dir)
    
    # Print summary
    print_test_header("Test Summary")
    
    total_tests = 0
    passed_tests = 0
    
    def count_results(results, prefix=""):
        nonlocal total_tests, passed_tests
        for key, value in results.items():
            if isinstance(value, dict):
                count_results(value, f"{prefix}{key}.")
            else:
                total_tests += 1
                if value:
                    passed_tests += 1
                    print_success(f"{prefix}{key}")
                else:
                    print_error(f"{prefix}{key}")
    
    count_results(test_results)
    
    print(f"\n{Colors.BOLD}Overall Results: {passed_tests}/{total_tests} tests passed{Colors.END}")
    
    if passed_tests == total_tests:
        print_success("All tests passed! MCP servers should work correctly.")
        return True
    else:
        print_error(f"{total_tests - passed_tests} test(s) failed. Check the output above for details.")
        return False

def main():
    """Main function"""
    parser = argparse.ArgumentParser(description="Test MCP server functionality")
    parser.add_argument(
        "--results-dir",
        help="Path to UORCA results directory (required for some tests)"
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable verbose logging"
    )
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Run tests
    success = asyncio.run(run_all_tests(args.results_dir))
    
    if not success:
        print(f"\n{Colors.YELLOW}Troubleshooting Tips:{Colors.END}")
        print("1. Make sure you're running this from the reporting directory")
        print("2. Check that all required Python packages are installed")
        print("3. Verify that server scripts exist in mcp_servers/ directory")
        print("4. Set OPENAI_API_KEY environment variable if using AI features")
        print("5. Provide --results-dir argument for full testing")
        
        sys.exit(1)
    else:
        print(f"\n{Colors.GREEN}✅ All tests passed! MCP server system is working correctly.{Colors.END}")

if __name__ == "__main__":
    main()
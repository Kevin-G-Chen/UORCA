#!/usr/bin/env python3
"""
MCP Server Recovery Utility for UORCA

This script helps diagnose and fix common MCP server issues in the UORCA
RNA-seq analysis system. It can detect stuck processes, clean up orphaned
servers, and restart the MCP server infrastructure.

Usage:
    python recover_mcp_servers.py [options]

Options:
    --kill-all          Force kill all MCP server processes
    --restart           Restart MCP servers after cleanup
    --check-only        Only check status, don't make changes
    --results-dir PATH  Specify results directory for testing
    --verbose           Enable detailed logging
"""

import os
import sys
import signal
import psutil
import argparse
import logging
import asyncio
import time
from pathlib import Path
from typing import List, Dict, Optional

# Add parent directory to path for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger("mcp_recovery")

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

def print_header(message):
    """Print a formatted header"""
    print(f"\n{Colors.CYAN}{Colors.BOLD}{'='*60}{Colors.END}")
    print(f"{Colors.CYAN}{Colors.BOLD}{message}{Colors.END}")
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

def find_mcp_processes() -> List[Dict]:
    """Find all running MCP server processes"""
    mcp_processes = []
    
    try:
        for proc in psutil.process_iter(['pid', 'name', 'cmdline', 'create_time', 'status']):
            try:
                cmdline = proc.info['cmdline']
                if cmdline and any('mcp_' in arg for arg in cmdline):
                    # Check if it's one of our MCP servers
                    if any(server in ' '.join(cmdline) for server in ['mcp_data_extractor', 'mcp_analysis']):
                        mcp_processes.append({
                            'pid': proc.info['pid'],
                            'name': proc.info['name'],
                            'cmdline': cmdline,
                            'create_time': proc.info['create_time'],
                            'status': proc.info['status'],
                            'process': proc
                        })
            except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
                continue
                
    except Exception as e:
        print_error(f"Error scanning for processes: {e}")
    
    return mcp_processes

def find_python_processes_with_mcp() -> List[Dict]:
    """Find Python processes that might be running MCP servers"""
    python_processes = []
    
    try:
        for proc in psutil.process_iter(['pid', 'name', 'cmdline', 'create_time', 'status']):
            try:
                if proc.info['name'] and 'python' in proc.info['name'].lower():
                    cmdline = proc.info['cmdline']
                    if cmdline and any('mcp' in arg.lower() for arg in cmdline):
                        python_processes.append({
                            'pid': proc.info['pid'],
                            'name': proc.info['name'],
                            'cmdline': cmdline,
                            'create_time': proc.info['create_time'],
                            'status': proc.info['status'],
                            'process': proc
                        })
            except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
                continue
                
    except Exception as e:
        print_error(f"Error scanning for Python processes: {e}")
    
    return python_processes

def check_server_status():
    """Check the current status of MCP servers"""
    print_header("MCP Server Status Check")
    
    # Find MCP processes
    mcp_processes = find_mcp_processes()
    python_mcp_processes = find_python_processes_with_mcp()
    
    all_processes = mcp_processes + python_mcp_processes
    
    if not all_processes:
        print_info("No MCP server processes found")
        return []
    
    print_info(f"Found {len(all_processes)} MCP-related processes:")
    
    for proc_info in all_processes:
        pid = proc_info['pid']
        status = proc_info['status']
        create_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(proc_info['create_time']))
        cmdline = ' '.join(proc_info['cmdline'][:3])  # First 3 args
        
        if status == 'zombie':
            print_error(f"  PID {pid}: ZOMBIE - {cmdline} (created: {create_time})")
        elif status == 'sleeping':
            print_warning(f"  PID {pid}: SLEEPING - {cmdline} (created: {create_time})")
        else:
            print_success(f"  PID {pid}: {status.upper()} - {cmdline} (created: {create_time})")
    
    return all_processes

def kill_process_safely(proc_info: Dict, force: bool = False) -> bool:
    """Safely kill a process"""
    pid = proc_info['pid']
    process = proc_info['process']
    
    try:
        if not process.is_running():
            print_info(f"Process {pid} is already dead")
            return True
        
        if force:
            print_warning(f"Force killing process {pid}")
            process.kill()
        else:
            print_info(f"Terminating process {pid} gracefully")
            process.terminate()
        
        # Wait for process to die
        try:
            process.wait(timeout=5)
            print_success(f"Process {pid} terminated successfully")
            return True
        except psutil.TimeoutExpired:
            if not force:
                print_warning(f"Process {pid} didn't terminate gracefully, force killing")
                process.kill()
                process.wait(timeout=5)
                print_success(f"Process {pid} force killed")
                return True
            else:
                print_error(f"Failed to kill process {pid}")
                return False
                
    except psutil.NoSuchProcess:
        print_info(f"Process {pid} is already gone")
        return True
    except psutil.AccessDenied:
        print_error(f"Access denied when trying to kill process {pid}")
        return False
    except Exception as e:
        print_error(f"Error killing process {pid}: {e}")
        return False

def cleanup_mcp_processes(force: bool = False) -> bool:
    """Clean up all MCP server processes"""
    print_header("Cleaning Up MCP Server Processes")
    
    all_processes = find_mcp_processes() + find_python_processes_with_mcp()
    
    if not all_processes:
        print_success("No MCP processes to clean up")
        return True
    
    print_info(f"Found {len(all_processes)} processes to clean up")
    
    success_count = 0
    for proc_info in all_processes:
        if kill_process_safely(proc_info, force):
            success_count += 1
    
    if success_count == len(all_processes):
        print_success(f"Successfully cleaned up all {success_count} processes")
        return True
    else:
        print_error(f"Only {success_count}/{len(all_processes)} processes cleaned up successfully")
        return False

def test_server_startup(results_dir: Optional[str] = None) -> bool:
    """Test starting MCP servers"""
    print_header("Testing MCP Server Startup")
    
    try:
        # Import here to avoid import errors if modules aren't available
        from auto_start_manager import AutoStartManager
        
        if not results_dir:
            print_warning("No results directory provided, using dummy directory")
            results_dir = "/tmp/dummy_results"
            os.makedirs(results_dir, exist_ok=True)
        
        print_info(f"Using results directory: {results_dir}")
        
        manager = AutoStartManager(results_dir)
        
        # Test server setup
        print_info("Setting up MCP servers...")
        manager.setup_servers()
        
        if len(manager.servers) > 0:
            print_success(f"Successfully started {len(manager.servers)} servers")
            
            # Test server health
            if manager.verify_all_servers_healthy():
                print_success("All servers are healthy")
                
                # Cleanup
                print_info("Cleaning up test servers...")
                asyncio.run(manager.cleanup())
                print_success("Test servers cleaned up")
                return True
            else:
                print_error("Some servers failed health check")
                asyncio.run(manager.cleanup())
                return False
        else:
            print_error("No servers were started")
            return False
            
    except ImportError as e:
        print_error(f"Cannot import required modules: {e}")
        return False
    except Exception as e:
        print_error(f"Error testing server startup: {e}")
        return False

def check_environment():
    """Check the environment for common issues"""
    print_header("Environment Check")
    
    issues = []
    
    # Check Python version
    python_version = sys.version_info
    if python_version < (3, 9):
        issues.append(f"Python version too old: {python_version.major}.{python_version.minor} (need >= 3.9)")
    else:
        print_success(f"Python version OK: {python_version.major}.{python_version.minor}.{python_version.micro}")
    
    # Check for required modules
    required_modules = ['pydantic_ai', 'asyncio']
    for module in required_modules:
        try:
            __import__(module)
            print_success(f"Module {module} available")
        except ImportError:
            issues.append(f"Missing module: {module}")
    
    # Check for tomllib/tomli
    try:
        import tomllib
        print_success("Module tomllib available")
    except ImportError:
        try:
            import tomli
            print_success("Module tomli available (fallback)")
        except ImportError:
            issues.append("Neither tomllib nor tomli available")
    
    # Check OpenAI API key
    api_key = os.getenv("OPENAI_API_KEY")
    if api_key:
        print_success("OpenAI API key is set")
    else:
        print_warning("OpenAI API key not set (some features may not work)")
    
    # Check servers.toml
    servers_toml = Path(__file__).parent.parent / "servers.toml"
    if servers_toml.exists():
        print_success(f"servers.toml found")
    else:
        issues.append("servers.toml not found")
    
    # Check server scripts
    base_path = Path(__file__).parent.parent
    server_scripts = [
        "mcp_servers/mcp_data_extractor.py",
        "mcp_servers/mcp_analysis.py"
    ]
    
    for script in server_scripts:
        script_path = base_path / script
        if script_path.exists():
            print_success(f"Server script found: {script}")
        else:
            issues.append(f"Server script missing: {script}")
    
    if issues:
        print_error("Environment issues found:")
        for issue in issues:
            print_error(f"  - {issue}")
        return False
    else:
        print_success("Environment check passed")
        return True

def provide_recommendations():
    """Provide recommendations for fixing issues"""
    print_header("Troubleshooting Recommendations")
    
    print_info("If you're experiencing MCP server issues, try these steps:")
    print("1. Run this script with --kill-all to clean up stuck processes")
    print("2. Check that all required Python packages are installed:")
    print("   pip install pydantic-ai asyncio psutil")
    print("3. Make sure your OpenAI API key is set:")
    print("   export OPENAI_API_KEY=your_key_here")
    print("4. Verify server scripts exist in mcp_servers/ directory")
    print("5. Try restarting your terminal/IDE to clear environment issues")
    print("6. If using Streamlit, try refreshing the page")
    
    print_info("\nFor persistent issues:")
    print("- Check system resources (memory, disk space)")
    print("- Look for firewall/antivirus blocking Python processes")
    print("- Check file permissions on server scripts")
    print("- Review logs for specific error messages")

def main():
    """Main function"""
    parser = argparse.ArgumentParser(description="MCP Server Recovery Utility")
    parser.add_argument(
        "--kill-all",
        action="store_true",
        help="Force kill all MCP server processes"
    )
    parser.add_argument(
        "--restart",
        action="store_true",
        help="Restart MCP servers after cleanup"
    )
    parser.add_argument(
        "--check-only",
        action="store_true",
        help="Only check status, don't make changes"
    )
    parser.add_argument(
        "--results-dir",
        help="Specify results directory for testing"
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable verbose logging"
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Use force when killing processes"
    )
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    print(f"{Colors.BOLD}{Colors.WHITE}MCP Server Recovery Utility{Colors.END}")
    print(f"{Colors.WHITE}{'='*40}{Colors.END}")
    
    # Always check environment first
    env_ok = check_environment()
    
    # Check current status
    current_processes = check_server_status()
    
    if args.check_only:
        if current_processes:
            print_info(f"Found {len(current_processes)} MCP-related processes")
        else:
            print_info("No MCP processes currently running")
        provide_recommendations()
        return
    
    # Kill processes if requested or if we found stuck ones
    if args.kill_all or current_processes:
        if args.kill_all:
            print_warning("--kill-all specified, terminating all MCP processes")
        else:
            print_info("Found existing MCP processes, cleaning up...")
        
        cleanup_success = cleanup_mcp_processes(force=args.force)
        
        # Wait a moment for cleanup to complete
        time.sleep(2)
        
        if not cleanup_success:
            print_error("Cleanup failed, some processes may still be running")
            if not args.force:
                print_info("Try running with --force flag for aggressive cleanup")
    
    # Restart servers if requested
    if args.restart:
        if env_ok:
            print_info("Attempting to restart MCP servers...")
            restart_success = test_server_startup(args.results_dir)
            
            if restart_success:
                print_success("MCP servers restarted successfully!")
            else:
                print_error("Failed to restart MCP servers")
                provide_recommendations()
        else:
            print_error("Environment issues prevent server restart")
            provide_recommendations()
    
    # Final status check
    final_processes = find_mcp_processes() + find_python_processes_with_mcp()
    
    if not final_processes:
        print_success("No MCP processes currently running")
    else:
        print_info(f"Current MCP processes: {len(final_processes)}")
        for proc in final_processes:
            print(f"  PID {proc['pid']}: {proc['status']}")
    
    if not env_ok:
        provide_recommendations()

if __name__ == "__main__":
    main()
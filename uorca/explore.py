"""
UORCA Explorer CLI Wrapper

This module provides a simplified CLI wrapper for the UORCA Explorer Streamlit application
that runs directly without container dependency.
"""

import sys
import os
import socket
import subprocess
from pathlib import Path
from dotenv import load_dotenv, find_dotenv


def is_port_available(host, port):
    """
    Check if a port is available on the given host.

    Args:
        host: Host address to check
        port: Port number to check

    Returns:
        bool: True if port is available, False otherwise
    """
    try:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
            sock.settimeout(1)
            result = sock.connect_ex((host, port))
            return result != 0
    except Exception:
        return False


def find_available_port(host, start_port, max_attempts=10):
    """
    Find an available port starting from start_port.

    Args:
        host: Host address to check
        start_port: Starting port number
        max_attempts: Maximum number of ports to try

    Returns:
        int: Available port number, or None if none found
    """
    for port in range(start_port, start_port + max_attempts):
        if is_port_available(host, port):
            return port
    return None


def attempt_port_cleanup(host, port):
    """
    Attempt to find and terminate processes using the specified port.

    Args:
        host: Host address
        port: Port number to clean up

    Returns:
        bool: True if cleanup was attempted, False otherwise
    """
    try:
        import signal
        import time

        # Try to find processes using the port on Unix-like systems
        if os.name == 'posix':
            try:
                # Use lsof to find processes using the port
                result = subprocess.run(['lsof', '-ti', f':{port}'],
                                      capture_output=True, text=True, timeout=5)
                if result.returncode == 0 and result.stdout.strip():
                    pids = result.stdout.strip().split('\n')
                    print(f"Found {len(pids)} process(es) using port {port}. Attempting cleanup...")

                    for pid in pids:
                        try:
                            pid = int(pid.strip())
                            os.kill(pid, signal.SIGTERM)
                            print(f"Sent SIGTERM to process {pid}")
                        except (ValueError, ProcessLookupError, PermissionError) as e:
                            print(f"Could not terminate process {pid}: {e}")

                    # Give processes time to terminate gracefully
                    time.sleep(2)
                    return True

            except (subprocess.TimeoutExpired, FileNotFoundError):
                # lsof not available or timed out
                pass

        # On Windows, try netstat approach
        elif os.name == 'nt':
            try:
                result = subprocess.run(['netstat', '-ano'],
                                      capture_output=True, text=True, timeout=5)
                if result.returncode == 0:
                    lines = result.stdout.split('\n')
                    for line in lines:
                        if f':{port} ' in line and 'LISTENING' in line:
                            parts = line.split()
                            if len(parts) >= 5:
                                try:
                                    pid = int(parts[-1])
                                    subprocess.run(['taskkill', '/F', '/PID', str(pid)],
                                                 capture_output=True, timeout=5)
                                    print(f"Terminated process {pid} using port {port}")
                                    time.sleep(1)
                                    return True
                                except (ValueError, subprocess.TimeoutExpired):
                                    pass
            except (subprocess.TimeoutExpired, FileNotFoundError):
                pass

    except Exception as e:
        print(f"Port cleanup attempt failed: {e}")

    return False


def main(results_dir=None, port=8501, host="127.0.0.1", headless=False):
    """
    Launch UORCA Explorer Streamlit application directly.

    Args:
        results_dir: Path to UORCA results directory
        port: Port number for the web application
        host: Host address to bind to
        headless: Run in headless mode (no browser auto-open)
    """

    # Load environment variables from .env file
    current_dir = Path(__file__).parent
    project_root = current_dir.parent
    env_file = project_root / ".env"
    if env_file.exists():
        load_dotenv(env_file)
    else:
        # Try to find .env file automatically
        load_dotenv(find_dotenv())

    # Check for OpenAI API key (required for AI features in UORCA Explorer)
    if not os.getenv("OPENAI_API_KEY"):
        print("Warning: OPENAI_API_KEY not found.")
        print("AI-powered analysis features will be disabled in UORCA Explorer.")
        print("To enable AI features, add OPENAI_API_KEY=your_key to your .env file.")
        print("")

    # Validate results directory if provided
    if results_dir:
        results_path = Path(results_dir).resolve()
        if not results_path.exists():
            print(f"Error: Results directory does not exist: {results_dir}")
            sys.exit(1)
        if not results_path.is_dir():
            print(f"Error: Results path is not a directory: {results_dir}")
            sys.exit(1)

        # Check if it looks like a UORCA results directory
        # Look for the specific UORCA structure: GSE*/metadata/analysis_info.json, etc.
        has_analysis_info = list(results_path.glob("*/metadata/analysis_info.json"))
        has_contrasts = list(results_path.glob("*/metadata/contrasts.csv"))
        has_cpm_data = list(results_path.glob("*/RNAseqAnalysis/CPM.csv"))

        if not (has_analysis_info or has_contrasts or has_cpm_data):
            print(f"Error: Directory does not contain UORCA results: {results_dir}")
            print("Expected structure: GSE*/metadata/analysis_info.json, GSE*/metadata/contrasts.csv, GSE*/RNAseqAnalysis/CPM.csv")
            print("Please provide a valid UORCA results directory.")
            sys.exit(1)

        # Count analysis folders using same logic as ResultsIntegrator
        analysis_folders = []

        # Check if the results directory itself contains RNAseqAnalysis
        rnaseq_dir = results_path / "RNAseqAnalysis"
        if rnaseq_dir.is_dir():
            analysis_folders.append(results_path.name)

        # If no analysis folders found in current directory, look for subdirectories
        if not analysis_folders:
            for item in results_path.iterdir():
                if item.is_dir():
                    # Check if it contains an RNAseqAnalysis subdirectory
                    rnaseq_subdir = item / "RNAseqAnalysis"
                    if rnaseq_subdir.is_dir():
                        analysis_folders.append(item.name)

        # Validation complete - directory contains UORCA results
        pass

    # Set up environment for Streamlit
    env = os.environ.copy()

    # Ensure .env variables are available to the subprocess
    # Re-load to make sure all variables are in the current environment
    if env_file.exists():
        load_dotenv(env_file, override=False)  # Don't override existing env vars

    # Copy all current environment variables to subprocess environment
    env.update(os.environ)

    # Configure UORCA Explorer
    if results_dir:
        env['UORCA_DEFAULT_RESULTS_DIR'] = str(Path(results_dir).resolve())

    # Configure Streamlit
    env['STREAMLIT_BROWSER_GATHER_USAGE_STATS'] = 'false'
    env['STREAMLIT_SERVER_HEADLESS'] = 'false'
    env['STREAMLIT_LOGGER_LEVEL'] = 'ERROR'
    env['STREAMLIT_CLIENT_SHOW_ERROR_DETAILS'] = 'false'
    env['STREAMLIT_SERVER_ENABLE_CORS'] = 'false'
    env['STREAMLIT_SERVER_ENABLE_XSRF_PROTECTION'] = 'false'
    # Ensure browser opens to 127.0.0.1 regardless of server binding address
    env['STREAMLIT_BROWSER_SERVER_ADDRESS'] = '127.0.0.1'

    # Get path to uorca_explorer.py
    current_dir = Path(__file__).parent
    project_root = current_dir.parent
    explorer_script = project_root / "main_workflow" / "reporting" / "uorca_explorer.py"

    if not explorer_script.exists():
        print(f"Error: Could not find uorca_explorer.py at {explorer_script}")
        sys.exit(1)

    # Check if the requested port is available
    if not is_port_available(host, port):
        print(f"Port {port} is already in use.")

        # Attempt to clean up the port
        print("Attempting to free up the port...")
        cleanup_success = attempt_port_cleanup(host, port)

        # Check if port is now available after cleanup attempt
        if cleanup_success and is_port_available(host, port):
            print(f"Successfully freed up port {port}")
        else:
            print("Could not free up the port. Finding an available port...")
            available_port = find_available_port(host, port)
            if available_port is None:
                print(f"Error: Could not find an available port starting from {port}")
                print("Please try a different port or manually stop the existing service.")
                sys.exit(1)
            else:
                print(f"Using port {available_port} instead of {port}")
                port = available_port

    # Build streamlit command for direct execution
    cmd = [
        'uv', 'run', 'streamlit', 'run',
        str(explorer_script),
        '--server.port', str(port),
        '--server.address', host,
        '--server.headless', str(not headless).lower(),
        '--browser.serverAddress', '127.0.0.1',
        '--logger.level', 'error'
    ]

    print("=" * 50)
    print("UORCA Explorer")
    print("=" * 50)

    # Determine which URL(s) to show the user
    if host == "0.0.0.0":
        # When binding to all interfaces, show localhost (most common) and mention network access
        primary_url = f"http://127.0.0.1:{port}"
        print(f"Starting UORCA Explorer on all network interfaces")
        if results_dir:
            print(f"Using results directory: {results_dir}")
        print("")
        print(f"Please navigate to: {primary_url}")
        print("(Also accessible from other devices on your network via your machine's IP)")
    else:
        # Show the specific host they requested
        primary_url = f"http://{host}:{port}"
        print(f"Starting UORCA Explorer on {primary_url}")
        if results_dir:
            print(f"Using results directory: {results_dir}")
        print("")
        print(f"Please navigate to {primary_url} in your browser")
    print("")
    print("Press Ctrl+C to stop the application")
    print("=" * 50)

    try:
        # Change to project root for proper imports
        os.chdir(project_root)

        # Run streamlit directly with selective output filtering
        process = subprocess.Popen(cmd, env=env, stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT, text=True, bufsize=1)

        # Filter out unwanted Streamlit messages while preserving important ones
        for line in process.stdout:
            line = line.strip()
            # Skip common Streamlit noise
            if any(skip_phrase in line.lower() for skip_phrase in [
                'welcome to streamlit',
                'if you\'d like to receive helpful onboarding emails',
                'please enter your email address',
                'you can find our privacy policy',
                'this open source library collects usage statistics',
                'telemetry data is stored',
                'if you\'d like to opt out',
                'creating that file if necessary',
                'you can now view your streamlit app',
                'for better performance, install the watchdog'
            ]):
                continue
            # Show important messages
            if line and not line.startswith('  '):
                print(f"  {line}")

        process.wait()
        if process.returncode != 0:
            raise subprocess.CalledProcessError(process.returncode, cmd)

    except subprocess.CalledProcessError as e:
        print(f"\nError: Failed to start UORCA Explorer: {e}")
        sys.exit(1)
    except FileNotFoundError:
        print("\nError: Could not find 'uv' command.")
        print("Make sure you're in the UORCA environment with uv installed.")
        sys.exit(1)
    except KeyboardInterrupt:
        print("\n\nStopping UORCA Explorer...")
        print("Thanks for using UORCA Explorer! Hopefully you found the analyses to be well orca-strated!")
        sys.exit(0)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description='Launch UORCA Explorer Streamlit application'
    )
    parser.add_argument('results_dir', nargs='?',
                       help='Path to UORCA results directory')
    parser.add_argument('--port', type=int, default=8501,
                       help='Port number for the web application')
    parser.add_argument('--host', default='127.0.0.1',
                       help='Host address to bind to')
    parser.add_argument('--headless', action='store_true',
                       help='Run in headless mode (no browser auto-open)')

    args = parser.parse_args()
    main(args.results_dir, args.port, args.host, args.headless)

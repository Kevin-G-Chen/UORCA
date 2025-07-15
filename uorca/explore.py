"""
UORCA Explorer CLI Wrapper

This module provides a simplified CLI wrapper for the UORCA Explorer Streamlit application
that runs directly without container dependency.
"""

import sys
import os
import subprocess
from pathlib import Path
from dotenv import load_dotenv, find_dotenv


def main(results_dir=None, port=8501, host="0.0.0.0", headless=False):
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

        if analysis_folders:
            print(f"Found {len(analysis_folders)} datasets with complete analysis results")
        else:
            print(f"Found partial results: {len(has_analysis_info)} analysis_info, {len(has_contrasts)} contrasts, {len(has_cpm_data)} CPM files")

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
    env['STREAMLIT_SERVER_HEADLESS'] = 'true' if headless else 'false'
    env['STREAMLIT_LOGGER_LEVEL'] = 'ERROR'
    env['STREAMLIT_CLIENT_SHOW_ERROR_DETAILS'] = 'false'
    env['STREAMLIT_SERVER_ENABLE_CORS'] = 'false'
    env['STREAMLIT_SERVER_ENABLE_XSRF_PROTECTION'] = 'false'

    # Get path to uorca_explorer.py
    current_dir = Path(__file__).parent
    project_root = current_dir.parent
    explorer_script = project_root / "main_workflow" / "reporting" / "uorca_explorer.py"

    if not explorer_script.exists():
        print(f"Error: Could not find uorca_explorer.py at {explorer_script}")
        sys.exit(1)

    # Build streamlit command for direct execution
    cmd = [
        'uv', 'run', 'streamlit', 'run',
        str(explorer_script),
        '--server.port', str(port),
        '--server.address', host,
        '--server.headless', 'true' if headless else 'false',
        '--logger.level', 'error'
    ]

    print("=" * 50)
    print("UORCA Explorer - Direct Python Execution")
    print("=" * 50)
    # Always show localhost URL for user convenience
    local_url = f"http://127.0.0.1:{port}"
    print(f"Starting UORCA Explorer on {local_url}")
    if results_dir:
        print(f"Using results directory: {results_dir}")
    print("")
    print("The application should open automatically in your browser.")
    print(f"If it doesn't, manually navigate to: {local_url}")
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
        print("\n\nUORCA Explorer stopped.")
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
    parser.add_argument('--host', default='0.0.0.0',
                       help='Host address to bind to')
    parser.add_argument('--headless', action='store_true',
                       help='Run in headless mode (no browser auto-open)')

    args = parser.parse_args()
    main(args.results_dir, args.port, args.host, args.headless)

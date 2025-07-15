"""
UORCA Explorer CLI Wrapper

This module provides a simplified CLI wrapper for the UORCA Explorer Streamlit application
that runs directly without container dependency.
"""

import sys
import os
import subprocess
from pathlib import Path


def main(results_dir=None, port=8501, host="0.0.0.0", headless=False):
    """
    Launch UORCA Explorer Streamlit application directly.

    Args:
        results_dir: Path to UORCA results directory
        port: Port number for the web application
        host: Host address to bind to
        headless: Run in headless mode (no browser auto-open)
    """

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
        has_results = any([
            list(results_path.glob("*/analysis_info.json")),
            list(results_path.glob("*/contrasts.csv")),
            list(results_path.glob("*/*.csv"))
        ])

        if not has_results:
            print(f"Warning: Directory may not contain UORCA results: {results_dir}")
            print("Proceeding anyway...")

    # Set up environment for Streamlit
    env = os.environ.copy()

    # Configure UORCA Explorer
    if results_dir:
        env['UORCA_DEFAULT_RESULTS_DIR'] = str(Path(results_dir).resolve())

    # Configure Streamlit
    env['STREAMLIT_BROWSER_GATHER_USAGE_STATS'] = 'false'
    env['STREAMLIT_SERVER_HEADLESS'] = 'true' if headless else 'false'

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
        '--server.headless', 'true' if headless else 'false'
    ]

    print("=" * 50)
    print("UORCA Explorer - Direct Python Execution")
    print("=" * 50)
    print(f"Starting UORCA Explorer on http://{host}:{port}")
    if results_dir:
        print(f"Using results directory: {results_dir}")
    print("Press Ctrl+C to stop the application")
    print("=" * 50)

    try:
        # Change to project root for proper imports
        os.chdir(project_root)

        # Run streamlit directly
        subprocess.run(cmd, env=env, check=True)

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

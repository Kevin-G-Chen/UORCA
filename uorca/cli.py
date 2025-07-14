"""
UORCA Unified Command Line Interface

This module provides a unified CLI for all UORCA workflows:
- identify: Dataset identification from GEO
- run: RNA-seq analysis pipeline
- explore: Interactive results explorer
"""

import argparse
import sys
import os
from typing import Optional

def main():
    """Main CLI entry point with subcommands."""
    parser = argparse.ArgumentParser(
        prog="uorca",
        description="UORCA - Unified Omics Reference Corpus of Analyses",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  uorca identify -q "cancer stem cell differentiation" -o results.csv
  uorca run --accession GSE123456 --output_dir ../UORCA_results
  uorca explore --results-dir ../UORCA_results --port 8501

For more help on a specific command, use:
  uorca COMMAND --help
        """
    )

    # Add version argument
    parser.add_argument('--version', action='version', version='UORCA 0.1.0')

    # Create subparsers for the three main commands
    subparsers = parser.add_subparsers(
        dest='command',
        help='Available commands',
        metavar='COMMAND'
    )

    # ========================================================================
    # IDENTIFY subcommand - Dataset identification
    # ========================================================================
    identify_parser = subparsers.add_parser(
        'identify',
        help='Identify relevant RNA-seq datasets from GEO',
        description='Identify and evaluate relevant RNA-seq datasets from GEO for a biological research query',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Core parameters for identify
    identify_core = identify_parser.add_argument_group('Core Options', 'Essential parameters for dataset identification')
    identify_core.add_argument('-q', '--query', required=True,
                              help='Biological research query (e.g., "neuroblastoma tumor vs normal")')
    identify_core.add_argument('-o', '--output', default='./dataset_identification_results',
                              help='Output directory for results')
    identify_core.add_argument('-t', '--threshold', type=float, default=7.0,
                              help='Relevance score threshold (0-10) for including datasets')

    # Search parameters for identify
    identify_search = identify_parser.add_argument_group('Search Options', 'Control dataset search and evaluation')
    identify_search.add_argument('--max-datasets', type=int, default=2000,
                                help='Maximum datasets to retrieve per search term')

    # Advanced parameters for identify
    identify_advanced = identify_parser.add_argument_group('Advanced Options', 'Fine-tune algorithm behavior (expert users)')
    identify_advanced.add_argument('--scoring-rounds', type=int, default=3,
                                  help='Number of independent relevance scoring rounds for reliability')
    identify_advanced.add_argument('--batch-size', type=int, default=20,
                                  help='Datasets per AI evaluation batch (affects memory usage)')
    identify_advanced.add_argument('--verbose', action='store_true',
                                   help='Enable verbose logging (DEBUG level)')

    identify_parser.set_defaults(func=run_identify)

    # ========================================================================
    # RUN subcommand - RNA-seq analysis pipeline
    # ========================================================================
    run_parser = subparsers.add_parser(
        'run',
        help='Run RNA-seq analysis pipeline on a GEO dataset',
        description='UORCA RNA-seq Analysis Pipeline',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    run_parser.add_argument("--accession", required=True,
                           help="GEO accession (e.g., GSE123456)")
    run_parser.add_argument("--output_dir", default="../UORCA_results",
                           help="Output directory")
    run_parser.add_argument("--resource_dir", default="./data/kallisto_indices/",
                           help="Resource directory")
    run_parser.add_argument("--cleanup", action="store_true",
                           help="Clean up FASTQ and SRA files after analysis")

    run_parser.set_defaults(func=run_analysis)

    # ========================================================================
    # EXPLORE subcommand - Interactive results explorer
    # ========================================================================
    explore_parser = subparsers.add_parser(
        'explore',
        help='Launch interactive results explorer web app',
        description='Launch the UORCA Explorer Streamlit web application',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    explore_parser.add_argument('--results-dir',
                               help='Results directory to explore (if not specified, will use default detection)')
    explore_parser.add_argument('--port', type=int, default=8501,
                               help='Port number for the web application')
    explore_parser.add_argument('--host', default='0.0.0.0',
                               help='Host address to bind to')
    explore_parser.add_argument('--headless', action='store_true',
                               help='Run in headless mode (no browser auto-open)')

    explore_parser.set_defaults(func=run_explore)

    # Parse arguments
    args = parser.parse_args()

    # If no command specified, show help
    if not args.command:
        parser.print_help()
        sys.exit(1)

    # Call the appropriate function
    try:
        args.func(args)
    except KeyboardInterrupt:
        print("\nOperation cancelled by user.")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


def run_identify(args):
    """Run the dataset identification workflow."""
    # Import here to avoid startup overhead when not needed
    from uorca.identify import main as identify_main

    # Check for required environment variable
    if not os.getenv("ENTREZ_EMAIL"):
        print("ERROR: Email is required for NCBI Entrez API access.")
        print("Please set the ENTREZ_EMAIL environment variable with your email address.")
        print("This is required by NCBI guidelines for API usage.")
        sys.exit(1)

    # Rebuild sys.argv to match what the original script expects
    sys.argv = ['identify']
    sys.argv.extend(['-q', args.query])
    sys.argv.extend(['-o', args.output])
    sys.argv.extend(['-t', str(args.threshold)])
    sys.argv.extend(['--max-datasets', str(args.max_datasets)])
    sys.argv.extend(['--scoring-rounds', str(args.scoring_rounds)])
    sys.argv.extend(['--batch-size', str(args.batch_size)])
    if args.verbose:
        sys.argv.append('--verbose')

    # Call the original main function
    identify_main()


def run_analysis(args):
    """Run the RNA-seq analysis workflow."""
    # Import here to avoid startup overhead when not needed
    from uorca.run import main as run_main

    # Rebuild sys.argv to match what the original script expects
    sys.argv = ['run']
    sys.argv.extend(['--accession', args.accession])
    sys.argv.extend(['--output_dir', args.output_dir])
    sys.argv.extend(['--resource_dir', args.resource_dir])
    if args.cleanup:
        sys.argv.append('--cleanup')

    # Call the original main function
    run_main()


def run_explore(args):
    """Run the interactive results explorer."""
    # Import here to avoid startup overhead when not needed
    import subprocess

    # Set up environment variables for the Streamlit app
    env = os.environ.copy()

    if args.results_dir:
        env['UORCA_DEFAULT_RESULTS_DIR'] = args.results_dir

    # Set Streamlit configuration
    env['STREAMLIT_BROWSER_GATHER_USAGE_STATS'] = 'false'
    env['STREAMLIT_SERVER_HEADLESS'] = 'true' if args.headless else 'false'

    # Build the streamlit command
    cmd = [
        'uv', 'run', 'streamlit', 'run',
        'main_workflow/reporting/uorca_explorer.py',
        '--server.port', str(args.port),
        '--server.address', args.host,
        '--server.headless', 'true' if args.headless else 'false'
    ]

    print(f"Starting UORCA Explorer on http://{args.host}:{args.port}")
    if args.results_dir:
        print(f"Using results directory: {args.results_dir}")
    print("Press Ctrl+C to stop the application")

    try:
        # Run streamlit
        subprocess.run(cmd, env=env, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Failed to start UORCA Explorer: {e}")
        sys.exit(1)
    except FileNotFoundError:
        print("Error: Could not find 'uv' command. Make sure you're in the UORCA environment.")
        sys.exit(1)


if __name__ == "__main__":
    main()

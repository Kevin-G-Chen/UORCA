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
from pathlib import Path
from typing import Optional
from dotenv import load_dotenv, find_dotenv

def main():
    """Main CLI entry point with subcommands."""
    parser = argparse.ArgumentParser(
        prog="uorca",
        description="UORCA - Unified Omics Reference Corpus of Analyses",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  uorca identify -q "cancer stem cell differentiation" -o results.csv
  uorca run slurm --input datasets.csv --output_dir ../UORCA_results
  uorca run local --input identification_results/ --output_dir ../UORCA_results
  uorca explore ../UORCA_results --port 8501

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
    identify_search.add_argument('-m', '--max-per-term', type=int, default=500,
                                help='Maximum datasets to retrieve per search term')
    identify_search.add_argument('-a', '--max-assess', type=int, default=300,
                                help='Maximum number of datasets to assess for relevance (distributed across clusters)')

    # Advanced parameters for identify
    identify_advanced = identify_parser.add_argument_group('Advanced Options', 'Fine-tune algorithm behavior (expert users)')
    identify_advanced.add_argument('--cluster-divisor', type=int, default=10,
                                  help='Divisor for cluster count (total_datasets / divisor). Smaller values = more clusters.')
    identify_advanced.add_argument('-r', '--rounds', type=int, default=3,
                                  help='Number of independent relevance scoring rounds for reliability')
    identify_advanced.add_argument('-b', '--batch-size', type=int, default=20,
                                  help='Datasets per AI evaluation batch (affects memory usage)')
    identify_advanced.add_argument('--model', type=str, default='gpt-5-mini',
                                  help='OpenAI model to use for relevance assessment')
    identify_advanced.add_argument('-v', '--verbose', action='store_true',
                                   help='Enable verbose logging (DEBUG level)')

    identify_parser.set_defaults(func=run_identify)

    # ========================================================================
    # RUN subcommand - Batch RNA-seq analysis pipeline
    # ========================================================================
    run_parser = subparsers.add_parser(
        'run',
        help='Run batch RNA-seq analysis pipeline',
        description='UORCA Batch RNA-seq Analysis Pipeline',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Create subparsers for batch systems
    run_subparsers = run_parser.add_subparsers(
        dest='batch_system',
        help='Batch processing system',
        metavar='SYSTEM'
    )

    # SLURM batch processing
    slurm_parser = run_subparsers.add_parser(
        'slurm',
        help='Run batch processing on SLURM cluster',
        description='Submit batch jobs to SLURM for dataset processing. SLURM-specific parameters (partition, constraint, CPUs, memory, time limits) should be configured in slurm_config.yaml.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
Configuration:
  SLURM job parameters are configured via YAML file (default: slurm_config.yaml):
    - partition, constraint, cpus_per_task, memory, time_limit
    - container settings (engine, image paths)
    - resource management (max_parallel, max_storage_gb)

  CLI arguments override workflow-level settings but not SLURM job parameters.

Example slurm_config.yaml:
  slurm:
    partition: "tki_agpdev"
    constraint: "clx"
    cpus_per_task: 12
    memory: "16G"
    time_limit: "6:00:00"
        """
    )

    slurm_parser.add_argument("--input", required=True,
                             help="Either a CSV file with dataset information (must have 'Accession' column) or a directory containing both CSV and identification metadata")
    slurm_parser.add_argument("--config",
                             help="YAML configuration file for SLURM settings (default: slurm_config.yaml if exists)")
    slurm_parser.add_argument("--output_dir", default="../UORCA_results",
                             help="Output directory for results")
    slurm_parser.add_argument("--resource_dir", default="./data/kallisto_indices/",
                             help="Resource directory for Kallisto indices")
    slurm_parser.add_argument("--max_parallel", type=int, default=10,
                             help="Maximum number of parallel jobs")
    slurm_parser.add_argument("--max_storage_gb", type=float, default=500,
                             help="Maximum storage usage in GB")
    slurm_parser.add_argument("--no-cleanup", action="store_true",
                             help="Skip cleanup of FASTQ and SRA files after analysis")

    slurm_parser.set_defaults(func=run_batch_slurm)

    # Local batch processing
    local_parser = run_subparsers.add_parser(
        'local',
        help='Run batch processing locally with parallel jobs',
        description='Process datasets locally using multiprocessing',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    local_parser.add_argument("--input", required=True,
                             help="Either a CSV file with dataset information (must have 'Accession' column) or a directory containing both CSV and identification metadata")
    local_parser.add_argument("--output_dir", default="../UORCA_results",
                             help="Output directory for results")
    local_parser.add_argument("--resource_dir", default="./data/kallisto_indices/",
                             help="Resource directory for Kallisto indices")
    local_parser.add_argument("--max_workers", type=int,
                             help="Maximum number of parallel workers (default: auto-detect 75%% of CPU cores)")
    local_parser.add_argument("--max_storage_gb", type=float,
                             help="Maximum storage usage in GB (default: min(50 GB, 75%% of available disk free space))")
    local_parser.add_argument("--no-cleanup", action="store_true",
                             help="Skip cleanup of FASTQ and SRA files after analysis")
    local_parser.add_argument("--timeout_hours", type=float, default=6,
                             help="Timeout per job in hours")
    local_parser.add_argument("--container_tmpfs_gb", type=int,
                             help="Size of container /tmp tmpfs in GB (default: 20). Set 0 to disable tmpfs mount.")


    local_parser.set_defaults(func=run_batch_local)

    # ========================================================================
    # EXPLORE subcommand - Interactive results explorer
    # ========================================================================
    explore_parser = subparsers.add_parser(
        'explore',
        help='Launch interactive results explorer web app',
        description='Launch the UORCA Explorer Streamlit web application',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    explore_parser.add_argument('results_dir', nargs='?',
                               help='Results directory to explore (if not specified, will use default detection)')
    explore_parser.add_argument('--port', type=int, default=8501,
                               help='Port number for the web application')
    explore_parser.add_argument('--host', default='127.0.0.1',
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

    # Special handling for run command without batch system
    if args.command == 'run' and not hasattr(args, 'batch_system') or (hasattr(args, 'batch_system') and not args.batch_system):
        run_parser.print_help()
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


def _load_environment_variables():
    """Load environment variables from .env file and validate requirements."""
    # Load environment variables from the project root .env file
    project_root = Path(__file__).resolve().parent.parent
    env_file = project_root / ".env"
    if env_file.exists():
        load_dotenv(env_file)
        print(f"Loaded environment variables from: {env_file}")
    else:
        load_dotenv(find_dotenv())
        print("Loaded environment variables from system .env file")

def _check_environment_requirements(require_openai: bool = False):
    """Check environment requirements and provide clear instructions if missing."""
    missing = []

    # Check required variables
    if not os.getenv("ENTREZ_EMAIL"):
        missing.append("ENTREZ_EMAIL")

    if require_openai and not os.getenv("OPENAI_API_KEY"):
        missing.append("OPENAI_API_KEY")

    # Inform about API key status
    entrez_api_key_set = bool(os.getenv("ENTREZ_API_KEY"))
    openai_key_set = bool(os.getenv("OPENAI_API_KEY"))

    if not entrez_api_key_set:
        print("ℹ️  ENTREZ_API_KEY not set - using slower rate limits (3 requests/second).")
        print("   For faster processing, consider getting a free NCBI API key.")

    if not require_openai and not openai_key_set:
        print("ℹ️  OPENAI_API_KEY not set - AI-powered analysis features will be disabled.")
    
    # If missing required variables, show comprehensive instructions and exit
    if missing:
        print("\n" + "="*70)
        print("❌ SETUP REQUIRED: Missing essential API credentials")
        print("="*70)

        print("\n📋 What's missing:")
        for var in missing:
            if var == "ENTREZ_EMAIL":
                print("  • ENTREZ_EMAIL - Required by NCBI for accessing biological databases")
            elif var == "OPENAI_API_KEY":
                print("  • OPENAI_API_KEY - Required for AI-powered dataset identification")

        env_file_path = Path(__file__).resolve().parent.parent / '.env'
        print(f"\n🔧 How to fix this:")
        print(f"   Create or edit the file: {env_file_path}")

        print(f"\n📝 Add these lines to your .env file:")
        for var in missing:
            if var == "ENTREZ_EMAIL":
                print(f"   {var}=your.email@institution.edu")
            elif var == "OPENAI_API_KEY":
                print(f"   {var}=sk-proj-your-key-here")

        print(f"\n🔗 Where to get API keys:")
        if "OPENAI_API_KEY" in missing:
            print(f"   • OpenAI API key: https://platform.openai.com/api-keys")
        
        print("="*70)
        sys.exit(1)


def _check_kallisto_indices(resource_dir: str):
    """Check if Kallisto indices are available and provide clear instructions if missing."""
    from pathlib import Path
    
    # Convert to Path object
    indices_path = Path(resource_dir).resolve()
    
    # Define expected species
    expected_species = ['human', 'mouse', 'dog', 'monkey', 'zebrafish']
    
    # Check if the directory exists
    if not indices_path.exists():
        print("\n" + "="*70)
        print("⚠️  ERROR: Kallisto indices directory not found!")
        print("="*70)
        print(f"\nExpected directory: {indices_path}")
        print("\nKallisto indices are REQUIRED for RNA-seq quantification.")
        print("\nTo download the indices, run:")
        print("\n  ./download_kallisto_indices.sh")
        print("\nThis will download indices for: human, mouse, dog, monkey, zebrafish")
        print("\nAlternatively, download specific species:")
        print("  ./download_kallisto_indices.sh human")
        print("  ./download_kallisto_indices.sh mouse")
        print("\nIf you have indices in a different location, specify with --resource_dir")
        print("="*70)
        sys.exit(1)
    
    # Check which species indices are present
    found_species = []
    missing_species = []
    
    for species in expected_species:
        species_dir = indices_path / species
        index_file = species_dir / "index.idx"
        t2g_file = species_dir / "t2g.txt"
        
        if species_dir.exists() and index_file.exists() and t2g_file.exists():
            found_species.append(species)
        else:
            missing_species.append(species)
    
    # Report status
    if not found_species:
        print("\n" + "="*70)
        print("⚠️  ERROR: No Kallisto indices found!")
        print("="*70)
        print(f"\nChecked directory: {indices_path}")
        print("\nNo valid indices found for any species.")
        print("\nTo download ALL indices, run:")
        print("\n  ./download_kallisto_indices.sh")
        print("\nTo download specific species:")
        print("  ./download_kallisto_indices.sh human")
        print("  ./download_kallisto_indices.sh mouse")
        print("\nExpected structure:")
        print(f"  {indices_path}/")
        print("    ├── human/")
        print("    │   ├── index.idx")
        print("    │   └── t2g.txt")
        print("    ├── mouse/")
        print("    │   ├── index.idx")
        print("    │   └── t2g.txt")
        print("    └── ...")
        print("="*70)
        sys.exit(1)
    
    # Warn if some species are missing
    if missing_species:
        print(f"\n✓ Found Kallisto indices for: {', '.join(found_species)}")
        print(f"ℹ️  Missing indices for: {', '.join(missing_species)}")
        print(f"   To download missing species, run: ./download_kallisto_indices.sh {missing_species[0]}")
    else:
        print(f"✓ All Kallisto indices found in: {indices_path}")


def run_identify(args):
    """Run the dataset identification workflow."""
    # Import here to avoid startup overhead when not needed
    from uorca.identify import main as identify_main

    # Load and validate environment variables
    _load_environment_variables()
    _check_environment_requirements(require_openai=True)  # Identify needs OpenAI

    # Rebuild sys.argv to match what the original script expects
    sys.argv = ['identify']
    sys.argv.extend(['-q', args.query])
    sys.argv.extend(['-o', args.output])
    sys.argv.extend(['-t', str(args.threshold)])
    sys.argv.extend(['-m', str(args.max_per_term)])
    sys.argv.extend(['--cluster-divisor', str(args.cluster_divisor)])
    sys.argv.extend(['-a', str(args.max_assess)])
    sys.argv.extend(['-r', str(args.rounds)])
    sys.argv.extend(['-b', str(args.batch_size)])
    sys.argv.extend(['--model', args.model])
    if args.verbose:
        sys.argv.append('-v')

    # Call the original main function
    identify_main()


def run_batch_slurm(args):
    """Run batch processing using SLURM."""
    from uorca.batch.slurm import SlurmBatchProcessor
    from pathlib import Path

    print("Starting SLURM batch processing...")

    # Load and validate environment variables
    _load_environment_variables()
    _check_environment_requirements(require_openai=True)  # Pipeline needs OpenAI
    
    # Check for Kallisto indices
    _check_kallisto_indices(args.resource_dir)

    try:
        # Handle config file auto-detection
        config_file = args.config
        if not config_file:
            # Auto-detect slurm_config.yaml in current directory
            default_config = Path("slurm_config.yaml")
            if default_config.exists():
                config_file = str(default_config)
                print(f"Auto-detected config file: {config_file}")

        # Initialize processor with config file
        processor = SlurmBatchProcessor(config_file=config_file)

        # Prepare parameters - CLI args override config file settings
        params = {
            'max_parallel': args.max_parallel,
            'max_storage_gb': args.max_storage_gb,
            'cleanup': not args.no_cleanup,
            'resource_dir': args.resource_dir
        }

        # Submit jobs
        jobs_submitted = processor.submit_datasets(args.input, args.output_dir, **params)

        if jobs_submitted > 0:
            print(f"\nSuccessfully submitted {jobs_submitted} jobs to SLURM")
        else:
            print("\nNo jobs were submitted")

    except Exception as e:
        print(f"\nError: {e}")
        sys.exit(1)


def run_batch_local(args):
    """Run batch processing locally."""
    from uorca.batch import get_batch_processor

    print("Starting local batch processing...")

    # Load and validate environment variables
    _load_environment_variables()
    _check_environment_requirements(require_openai=True)  # Pipeline needs OpenAI
    
    # Check for Kallisto indices
    _check_kallisto_indices(args.resource_dir)

    try:
        processor = get_batch_processor('local')

        # Prepare parameters
        params = {
            'cleanup': not args.no_cleanup,
            'resource_dir': args.resource_dir,
            'timeout_hours': args.timeout_hours
        }

        # Add optional parameters if specified
        if args.max_workers:
            params['max_workers'] = args.max_workers
        if args.max_storage_gb:
            params['max_storage_gb'] = args.max_storage_gb
        if args.container_tmpfs_gb is not None:
            params['container_tmpfs_gb'] = args.container_tmpfs_gb

        # Submit jobs
        jobs_submitted = processor.submit_datasets(args.input, args.output_dir, **params)

        if jobs_submitted > 0:
            print(f"\nSuccessfully processed {jobs_submitted} datasets")
        else:
            print("\nNo datasets were processed")

    except Exception as e:
        print(f"\nError: {e}")
        sys.exit(1)


def run_explore(args):
    """Run the interactive results explorer."""
    # Import here to avoid startup overhead when not needed
    from uorca.explore import main as explore_main

    # Call the simplified explore main function directly
    explore_main(
        results_dir=args.results_dir,
        port=args.port,
        host=args.host,
        headless=args.headless
    )


if __name__ == "__main__":
    main()

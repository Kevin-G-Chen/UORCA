"""
UORCA Dataset Identification CLI Wrapper

This module provides a CLI wrapper for the dataset identification functionality.
"""

import sys
import os
from pathlib import Path
from dotenv import load_dotenv, find_dotenv

def main():
    """
    Launch UORCA dataset identification functionality.
    """
    # Get path to the main workflow directory
    current_dir = Path(__file__).parent
    project_root = current_dir.parent

    # Load environment variables from .env file in project root
    env_file = project_root / ".env"
    if env_file.exists():
        load_dotenv(env_file)
    else:
        # Try to find .env file automatically
        load_dotenv(find_dotenv())

    # Check for required environment variables
    if not os.getenv("ENTREZ_EMAIL"):
        print("ERROR: Email is required for NCBI Entrez API access.")
        print("Please set the ENTREZ_EMAIL environment variable with your email address.")
        print("This is required by NCBI guidelines for API usage.")
        print("")
        print("You can either:")
        print("1. Add ENTREZ_EMAIL=your.email@example.com to your .env file in the project root")
        print("2. Set the environment variable: export ENTREZ_EMAIL=your.email@example.com")
        sys.exit(1)

    # Inform about OpenAI API key if missing (optional but recommended)
    if not os.getenv("OPENAI_API_KEY"):
        print("Warning: OPENAI_API_KEY not found.")
        print("AI-powered dataset evaluation will be disabled.")
        print("To enable AI features, add OPENAI_API_KEY=your_key to your .env file.")
        print("")
    main_workflow_dir = project_root / "main_workflow"
    dataset_id_script = main_workflow_dir / "dataset_identification" / "DatasetIdentification.py"

    if not dataset_id_script.exists():
        print(f"Error: Could not find DatasetIdentification.py at {dataset_id_script}")
        sys.exit(1)

    # Add main_workflow to path for imports
    sys.path.insert(0, str(main_workflow_dir))

    # Change to project root for proper execution context
    original_cwd = os.getcwd()
    os.chdir(project_root)

    try:
        # Import and run the main function from DatasetIdentification
        from dataset_identification.DatasetIdentification import main as identify_main
        identify_main()

    finally:
        # Restore original directory and clean up path
        os.chdir(original_cwd)
        if str(main_workflow_dir) in sys.path:
            sys.path.remove(str(main_workflow_dir))

# Re-export the main function
__all__ = ["main"]

if __name__ == "__main__":
    main()

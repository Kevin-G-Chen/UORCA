"""
UORCA Dataset Identification CLI Wrapper

This module provides a CLI wrapper for the dataset identification functionality.
"""

from main_workflow.dataset_identification.DatasetIdentification import main

# Re-export the main function
__all__ = ["main"]

if __name__ == "__main__":
    main()

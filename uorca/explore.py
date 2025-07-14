"""
UORCA Explorer CLI Wrapper

This module provides a CLI wrapper for the UORCA Explorer Streamlit application.
"""

import sys
import os

# Add main_workflow/reporting to path and change directory for relative imports
reporting_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'main_workflow', 'reporting')
sys.path.insert(0, reporting_dir)

# Change to reporting directory for relative imports to work
original_cwd = os.getcwd()
os.chdir(reporting_dir)

try:
    from uorca_explorer import main
finally:
    # Restore original directory
    os.chdir(original_cwd)
    # Remove from path to avoid pollution
    if reporting_dir in sys.path:
        sys.path.remove(reporting_dir)

# Re-export the main function
__all__ = ["main"]

if __name__ == "__main__":
    main()

"""
UORCA Dataset Analysis CLI Wrapper

This module provides a CLI wrapper for the main UORCA analysis functionality.
"""

import sys
import os

# Add main_workflow to path and change directory for relative imports
main_workflow_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'main_workflow')
sys.path.insert(0, main_workflow_dir)

# Change to main_workflow directory for relative imports to work
original_cwd = os.getcwd()
os.chdir(main_workflow_dir)

try:
    from master import main
finally:
    # Restore original directory
    os.chdir(original_cwd)
    # Remove from path to avoid pollution
    if main_workflow_dir in sys.path:
        sys.path.remove(main_workflow_dir)

# Re-export the main function
__all__ = ["main"]

if __name__ == "__main__":
    main()

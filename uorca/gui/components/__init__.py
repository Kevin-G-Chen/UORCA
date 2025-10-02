"""
UORCA Explorer Streamlit Tabs Package.

This package contains modular Streamlit tab implementations for the UORCA Explorer
application. Each tab is implemented as a separate module for better maintainability
and development workflow.
"""

# Import main render functions from each tab module
# Note: ai_assistant_tab is imported directly by uorca_explorer.py to avoid circular imports
# with uorca.gui.ai modules

from .heatmap_tab import render_heatmap_tab
from .expression_plots_tab import render_expression_plots_tab
from .analysis_plots_tab import render_analysis_plots_tab
from .datasets_info_tab import render_datasets_info_tab
from .contrasts_info_tab import render_contrasts_info_tab
from .sidebar_controls import render_sidebar_controls

# Import helper functions that might be useful externally
from .helpers import (
    get_integrator,
    cached_identify_important_genes,
    add_custom_css,
    setup_fragment_decorator
)

__all__ = [
    # Tab render functions
    'render_heatmap_tab',
    'render_expression_plots_tab',
    'render_analysis_plots_tab',
    'render_datasets_info_tab',
    'render_contrasts_info_tab',
    'render_sidebar_controls',

    # Helper functions
    'get_integrator',
    'cached_identify_important_genes',
    'add_custom_css',
    'setup_fragment_decorator'
]

__version__ = "1.0.0"
__author__ = "UORCA Development Team"
__description__ = "Modular Streamlit tabs for UORCA RNA-seq analysis exploration"

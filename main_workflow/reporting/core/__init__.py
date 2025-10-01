"""
Backwards compatibility layer for main_workflow.reporting.core.

DEPRECATED: This module has been moved to uorca.core.
Please update your imports to use `from uorca.core import ...` instead.
"""

import warnings

# Issue deprecation warning
warnings.warn(
    "main_workflow.reporting.core is deprecated. Use uorca.core instead.",
    DeprecationWarning,
    stacklevel=2
)

# Re-export from new location
from uorca.core import TaskManager, TaskStatus

__all__ = [
    "TaskManager",
    "TaskStatus",
]

"""
UORCA Batch Processing Module

This module provides batch processing capabilities for UORCA with support for
different HPC systems and local parallel processing.

Supported systems:
- SLURM (via sbatch)
- Local parallel processing
"""

import shutil
import subprocess
from typing import Optional, Dict, Any
from pathlib import Path


def detect_batch_system() -> str:
    """
    Detect available batch systems on the current machine.

    Returns:
        str: The name of the detected system ('slurm', 'local')
    """
    # Check for SLURM
    if shutil.which('sbatch'):
        try:
            result = subprocess.run(['sinfo'], capture_output=True, timeout=5)
            if result.returncode == 0:
                return 'slurm'
        except (subprocess.TimeoutExpired, FileNotFoundError):
            pass

    # Default to local processing
    return 'local'


def get_batch_processor(system: Optional[str] = None):
    """
    Get the appropriate batch processor for the given system.

    Args:
        system: Batch system name ('slurm', 'local', or None for auto-detect)

    Returns:
        Batch processor class instance
    """
    if system is None:
        system = detect_batch_system()

    if system == 'slurm':
        from .slurm import SlurmBatchProcessor
        return SlurmBatchProcessor()
    elif system == 'local':
        from .local import LocalBatchProcessor
        return LocalBatchProcessor()
    else:
        raise ValueError(f"Unsupported batch system: {system}")


def list_available_systems() -> Dict[str, bool]:
    """
    List all supported batch systems and their availability.

    Returns:
        Dict mapping system names to availability status
    """
    systems = {}

    # Check SLURM
    systems['slurm'] = shutil.which('sbatch') is not None

    # Local is always available
    systems['local'] = True

    return systems


def submit_batch_job(csv_file: str, output_dir: str, system: Optional[str] = None, **kwargs) -> int:
    """
    Submit a batch job using the specified or auto-detected system.

    Args:
        csv_file: Path to CSV file with dataset information
        output_dir: Output directory for results
        system: Batch system to use (None for auto-detect)
        **kwargs: Additional system-specific parameters

    Returns:
        Number of jobs submitted
    """
    processor = get_batch_processor(system)
    return processor.submit_datasets(csv_file, output_dir, **kwargs)


# Export main functions
__all__ = [
    'detect_batch_system',
    'get_batch_processor',
    'list_available_systems',
    'submit_batch_job'
]

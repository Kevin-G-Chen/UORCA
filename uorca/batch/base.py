"""
Base Batch Processor Interface

This module defines the abstract base class for all batch processing systems.
All concrete batch processors (SLURM, local, etc.) must inherit from this class.
"""

import json
import os
import pandas as pd
from abc import ABC, abstractmethod
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple


class BatchProcessor(ABC):
    """
    Abstract base class for batch processing systems.

    This class defines the common interface that all batch processors must implement.
    """

    def __init__(self):
        """Initialize the batch processor."""
        self.name = self.__class__.__name__.replace('BatchProcessor', '').lower()

    @abstractmethod
    def submit_datasets(self, csv_file: str, output_dir: str, **kwargs) -> int:
        """
        Submit batch jobs for processing datasets.

        Args:
            csv_file: Path to CSV file containing dataset information
            output_dir: Output directory for results
            **kwargs: Additional system-specific parameters

        Returns:
            Number of jobs submitted

        Raises:
            FileNotFoundError: If CSV file doesn't exist
            ValueError: If CSV format is invalid
        """
        pass

    @abstractmethod
    def check_job_status(self, job_id: str) -> Dict[str, Any]:
        """
        Check the status of a specific job.

        Args:
            job_id: Job identifier

        Returns:
            Dictionary containing job status information
        """
        pass

    @abstractmethod
    def cancel_job(self, job_id: str) -> bool:
        """
        Cancel a running job.

        Args:
            job_id: Job identifier

        Returns:
            True if job was successfully cancelled, False otherwise
        """
        pass

    @abstractmethod
    def get_available_resources(self) -> Dict[str, Any]:
        """
        Get information about available computational resources.

        Returns:
            Dictionary containing resource information (CPUs, memory, etc.)
        """
        pass

    def validate_csv_file(self, csv_file: str) -> pd.DataFrame:
        """
        Validate and load the CSV file containing dataset information.

        Args:
            csv_file: Path to CSV file

        Returns:
            Validated DataFrame

        Raises:
            FileNotFoundError: If CSV file doesn't exist
            ValueError: If required columns are missing
        """
        csv_path = Path(csv_file)
        if not csv_path.exists():
            raise FileNotFoundError(f"CSV file not found: {csv_file}")

        try:
            df = pd.read_csv(csv_path)
        except Exception as e:
            raise ValueError(f"Failed to read CSV file: {e}")

        # Check required columns
        required_columns = ['Accession']
        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            raise ValueError(f"Missing required columns: {missing_columns}")

        # Handle dataset size columns (prefer DatasetSizeGB, fallback to DatasetSizeBytes)
        if 'DatasetSizeGB' in df.columns:
            # Size already in GB, apply 2x safety factor
            df['SafeStorageGB'] = df['DatasetSizeGB'] * 2.0
            print("Using recorded dataset size × 2 as approximation for required storage space")
        elif 'DatasetSizeBytes' in df.columns:
            # Convert from bytes to GB and apply 2x safety factor
            df['SafeStorageGB'] = (df['DatasetSizeBytes'] / (1024**3)) * 2.0
            print("Using recorded dataset size × 2 as approximation for required storage space")
        else:
            # No size information available
            df['SafeStorageGB'] = 0
            print("Warning: Neither DatasetSizeGB nor DatasetSizeBytes column found, using default value of 0")

        return df

    def setup_job_tracking(self, output_dir: str) -> Path:
        """
        Set up job status tracking directory.

        Args:
            output_dir: Base output directory

        Returns:
            Path to job status directory
        """
        status_dir = Path(output_dir) / "job_status"
        status_dir.mkdir(parents=True, exist_ok=True)
        return status_dir

    def create_job_status_file(self, status_dir: Path, accession: str,
                              job_id: str, **extra_info) -> Path:
        """
        Create a job status tracking file.

        Args:
            status_dir: Directory for status files
            accession: Dataset accession
            job_id: Job identifier
            **extra_info: Additional information to store

        Returns:
            Path to created status file
        """
        status_file = status_dir / f"{accession}_status.json"

        status_info = {
            'accession': accession,
            'job_id': job_id,
            'state': 'submitted',
            'submitted_time': datetime.now().isoformat(),
            'system': self.name,
            **extra_info
        }

        with open(status_file, 'w') as f:
            json.dump(status_info, f, indent=2)

        return status_file

    def update_job_status(self, status_file: Path, **updates):
        """
        Update job status information.

        Args:
            status_file: Path to status file
            **updates: Fields to update
        """
        if status_file.exists():
            with open(status_file, 'r') as f:
                status = json.load(f)
        else:
            status = {}

        status.update(updates)
        status['last_updated'] = datetime.now().isoformat()

        with open(status_file, 'w') as f:
            json.dump(status, f, indent=2)

    def get_running_jobs_storage(self, output_dir: str) -> float:
        """
        Calculate total storage being used by currently running jobs.

        Args:
            output_dir: Output directory containing job status

        Returns:
            Total storage in GB
        """
        total_storage = 0.0
        status_dir = Path(output_dir) / "job_status"

        if not status_dir.exists():
            return total_storage

        for status_file in status_dir.glob("*_status.json"):
            try:
                with open(status_file, 'r') as f:
                    status = json.load(f)

                if status.get('state') in ['submitted', 'running']:
                    total_storage += status.get('storage_gb', 0)
            except (json.JSONDecodeError, FileNotFoundError):
                continue

        return total_storage

    def display_scheduling_order(self, df: pd.DataFrame, max_storage_gb: float):
        """
        Display the dataset scheduling order and storage analysis.

        Args:
            df: DataFrame with dataset information
            max_storage_gb: Maximum storage limit
        """
        df_sorted = df.sort_values('SafeStorageGB', ascending=False)

        print(f"\n{'='*60}")
        print(f"DATASET SCHEDULING ORDER (Largest First)")
        print(f"{'='*60}")
        print(f"Batch system: {self.name.upper()}")
        print(f"Maximum storage limit: {max_storage_gb:.2f} GB")
        print(f"Total datasets: {len(df_sorted)}")

        total_storage = df_sorted['SafeStorageGB'].sum()
        print(f"Total storage needed: {total_storage:.2f} GB")

        # Show largest 10 datasets
        print(f"\nLargest 10 datasets (will be prioritized first):")
        for i, row in df_sorted.head(10).iterrows():
            print(f"  {i+1:2d}. {row['Accession']:12s} - {row['SafeStorageGB']:6.2f} GB")

        if len(df_sorted) > 10:
            print(f"  ... and {len(df_sorted)-10} more datasets")

        # Show smallest 5 datasets
        if len(df_sorted) > 5:
            print(f"\nSmallest 5 datasets (will be processed last):")
            for i, row in df_sorted.tail(5).iterrows():
                pos = len(df_sorted) - (len(df_sorted.tail(5)) - list(df_sorted.tail(5).index).index(i))
                print(f"  {pos:2d}. {row['Accession']:12s} - {row['SafeStorageGB']:6.2f} GB")

        # Storage analysis
        datasets_that_fit = df_sorted[df_sorted['SafeStorageGB'] <= max_storage_gb]
        datasets_too_large = df_sorted[df_sorted['SafeStorageGB'] > max_storage_gb]

        print(f"\nStorage Analysis:")
        print(f"  Datasets that fit individually: {len(datasets_that_fit)}/{len(df_sorted)}")
        if len(datasets_too_large) > 0:
            print(f"  Datasets too large for storage limit:")
            for _, row in datasets_too_large.iterrows():
                print(f"    {row['Accession']:12s} - {row['SafeStorageGB']:6.2f} GB (exceeds {max_storage_gb:.2f} GB limit)")

        print(f"{'='*60}\n")

    def check_environment_requirements(self) -> List[str]:
        """
        Check if all required environment variables and tools are available.

        Note: Environment variable validation is now handled at the CLI level.

        Returns:
            List of missing requirements (empty if all satisfied)
        """
        # Environment variables are now validated in CLI before reaching this point
        return []

    @property
    @abstractmethod
    def default_parameters(self) -> Dict[str, Any]:
        """
        Get default parameters for this batch system.

        Returns:
            Dictionary of default parameters
        """
        pass

    def __str__(self) -> str:
        """String representation of the batch processor."""
        return f"{self.name.capitalize()} Batch Processor"

    def __repr__(self) -> str:
        """Developer representation of the batch processor."""
        return f"<{self.__class__.__name__}()>"

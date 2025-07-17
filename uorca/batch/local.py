"""
Local Batch Processor

This module provides local parallel processing capabilities for UORCA datasets
using Docker containers. It uses multiprocessing to run jobs in parallel
while respecting resource constraints.
"""

import json
import os
import psutil
import subprocess
import time
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, Future, as_completed
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
import multiprocessing as mp
from dotenv import load_dotenv, find_dotenv

from .base import BatchProcessor


def run_single_dataset_local(accession: str, output_dir: str, resource_dir: str,
                           cleanup: bool, status_dir: str, docker_image: str = "kevingchen/uorca:0.1.0") -> Dict[str, Any]:
    """
    Run a single dataset analysis locally using Docker.

    This function is designed to be run in a separate process.

    Args:
        accession: Dataset accession ID
        output_dir: Output directory on host
        resource_dir: Resource directory for Kallisto indices on host
        cleanup: Whether to cleanup intermediate files
        status_dir: Directory for status tracking on host
        docker_image: Docker image to use

    Returns:
        Dictionary with job results
    """
    import subprocess
    import json
    from datetime import datetime
    from pathlib import Path

    process_id = os.getpid()
    start_time = datetime.now()

    # Update status to running
    status_file = Path(status_dir) / f"{accession}_status.json"
    if status_file.exists():
        with open(status_file, 'r') as f:
            status = json.load(f)
        status.update({
            'state': 'running',
            'started_time': start_time.isoformat(),
            'process_id': process_id
        })
        with open(status_file, 'w') as f:
            json.dump(status, f, indent=2)

    try:
        # Setup paths
        output_path = Path(output_dir).absolute()
        resource_path = Path(resource_dir).absolute()

        # Create directories if they don't exist
        output_path.mkdir(parents=True, exist_ok=True)
        resource_path.mkdir(parents=True, exist_ok=True)

        # Docker paths (inside container)
        docker_output_dir = "/app/output"
        docker_resource_dir = "/app/resources"

        # Build Docker command
        cmd = [
            'docker', 'run', '--rm',
            # Mount volumes
            '-v', f'{output_path}:{docker_output_dir}',
            '-v', f'{resource_path}:{docker_resource_dir}',
            # Pass environment variables
            '-e', f'ENTREZ_EMAIL={os.getenv("ENTREZ_EMAIL", "")}',
            '-e', f'OPENAI_API_KEY={os.getenv("OPENAI_API_KEY", "")}',
            '-e', f'ENTREZ_API_KEY={os.getenv("ENTREZ_API_KEY", "")}',
            # Use the specified Docker image
            docker_image,
            # Command to run inside container
            'uv', 'run', 'python', 'main_workflow/master.py',
            '--accession', accession,
            '--output_dir', docker_output_dir,
            '--resource_dir', docker_resource_dir
        ]

        if cleanup:
            cmd.append('--cleanup')

        # Run the analysis
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=21600  # 6 hours timeout
        )

        end_time = datetime.now()
        success = result.returncode == 0

        # Update final status
        final_status = {
            'state': 'completed' if success else 'failed',
            'completed_time': end_time.isoformat(),
            'runtime_seconds': (end_time - start_time).total_seconds(),
            'return_code': result.returncode,
            'success': success,
            'stdout': result.stdout[-2000:] if result.stdout else '',  # Last 2000 chars
            'stderr': result.stderr[-2000:] if result.stderr else ''   # Last 2000 chars
        }

        if status_file.exists():
            with open(status_file, 'r') as f:
                status = json.load(f)
            status.update(final_status)
            with open(status_file, 'w') as f:
                json.dump(status, f, indent=2)

        return {
            'accession': accession,
            'success': success,
            'return_code': result.returncode,
            'runtime': (end_time - start_time).total_seconds(),
            'stdout': result.stdout,
            'stderr': result.stderr
        }

    except subprocess.TimeoutExpired:
        # Handle timeout
        final_status = {
            'state': 'timeout',
            'completed_time': datetime.now().isoformat(),
            'error': 'Job exceeded 6 hour timeout',
            'return_code': -1
        }

        if status_file.exists():
            with open(status_file, 'r') as f:
                status = json.load(f)
            status.update(final_status)
            with open(status_file, 'w') as f:
                json.dump(status, f, indent=2)

        return {
            'accession': accession,
            'success': False,
            'error': 'timeout',
            'runtime': 21600
        }

    except Exception as e:
        # Handle other errors
        final_status = {
            'state': 'failed',
            'completed_time': datetime.now().isoformat(),
            'error': str(e),
            'return_code': -1
        }

        if status_file.exists():
            with open(status_file, 'r') as f:
                status = json.load(f)
            status.update(final_status)
            with open(status_file, 'w') as f:
                json.dump(status, f, indent=2)

        return {
            'accession': accession,
            'success': False,
            'error': str(e),
            'runtime': 0
        }


class LocalBatchProcessor(BatchProcessor):
    """
    Local parallel batch processor for UORCA datasets using Docker.

    This processor handles parallel execution of dataset analyses on the local machine
    using Docker containers, with resource management and job tracking.
    """

    def __init__(self):
        """Initialize the local batch processor."""
        super().__init__()
        self.executor = None
        self.futures = {}  # Map accession -> Future
        self.job_counter = 0

    def validate_csv_file(self, csv_file: str) -> pd.DataFrame:
        """
        Validate and load the CSV file containing dataset information.
        Override base class to handle DatasetSizeGB column.

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

        # Handle storage size column - prefer DatasetSizeGB, fallback to DatasetSizeBytes
        if 'DatasetSizeGB' in df.columns:
            df['SafeStorageGB'] = df['DatasetSizeGB'] * 3.0  # Apply safety factor
        elif 'DatasetSizeBytes' in df.columns:
            df['SafeStorageGB'] = (df['DatasetSizeBytes'] / (1024**3)) * 3.0
            print("Using DatasetSizeBytes column for storage calculations")
        else:
            df['SafeStorageGB'] = 30.0  # Default safety value of 30GB per dataset
            print("Warning: Neither DatasetSizeGB nor DatasetSizeBytes column found, using default value of 30GB per dataset")

        return df

    @property
    def default_parameters(self) -> Dict[str, Any]:
        """Get default parameters for local batch processing."""
        # Use 75% of available CPU cores, minimum 1
        max_workers = max(1, int(mp.cpu_count() * 0.75))

        # Use 75% of available memory for storage limit
        total_memory_gb = psutil.virtual_memory().total / (1024**3)
        max_storage_gb = max(50, int(total_memory_gb * 0.75))

        return {
            'max_workers': max_workers,
            'max_storage_gb': max_storage_gb,
            'cleanup': True,  # Default to cleanup
            'resource_dir': './data/kallisto_indices/',
            'check_interval': 10,
            'timeout_hours': 6,
            'docker_image': 'kevingchen/uorca:0.1.0'
        }

    def load_environment_variables(self):
        """Load environment variables from .env file."""
        # Try to find and load .env file
        project_root = Path(__file__).resolve().parent.parent.parent
        env_file = project_root / ".env"

        if env_file.exists():
            load_dotenv(env_file)
            print(f"Loaded environment variables from {env_file}")
        else:
            # Try to find .env in current directory or parent directories
            load_dotenv(find_dotenv())

    def check_environment_requirements(self) -> List[str]:
        """
        Check if all required environment variables and tools are available.

        Returns:
            List of missing requirements (empty if all satisfied)
        """
        # Load environment variables first
        self.load_environment_variables()

        missing = []

        # Check for Docker
        try:
            result = subprocess.run(['docker', '--version'],
                                  capture_output=True, text=True, timeout=10)
            if result.returncode != 0:
                missing.append("Docker (not available or not running)")
        except (subprocess.TimeoutExpired, FileNotFoundError):
            missing.append("Docker (not installed or not in PATH)")

        # Check for required environment variables
        if not os.getenv('ENTREZ_EMAIL'):
            missing.append("Environment variable: ENTREZ_EMAIL")
            print("\nMissing required environment variables!")
            print("Please create a .env file in your project root with the following variables:")
            print("")
            print("  ENTREZ_EMAIL=your_value_here  # Required for NCBI API access")
            print("")
            print("Example .env file:")
            print("  ENTREZ_EMAIL=your.email@example.com")
            print("  OPENAI_API_KEY=sk-your-openai-api-key  # Optional, for AI features")
            print("  ENTREZ_API_KEY=your_ncbi_api_key  # Optional, for higher rate limits")

        # Check for optional environment variables (warn but don't fail)
        if not os.getenv('OPENAI_API_KEY'):
            print("Note: OPENAI_API_KEY not set. AI-powered features will be disabled.")

        if not os.getenv('ENTREZ_API_KEY'):
            print("Note: ENTREZ_API_KEY not set. Using default NCBI rate limits.")

        return missing

    def submit_datasets(self, csv_file: str, output_dir: str, **kwargs) -> int:
        """
        Submit datasets for local parallel processing using Docker.

        Args:
            csv_file: Path to CSV file containing dataset information
            output_dir: Output directory for results
            **kwargs: Additional local processing parameters

        Returns:
            Number of jobs submitted
        """
        # Merge default parameters with provided kwargs
        params = {**self.default_parameters, **kwargs}

        # Validate inputs
        df = self.validate_csv_file(csv_file)
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        # Check environment requirements
        missing_reqs = self.check_environment_requirements()
        if missing_reqs:
            raise EnvironmentError(f"Missing requirements: {missing_reqs}")

        # Test Docker image availability
        try:
            print(f"Checking Docker image: {params['docker_image']}")
            result = subprocess.run([
                'docker', 'image', 'inspect', params['docker_image']
            ], capture_output=True, text=True, timeout=30)

            if result.returncode != 0:
                print(f"Pulling Docker image: {params['docker_image']}")
                pull_result = subprocess.run([
                    'docker', 'pull', params['docker_image']
                ], timeout=300)  # 5 minute timeout for pull

                if pull_result.returncode != 0:
                    raise EnvironmentError(f"Failed to pull Docker image: {params['docker_image']}")
                else:
                    print(f"Successfully pulled Docker image: {params['docker_image']}")
            else:
                print(f"Docker image {params['docker_image']} is available")

        except subprocess.TimeoutExpired:
            raise EnvironmentError("Docker operations timed out. Check Docker daemon status.")
        except Exception as e:
            raise EnvironmentError(f"Docker error: {e}")

        # Setup job tracking
        status_dir = self.setup_job_tracking(output_dir)
        logs_dir = output_path / "logs"
        logs_dir.mkdir(exist_ok=True)

        # Sort datasets by size (largest first for better resource utilization)
        df_sorted = df.sort_values('SafeStorageGB', ascending=False)

        # Display scheduling information
        self.display_scheduling_order(df_sorted, params['max_storage_gb'])

        # Filter datasets that fit within storage constraints
        datasets_to_process = []
        total_storage_needed = 0

        for _, row in df_sorted.iterrows():
            accession = row['Accession']
            storage_needed = row['SafeStorageGB']

            if storage_needed > params['max_storage_gb']:
                print(f"Skipping {accession}: exceeds storage limit ({storage_needed:.2f} > {params['max_storage_gb']:.2f} GB)")
                continue

            datasets_to_process.append(row)
            total_storage_needed += storage_needed

        if not datasets_to_process:
            print("No datasets to process after storage filtering.")
            return 0

        print(f"\nStarting local parallel processing...")
        print(f"   Docker image: {params['docker_image']}")
        print(f"   Max workers: {params['max_workers']}")
        print(f"   Total datasets: {len(datasets_to_process)}")
        print(f"   Total storage needed: {total_storage_needed:.2f} GB")
        print(f"   Storage limit: {params['max_storage_gb']:.2f} GB")
        print(f"   Cleanup enabled: {params['cleanup']}")

        # Initialize process pool executor
        self.executor = ProcessPoolExecutor(max_workers=params['max_workers'])
        jobs_submitted = 0

        try:
            # Submit all jobs
            for row in datasets_to_process:
                accession = row['Accession']
                storage_gb = row['SafeStorageGB']

                # Create job status file
                job_id = f"local_{self.job_counter:04d}"
                self.job_counter += 1

                self.create_job_status_file(
                    status_dir=status_dir,
                    accession=accession,
                    job_id=job_id,
                    storage_gb=storage_gb,
                    max_workers=params['max_workers'],
                    docker_image=params['docker_image']
                )

                # Submit job
                future = self.executor.submit(
                    run_single_dataset_local,
                    accession=accession,
                    output_dir=output_dir,
                    resource_dir=params['resource_dir'],
                    cleanup=params['cleanup'],
                    status_dir=str(status_dir),
                    docker_image=params['docker_image']
                )

                self.futures[accession] = future
                jobs_submitted += 1

                print(f"  Submitted {accession} (job {job_id})")

            print(f"\nAll {jobs_submitted} jobs submitted")
            print(f"Job logs will be written to: {logs_dir}")
            print(f"Job status tracking: {status_dir}")

            # Monitor job progress
            self._monitor_jobs(params['check_interval'])

        finally:
            # Clean up executor
            if self.executor:
                self.executor.shutdown(wait=False)
                self.executor = None

        return jobs_submitted

    def _monitor_jobs(self, check_interval: int):
        """
        Monitor running jobs and display progress.

        Args:
            check_interval: Check interval in seconds
        """
        print(f"\nMonitoring job progress (checking every {check_interval}s)...")
        print("   Press Ctrl+C to stop monitoring (jobs will continue running)")

        try:
            while self.futures:
                completed_jobs = []

                for accession, future in self.futures.items():
                    if future.done():
                        completed_jobs.append(accession)

                        try:
                            result = future.result()
                            if result['success']:
                                print(f"  {accession} completed successfully ({result['runtime']:.0f}s)")
                            else:
                                error_msg = result.get('error', 'Unknown error')
                                if 'stderr' in result and result['stderr']:
                                    print(f"  {accession} failed: {error_msg}")
                                    print(f"      stderr: {result['stderr'][:200]}...")  # First 200 chars
                                else:
                                    print(f"  {accession} failed: {error_msg}")
                        except Exception as e:
                            print(f"  {accession} failed with exception: {e}")

                # Remove completed jobs
                for accession in completed_jobs:
                    del self.futures[accession]

                if not self.futures:
                    break

                # Show status
                running_count = len(self.futures)
                print(f"  {running_count} jobs still running...")

                time.sleep(check_interval)

        except KeyboardInterrupt:
            print(f"\nMonitoring stopped. {len(self.futures)} jobs still running in background.")
            print("   Check status files for progress updates.")

    def check_job_status(self, job_id: str) -> Dict[str, Any]:
        """
        Check the status of a specific local job.

        Args:
            job_id: Local job ID

        Returns:
            Dictionary containing job status information
        """
        # For local jobs, we check the futures dict and process status
        return {
            'job_id': job_id,
            'state': 'UNKNOWN',
            'system': 'local',
            'found': False,
            'note': 'Use job status files for local job tracking'
        }

    def cancel_job(self, job_id: str) -> bool:
        """
        Cancel a local job.

        Args:
            job_id: Local job ID

        Returns:
            True if job was successfully cancelled, False otherwise
        """
        # Try to cancel futures if they exist
        cancelled_count = 0

        for accession, future in list(self.futures.items()):
            if future and not future.done():
                if future.cancel():
                    cancelled_count += 1
                    del self.futures[accession]
                    print(f"Cancelled job for {accession}")

        return cancelled_count > 0

    def get_available_resources(self) -> Dict[str, Any]:
        """
        Get information about available local resources.

        Returns:
            Dictionary containing resource information
        """
        # Get system information
        cpu_count = mp.cpu_count()
        memory = psutil.virtual_memory()
        disk = psutil.disk_usage('/')

        # Get load average if available (Unix systems)
        load_avg = None
        try:
            load_avg = os.getloadavg()
        except AttributeError:
            pass  # Windows doesn't have getloadavg

        # Check Docker status
        docker_status = "unknown"
        try:
            result = subprocess.run(['docker', 'info'],
                                  capture_output=True, text=True, timeout=10)
            docker_status = "running" if result.returncode == 0 else "error"
        except:
            docker_status = "not available"

        return {
            'system': 'local',
            'cpu_cores': cpu_count,
            'cpu_logical': psutil.cpu_count(logical=True),
            'memory_total_gb': memory.total / (1024**3),
            'memory_available_gb': memory.available / (1024**3),
            'memory_percent_used': memory.percent,
            'disk_total_gb': disk.total / (1024**3),
            'disk_free_gb': disk.free / (1024**3),
            'disk_percent_used': (disk.used / disk.total) * 100,
            'load_average': load_avg,
            'docker_status': docker_status,
            'active_jobs': len(self.futures) if self.futures else 0
        }

    def list_active_jobs(self) -> List[Dict[str, Any]]:
        """
        List currently active local jobs.

        Returns:
            List of job information dictionaries
        """
        jobs = []

        if self.futures:
            for accession, future in self.futures.items():
                if future and not future.done():
                    jobs.append({
                        'accession': accession,
                        'state': 'running',
                        'system': 'local'
                    })
                elif future and future.done():
                    try:
                        result = future.result()
                        jobs.append({
                            'accession': accession,
                            'state': 'completed' if result['success'] else 'failed',
                            'system': 'local'
                        })
                    except Exception:
                        jobs.append({
                            'accession': accession,
                            'state': 'failed',
                            'system': 'local'
                        })

        return jobs

    def cleanup(self):
        """Clean up resources and shutdown executor."""
        if self.executor:
            self.executor.shutdown(wait=True)
            self.executor = None
        self.futures.clear()

    def __del__(self):
        """Destructor to ensure cleanup."""
        self.cleanup()

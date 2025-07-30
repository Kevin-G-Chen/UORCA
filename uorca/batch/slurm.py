"""
SLURM Batch Processor

This module provides SLURM-specific batch processing capabilities for UORCA.
It adapts the existing submit_datasets.py functionality into the new batch processing architecture.
"""

import json
import logging
import os
import subprocess
import time
import pandas as pd
import yaml
from datetime import datetime
from jinja2 import Environment, FileSystemLoader
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple

from .base import BatchProcessor


class SlurmBatchProcessor(BatchProcessor):
    """
    SLURM-specific batch processor for UORCA datasets.

    This processor handles job submission, monitoring, and resource management
    specifically for SLURM-based HPC systems.
    """

    def __init__(self, config_file: Optional[str] = None):
        """Initialize the SLURM batch processor."""
        self.config_file = config_file
        self._previous_job_states = {}  # Track job states for change detection
        super().__init__()
        self._validate_slurm_environment()

    def _validate_slurm_environment(self):
        """Validate that SLURM tools are available."""
        required_commands = ['sbatch', 'squeue', 'scancel', 'sinfo']
        missing_commands = []

        for cmd in required_commands:
            if not subprocess.run(['which', cmd], capture_output=True).returncode == 0:
                missing_commands.append(cmd)

        if missing_commands:
            raise EnvironmentError(f"Missing SLURM commands: {missing_commands}")

    def _load_config(self) -> Dict[str, Any]:
        """Load configuration from YAML file if provided."""
        if not self.config_file:
            return {}

        try:
            config_path = Path(self.config_file)
            if not config_path.exists():
                logging.warning(f"Config file not found: {self.config_file}, using defaults")
                return {}

            with open(config_path, 'r') as f:
                config = yaml.safe_load(f)

            logging.info(f"Loaded SLURM configuration from {self.config_file}")
            return config or {}
        except Exception as e:
            logging.warning(f"Failed to load config file {self.config_file}: {e}, using defaults")
            return {}

    @property
    def default_parameters(self) -> Dict[str, Any]:
        """Get default parameters for SLURM batch processing, merging config file and defaults."""
        # Load config file
        config = self._load_config()

        # Extract relevant sections
        slurm_config = config.get('slurm', {})
        resource_config = config.get('resource_management', {})
        container_config = config.get('container', {})

        # Merge with defaults (CLI args will override these later)
        defaults = {
            'max_parallel': resource_config.get('max_parallel', 10),
            'max_storage_gb': resource_config.get('max_storage_gb', 500),
            'cleanup': True,
            'resource_dir': './data/kallisto_indices/',
            'partition': slurm_config.get('partition', 'tki_agpdev'),
            'constraint': slurm_config.get('constraint', 'icx'),
            'cpus_per_task': slurm_config.get('cpus_per_task', 8),
            'memory': slurm_config.get('memory', '16G'),
            'time_limit': slurm_config.get('time_limit', '6:00:00'),
            'check_interval': resource_config.get('check_interval', 300),
            'container_engine': container_config.get('engine', 'apptainer'),
            'apptainer_image': container_config.get('apptainer_image', '/data/tki_agpdev/kevin/phd/aim1/UORCA/scratch/container_testing/uorca_0.1.0.sif'),
            'docker_image': container_config.get('docker_image', 'kevingchen/uorca:0.1.0')
        }

        return defaults

    def submit_datasets(self, input_path: str, output_dir: str, **kwargs) -> int:
        """
        Submit SLURM batch jobs for processing datasets.

        Args:
            input_path: Either CSV file with dataset information or directory containing CSV and metadata
            output_dir: Output directory for results
            **kwargs: Additional SLURM-specific parameters

        Returns:
            Number of jobs submitted
        """
        # Merge default parameters with provided kwargs
        params = {**self.default_parameters, **kwargs}

        # Validate inputs and extract research question
        df, research_question = self.validate_csv_file(input_path)
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        # Save research question if available
        self.save_research_question(output_dir, research_question)

        # Environment requirements are now checked at CLI level

        # Setup job tracking
        status_dir = self.setup_job_tracking(output_dir)
        logs_dir = output_path / "logs"
        logs_dir.mkdir(exist_ok=True)

        # Sort datasets by size (largest first for better resource utilization)
        df_sorted = df.sort_values('SafeStorageGB', ascending=False)

        # Display scheduling information
        self.display_scheduling_order(df_sorted, params['max_storage_gb'])

        # Display SLURM configuration being used
        self._display_slurm_configuration(params)

        # Submit jobs with storage-aware scheduling
        jobs_submitted = 0

        # Create a list of datasets to process (excluding those that exceed storage limit)
        processable_datasets = []
        oversized_datasets = []

        for _, row in df_sorted.iterrows():
            if row['SafeStorageGB'] > params['max_storage_gb']:
                oversized_datasets.append(row)
            else:
                processable_datasets.append(row)

        # Report oversized datasets
        if oversized_datasets:
            print(f"\n⚠️  Skipping {len(oversized_datasets)} dataset(s) that exceed storage limit ({params['max_storage_gb']:.2f} GB):")
            for row in oversized_datasets:
                print(f"    {row['Accession']:12s} - {row['SafeStorageGB']:6.2f} GB")

        # Process datasets with intelligent scheduling
        remaining_datasets = processable_datasets.copy()

        while remaining_datasets:
            # Find the largest dataset that can fit current resource constraints
            submitted_this_round = False

            for i, row in enumerate(remaining_datasets):
                accession = row['Accession']
                storage_needed = row['SafeStorageGB']

                print(f"Processing {accession} (requires {storage_needed:.2f} GB)...")

                # Wait for available resources if needed
                self._wait_for_resources(output_dir, storage_needed, params['max_parallel'], params['max_storage_gb'], params['check_interval'])

                # Submit the job
                job_id = self._submit_single_dataset(
                    accession=accession,
                    output_dir=output_dir,
                    logs_dir=logs_dir,
                    status_dir=status_dir,
                    storage_gb=storage_needed,
                    **params
                )

                if job_id:
                    print(f"  ✅ Submitted {accession} as job {job_id}")
                    jobs_submitted += 1
                    remaining_datasets.pop(i)
                    submitted_this_round = True
                    break
                else:
                    print(f"  ❌ Failed to submit {accession}")
                    remaining_datasets.pop(i)
                    break

            # If no dataset could be submitted this round, break to avoid infinite loop
            if not submitted_this_round and remaining_datasets:
                print(f"  ⚠️  No remaining datasets can be submitted with current constraints")
                break

            # Brief pause between submissions
            if submitted_this_round:
                time.sleep(1)

        print(f"\n🎉 Successfully submitted {jobs_submitted} jobs to SLURM")
        print(f"📁 Job logs will be written to: {logs_dir}")
        print(f"📊 Job status tracking: {status_dir}")
        print(f"\nMonitor progress with: squeue -u $(whoami)")

        return jobs_submitted

    def _submit_single_dataset(self, accession: str, output_dir: str, logs_dir: Path,
                             status_dir: Path, storage_gb: float, **params) -> Optional[str]:
        """
        Submit a single dataset processing job to SLURM.

        Args:
            accession: Dataset accession ID
            output_dir: Output directory
            logs_dir: Directory for log files
            status_dir: Directory for status files
            storage_gb: Storage requirement in GB
            **params: SLURM parameters

        Returns:
            Job ID if successful, None otherwise
        """
        try:
            # Generate SLURM script from template
            script_content = self._generate_slurm_script(
                accession=accession,
                output_dir=output_dir,
                logs_dir=logs_dir,
                **params
            )

            # Write temporary script file
            script_file = Path(output_dir) / f"temp_submit_{accession}.sbatch"
            with open(script_file, 'w') as f:
                f.write(script_content)

            # Submit job
            result = subprocess.run(
                ['sbatch', str(script_file)],
                capture_output=True,
                text=True,
                timeout=30
            )

            # Clean up script file immediately
            script_file.unlink(missing_ok=True)

            if result.returncode == 0:
                # Parse job ID from sbatch output
                job_id = result.stdout.strip().split()[-1]

                # Create job status file
                self.create_job_status_file(
                    status_dir=status_dir,
                    accession=accession,
                    job_id=job_id,
                    storage_gb=storage_gb,
                    partition=params.get('partition'),
                    cpus=params.get('cpus_per_task'),
                    memory=params.get('memory')
                )

                return job_id
            else:
                print(f"    Error submitting {accession}: {result.stderr}")
                return None

        except Exception as e:
            print(f"    Exception submitting {accession}: {e}")
            return None

    def _generate_slurm_script(self, accession: str, output_dir: str, logs_dir: Path, **params) -> str:
        """
        Generate SLURM script content from Jinja2 template.

        Args:
            accession: Dataset accession
            output_dir: Output directory
            logs_dir: Log files directory
            **params: Template parameters

        Returns:
            Generated script content
        """
        # Find template directory
        current_dir = Path(__file__).parent
        project_root = current_dir.parent.parent
        template_dir = project_root / "main_workflow" / "run_helpers"

        if not template_dir.exists():
            raise FileNotFoundError(f"Template directory not found: {template_dir}")

        # Setup Jinja2 environment
        env = Environment(loader=FileSystemLoader(str(template_dir)))
        template = env.get_template('run_single_dataset.sbatch.j2')

        # Determine container engine and image
        container_engine = params.get('container_engine', 'apptainer')
        if container_engine == 'apptainer':
            container_image = params.get('apptainer_image', '/data/tki_agpdev/kevin/phd/aim1/UORCA/scratch/container_testing/uorca_0.1.0.sif')
        else:
            container_image = params.get('docker_image', 'kevingchen/uorca:0.1.0')

        # Render template with parameters
        script_content = template.render(
            accession=accession,
            output_dir=output_dir,
            logs_dir=logs_dir,
            project_root=project_root,
            resource_dir=params.get('resource_dir', './data/kallisto_indices/'),
            cleanup=params.get('cleanup', True),
            partition=params.get('partition', 'tki_agpdev'),
            constraint=params.get('constraint', 'icx'),
            cpus_per_task=params.get('cpus_per_task', 8),
            memory=params.get('memory', '16G'),
            time_limit=params.get('time_limit', '6:00:00'),
            container_engine=container_engine,
            container_image=container_image
        )

        return script_content

    def _display_slurm_configuration(self, params: Dict[str, Any]):
        """
        Display the SLURM configuration parameters being used for job submission.

        Args:
            params: Dictionary of parameters being used
        """
        print(f"\n{'='*60}")
        print(f"SLURM CONFIGURATION")
        print(f"{'='*60}")

        # SLURM job parameters
        print(f"Job Parameters:")
        print(f"  Partition:        {params.get('partition', 'N/A')}")
        print(f"  Constraint:       {params.get('constraint', 'N/A')}")
        print(f"  CPUs per task:    {params.get('cpus_per_task', 'N/A')}")
        print(f"  Memory per job:   {params.get('memory', 'N/A')}")
        print(f"  Time limit:       {params.get('time_limit', 'N/A')}")

        # Resource management
        print(f"\nResource Management:")
        print(f"  Max parallel jobs:    {params.get('max_parallel', 'N/A')}")
        print(f"  Max storage (GB):     {params.get('max_storage_gb', 'N/A')}")
        print(f"  Check interval (s):   {params.get('check_interval', 'N/A')}")

        # Container settings
        print(f"\nContainer Settings:")
        print(f"  Engine:           {params.get('container_engine', 'N/A')}")
        if params.get('container_engine') == 'apptainer':
            print(f"  Image:            {params.get('apptainer_image', 'N/A')}")
        else:
            print(f"  Image:            {params.get('docker_image', 'N/A')}")

        # Workflow settings
        print(f"\nWorkflow Settings:")
        print(f"  Resource dir:     {params.get('resource_dir', 'N/A')}")
        print(f"  Cleanup enabled:  {params.get('cleanup', 'N/A')}")

        # Configuration source
        config_source = "Configuration file" if self.config_file else "Built-in defaults"
        if self.config_file:
            print(f"\nConfiguration loaded from: {self.config_file}")
        else:
            print(f"\nUsing built-in default configuration")

        print(f"{'='*60}\n")

    def _wait_for_resources(self, output_dir: str, storage_needed: float,
                          max_parallel: int, max_storage_gb: float, check_interval: int):
        """
        Wait until sufficient resources are available for the next job.

        Args:
            output_dir: Output directory for checking running jobs
            storage_needed: Storage required for next job (GB)
            max_parallel: Maximum parallel jobs
            max_storage_gb: Maximum total storage (GB)
            check_interval: Check interval in seconds
        """
        while True:
            # Check running jobs
            running_jobs = self._count_running_jobs()
            running_storage = self.get_running_jobs_storage(output_dir)

            # Check if we can submit the next job
            storage_available = (running_storage + storage_needed) <= max_storage_gb
            parallel_available = running_jobs < max_parallel

            if storage_available and parallel_available:
                break

            # Show waiting message
            print(f"  ⏳ Waiting for resources...")
            print(f"    Running jobs: {running_jobs}/{max_parallel}")
            print(f"    Storage usage: {running_storage:.1f}/{max_storage_gb:.1f} GB")
            print(f"    Need: {storage_needed:.1f} GB")

            time.sleep(check_interval)

    def _count_running_jobs(self) -> int:
        """
        Count the number of currently running UORCA jobs.

        Returns:
            Number of running jobs
        """
        try:
            result = subprocess.run(
                ['squeue', '-u', os.getenv('USER', 'unknown'), '-h', '-o', '%j'],
                capture_output=True,
                text=True,
                timeout=10
            )

            if result.returncode == 0:
                # Count jobs with UORCA-related names
                job_names = result.stdout.strip().split('\n')
                uorca_jobs = [name for name in job_names if 'run_' in name or 'uorca' in name.lower()]
                return len(uorca_jobs)
            else:
                return 0

        except (subprocess.TimeoutExpired, FileNotFoundError):
            return 0

    def monitor_job_progress(self, output_dir: str) -> Dict[str, int]:
        """
        Monitor job progress and only print updates when status changes.

        Args:
            output_dir: Output directory containing job status files

        Returns:
            Dictionary with counts of jobs in each state
        """
        status_dir = Path(output_dir) / "job_status"

        if not status_dir.exists():
            return {'running': 0, 'completed': 0, 'failed': 0}

        current_states = {}
        state_counts = {'running': 0, 'completed': 0, 'failed': 0}

        # Read current job statuses
        for status_file in status_dir.glob("*_status.json"):
            try:
                with open(status_file, 'r') as f:
                    status = json.load(f)

                accession = status.get('accession', 'unknown')
                job_state = status.get('state', 'unknown')
                job_id = status.get('job_id', 'unknown')

                current_states[accession] = {
                    'state': job_state,
                    'job_id': job_id,
                    'storage_gb': status.get('storage_gb', 0)
                }

                # Count states
                if job_state in state_counts:
                    state_counts[job_state] += 1

            except Exception as e:
                logging.warning(f"Failed to read status file {status_file}: {e}")
                continue

        # Check for changes and print updates
        changes_detected = False
        for accession, current_info in current_states.items():
            previous_info = self._previous_job_states.get(accession)

            if previous_info is None:
                # New job
                print(f"  ▶ {accession} (Job {current_info['job_id']}) - {current_info['state']}")
                changes_detected = True
            elif previous_info['state'] != current_info['state']:
                # State change
                old_state = previous_info['state']
                new_state = current_info['state']

                if new_state == 'completed':
                    print(f"  ✅ {accession} - COMPLETED (was {old_state})")
                elif new_state == 'failed':
                    print(f"  ❌ {accession} - FAILED (was {old_state})")
                elif new_state == 'running':
                    print(f"  🔄 {accession} - STARTED (was {old_state})")
                else:
                    print(f"  ℹ️  {accession} - {old_state} → {new_state}")

                changes_detected = True

        # Check for jobs that disappeared (shouldn't happen normally)
        for accession in self._previous_job_states:
            if accession not in current_states:
                print(f"  ⚠️  {accession} - Status file removed")
                changes_detected = True

        # Update tracking
        self._previous_job_states = current_states.copy()

        # Print summary only if there were changes or this is the first check
        if changes_detected or not hasattr(self, '_monitoring_initialized'):
            running = state_counts['running']
            completed = state_counts['completed']
            failed = state_counts['failed']
            total = running + completed + failed

            if total > 0:
                print(f"  📊 Status: {running} running, {completed} completed, {failed} failed (total: {total})")

            self._monitoring_initialized = True

        return state_counts

    def check_job_status(self, job_id: str) -> Dict[str, Any]:
        """
        Check the status of a specific SLURM job.

        Args:
            job_id: SLURM job ID

        Returns:
            Dictionary containing job status information
        """
        try:
            result = subprocess.run(
                ['squeue', '-j', job_id, '-h', '-o', '%T,%S,%M,%L'],
                capture_output=True,
                text=True,
                timeout=10
            )

            if result.returncode == 0 and result.stdout.strip():
                state, state_code, time_used, time_left = result.stdout.strip().split(',')
                return {
                    'job_id': job_id,
                    'state': state,
                    'state_code': state_code,
                    'time_used': time_used,
                    'time_left': time_left,
                    'found': True
                }
            else:
                # Job not found in queue (completed or failed)
                return {
                    'job_id': job_id,
                    'state': 'COMPLETED_OR_FAILED',
                    'found': False
                }

        except (subprocess.TimeoutExpired, FileNotFoundError) as e:
            return {
                'job_id': job_id,
                'state': 'UNKNOWN',
                'error': str(e),
                'found': False
            }

    def cancel_job(self, job_id: str) -> bool:
        """
        Cancel a SLURM job.

        Args:
            job_id: SLURM job ID

        Returns:
            True if job was successfully cancelled, False otherwise
        """
        try:
            result = subprocess.run(
                ['scancel', job_id],
                capture_output=True,
                text=True,
                timeout=10
            )

            return result.returncode == 0

        except (subprocess.TimeoutExpired, FileNotFoundError):
            return False

    def get_available_resources(self) -> Dict[str, Any]:
        """
        Get information about available SLURM resources.

        Returns:
            Dictionary containing resource information
        """
        resources = {
            'system': 'slurm',
            'partitions': [],
            'total_nodes': 0,
            'available_nodes': 0,
            'total_cpus': 0,
            'available_cpus': 0
        }

        try:
            # Get partition information
            result = subprocess.run(
                ['sinfo', '-h', '-o', '%P,%N,%T,%c,%m'],
                capture_output=True,
                text=True,
                timeout=15
            )

            if result.returncode == 0:
                for line in result.stdout.strip().split('\n'):
                    if line:
                        parts = line.split(',')
                        if len(parts) >= 5:
                            partition = parts[0].rstrip('*')
                            nodes = parts[1]
                            state = parts[2]
                            cpus = parts[3]
                            memory = parts[4]

                            resources['partitions'].append({
                                'name': partition,
                                'nodes': nodes,
                                'state': state,
                                'cpus': cpus,
                                'memory': memory
                            })

            # Get queue information
            queue_result = subprocess.run(
                ['squeue', '-h', '-o', '%T'],
                capture_output=True,
                text=True,
                timeout=10
            )

            if queue_result.returncode == 0:
                states = queue_result.stdout.strip().split('\n')
                resources['running_jobs'] = len([s for s in states if s in ['RUNNING', 'R']])
                resources['pending_jobs'] = len([s for s in states if s in ['PENDING', 'PD']])

        except (subprocess.TimeoutExpired, FileNotFoundError):
            pass

        return resources

    def list_user_jobs(self) -> List[Dict[str, Any]]:
        """
        List all SLURM jobs for the current user.

        Returns:
            List of job information dictionaries
        """
        jobs = []

        try:
            result = subprocess.run(
                ['squeue', '-u', os.getenv('USER', 'unknown'), '-h', '-o', '%i,%j,%T,%M,%L,%P'],
                capture_output=True,
                text=True,
                timeout=15
            )

            if result.returncode == 0:
                for line in result.stdout.strip().split('\n'):
                    if line:
                        parts = line.split(',')
                        if len(parts) >= 6:
                            jobs.append({
                                'job_id': parts[0],
                                'name': parts[1],
                                'state': parts[2],
                                'time_used': parts[3],
                                'time_left': parts[4],
                                'partition': parts[5]
                            })

        except (subprocess.TimeoutExpired, FileNotFoundError):
            pass

        return jobs

    def get_job_output(self, job_id: str, output_dir: str) -> Tuple[Optional[str], Optional[str]]:
        """
        Get the stdout and stderr output files for a job.

        Args:
            job_id: SLURM job ID
            output_dir: Output directory containing log files

        Returns:
            Tuple of (stdout_content, stderr_content) or (None, None) if not found
        """
        logs_dir = Path(output_dir) / "logs"

        # Try to find log files matching the job pattern
        stdout_files = list(logs_dir.glob(f"*{job_id}*.out"))
        stderr_files = list(logs_dir.glob(f"*{job_id}*.err"))

        stdout_content = None
        stderr_content = None

        try:
            if stdout_files:
                with open(stdout_files[0], 'r') as f:
                    stdout_content = f.read()

            if stderr_files:
                with open(stderr_files[0], 'r') as f:
                    stderr_content = f.read()

        except FileNotFoundError:
            pass

        return stdout_content, stderr_content

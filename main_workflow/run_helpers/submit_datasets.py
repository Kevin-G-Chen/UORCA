#!/usr/bin/env python
import argparse, os, pathlib, subprocess, pandas as pd, time, json
from jinja2 import Environment, FileSystemLoader
from datetime import datetime

def display_scheduling_order(df_sorted, max_storage_gb):
    """Display the dataset scheduling order and storage analysis."""
    print(f"\n{'='*60}")
    print(f"DATASET SCHEDULING ORDER (Largest First)")
    print(f"{'='*60}")
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

def get_running_jobs_storage(output_dir):
    """Calculate total storage being used by currently running jobs."""
    total_storage = 0
    status_dir = os.path.join(output_dir, "job_status")

    if not os.path.exists(status_dir):
        return total_storage

    for status_file in os.listdir(status_dir):
        if status_file.endswith("_status.json"):
            status_path = os.path.join(status_dir, status_file)
            try:
                with open(status_path, 'r') as f:
                    status = json.load(f)
                if status.get('state') == 'running':
                    total_storage += status.get('storage_bytes', 0)
            except (json.JSONDecodeError, FileNotFoundError):
                continue

    return total_storage

def submit_single_dataset(accession, dataset_size_bytes, output_dir, resource_dir, cleanup, project_root):
    """Submit a single dataset for analysis."""
    logs_dir = os.path.join(output_dir, "logs")
    os.makedirs(logs_dir, exist_ok=True)

    status_dir = os.path.join(output_dir, "job_status")
    os.makedirs(status_dir, exist_ok=True)

    env = Environment(loader=FileSystemLoader("main_workflow/run_helpers"))

    script_txt = env.get_template("run_single_dataset.sbatch.j2").render(
        accession=accession,
        output_dir=output_dir,
        resource_dir=resource_dir,
        cleanup=cleanup,
        project_root=project_root,
        logs_dir=logs_dir,
    )

    script_path = pathlib.Path(f"run_{accession}.sbatch")
    script_path.write_text(script_txt)
    script_path.chmod(0o755)

    # Submit job and capture job ID
    result = subprocess.run(["sbatch", str(script_path)], capture_output=True, text=True, check=True)
    job_id = result.stdout.strip().split()[-1]  # Extract job ID from "Submitted batch job XXXXXX"

    # Clean up script file immediately after submission
    try:
        script_path.unlink()
    except OSError as e:
        print(f"Warning: Could not remove script file {script_path}: {e}")

    # Track job status (don't store script_path since we removed it)
    status_info = {
        'job_id': job_id,
        'accession': accession,
        'storage_bytes': dataset_size_bytes,
        'storage_gb': dataset_size_bytes / (1024**3),
        'state': 'running',
        'submitted_time': datetime.now().isoformat()
    }

    status_file = os.path.join(status_dir, f"{accession}_status.json")
    with open(status_file, 'w') as f:
        json.dump(status_info, f, indent=2)

    print(f"Submitted {accession} (Job ID: {job_id}, Size: {dataset_size_bytes/(1024**3):.2f} GB)")
    return job_id

def _checkpoint_summary(checkpoints):
    """Helper function to summarize checkpoint status without early breaks."""
    completed = failed = 0
    furthest = None
    for name, cp in checkpoints.items():
        status = cp.get("status", "not_started")
        err_msg = cp.get("error_message")
        if status == "completed" and not err_msg:
            completed += 1
            furthest = name
        elif status == "failed" or err_msg:
            failed += 1
    return completed, failed, furthest

def extract_analysis_results(accession, output_dir):
    """Extract analysis results from analysis_info.json using checkpoint system."""
    results = {
        'success': False,
        'reflection_iterations': 0,
        'checkpoints_completed': 0,
        'checkpoints_failed': 0,
        'furthest_checkpoint': 'none'
    }

    # Check analysis_info.json for success status - prioritize this as the authoritative source
    analysis_info_paths = [
        os.path.join(output_dir, accession, "metadata", "analysis_info.json"),
        os.path.join(output_dir, accession, "analysis_info.json")
    ]

    print(f"DEBUG [{accession}]: Checking analysis_info.json paths:")
    for path in analysis_info_paths:
        print(f"DEBUG [{accession}]:   - {path} (exists: {os.path.exists(path)})")

    analysis_info_found = False
    for info_path in analysis_info_paths:
        print(f"DEBUG [{accession}]: Checking path: {info_path}")
        if os.path.exists(info_path):
            print(f"DEBUG [{accession}]: File exists, attempting to read...")
            try:
                with open(info_path, 'r') as f:
                    analysis_info = json.load(f)

                print(f"DEBUG [{accession}]: Successfully loaded JSON with keys: {list(analysis_info.keys())}")
                print(
                    f"DEBUG [{accession}]: [{info_path}] 'analysis_success' -> {analysis_info.get('analysis_success', 'NOT_FOUND')}"
                )
                print(
                    f"DEBUG [{accession}]: [{info_path}] 'reflection_iterations' -> {analysis_info.get('reflection_iterations', 'NOT_FOUND')}"
                )

                # Trust analysis_success if present - this is the authoritative source
                if 'analysis_success' in analysis_info:
                    results['success'] = bool(analysis_info['analysis_success'])
                    print(
                        f"DEBUG [{accession}]: Using authoritative analysis_success from {info_path}: {results['success']}"
                    )

                analysis_info_found = True

                # Get reflection iterations from analysis_info.json
                results['reflection_iterations'] = analysis_info.get('reflection_iterations', 0)
                
                # Evaluate checkpoints if available (but don't override analysis_success unless missing)
                checkpoints = analysis_info.get('checkpoints', {})
                if checkpoints:
                    print(f"DEBUG [{accession}]: [{info_path}] Found checkpoints: {list(checkpoints.keys())}")

                    completed, failed, furthest = _checkpoint_summary(checkpoints)
                    results['checkpoints_completed'] = completed
                    results['checkpoints_failed'] = failed
                    results['furthest_checkpoint'] = furthest or 'none'

                    print(
                        f"DEBUG [{accession}]: [{info_path}] Checkpoint summary -> {completed} completed, {failed} failed, furthest: {furthest}"
                    )

                    # Only use checkpoint-based success if analysis_success was missing
                    if 'analysis_success' not in analysis_info:
                        results['success'] = failed == 0 and completed > 0
                        print(
                            f"DEBUG [{accession}]: Derived success from checkpoints in {info_path}: {results['success']}"
                        )
                    else:
                        # Check for contradiction: if analysis_success=True but we have failures, flag it
                        if results['success'] and failed > 0:
                            print(f"DEBUG [{accession}]: WARNING: analysis_success=True but {failed} checkpoint failures detected")
                
                print(
                    f"DEBUG [{accession}]: Final results from {info_path}: success={results['success']}, "
                    f"reflections={results['reflection_iterations']}, "
                    f"checkpoints={results['checkpoints_completed']}/{completed + failed + len([cp for cp in checkpoints.values() if cp.get('status') == 'not_started'])}"
                )
                break
            except (json.JSONDecodeError, FileNotFoundError, KeyError) as e:
                print(f"DEBUG [{accession}]: Error reading JSON: {str(e)}")
                continue
        else:
            print(f"DEBUG [{accession}]: File does not exist at {info_path}")
    
    if not analysis_info_found:
        print(f"DEBUG [{accession}]: No analysis_info.json found in any location")

    # Check log files for reflection iterations and additional error info
    log_dir = os.path.join(output_dir, "logs")
    log_file = os.path.join(log_dir, f"run_{accession}.out")

    if os.path.exists(log_file):
        try:
            print(f"DEBUG [{accession}]: Checking log file {log_file} for success indicators")
            with open(log_file, 'r') as f:
                log_content = f.read()

            # Count reflections from logs if analysis_info.json wasn't found, or take the maximum
            if not analysis_info_found:
                reflection_count = log_content.count("❌ Analysis attempt")
                results['reflection_iterations'] = max(0, reflection_count)
            else:
                # Take the maximum of JSON and log counts to handle manual restarts
                log_reflection_count = log_content.count("❌ Analysis attempt")
                results['reflection_iterations'] = max(results['reflection_iterations'], log_reflection_count)

            # If no analysis_info.json found, fall back to log-based success detection
            if not analysis_info_found:
                # Look for success indicators in logs
                if "✅ Analysis attempt" in log_content and "succeeded!" in log_content:
                    results['success'] = True
                    print(f"DEBUG [{accession}]: {log_file} -> success pattern found")
                elif "❌ Final analysis attempt failed" in log_content:
                    results['success'] = False
                    print(f"DEBUG [{accession}]: {log_file} -> failure pattern found")
                
        except Exception as e:
            print(f"DEBUG [{accession}]: Could not read log file: {str(e)}")

    # If we still don't have success info and no analysis_info.json was found, 
    # the analysis likely didn't complete properly
    if not analysis_info_found:
        print(f"DEBUG [{accession}]: No analysis_info.json found - analysis incomplete")

    return results

def update_job_statuses(output_dir):
    """Update status of all tracked jobs with robust error detection."""
    status_dir = os.path.join(output_dir, "job_status")

    if not os.path.exists(status_dir):
        return

    for status_file in os.listdir(status_dir):
        if status_file.endswith("_status.json"):
            status_path = os.path.join(status_dir, status_file)
            try:
                with open(status_path, 'r') as f:
                    status = json.load(f)

                if status.get('state') == 'running':
                    job_id = status.get('job_id')
                    accession = status.get('accession')

                    # First check job status using squeue with more detailed info
                    squeue_result = subprocess.run(
                        ["squeue", "-j", job_id, "-h", "-o", "%T,%R"],
                        capture_output=True, text=True
                    )

                    job_still_running = False
                    job_failed = False
                    failure_reason = None

                    if squeue_result.returncode == 0 and squeue_result.stdout.strip():
                        # Job is still in the queue
                        job_info = squeue_result.stdout.strip().split(',')
                        job_state = job_info[0] if job_info else "UNKNOWN"

                        if job_state in ["RUNNING", "PENDING", "CONFIGURING"]:
                            job_still_running = True
                        elif job_state in ["FAILED", "CANCELLED", "TIMEOUT", "NODE_FAIL"]:
                            job_failed = True
                            failure_reason = f"SLURM state: {job_state}"
                            if len(job_info) > 1:
                                failure_reason += f", Reason: {job_info[1]}"
                    else:
                        # Job not in queue - check if it completed or failed
                        # Use sacct to get completion status if available
                        sacct_result = subprocess.run(
                            ["sacct", "-j", job_id, "-n", "-o", "State,ExitCode", "--parsable2"],
                            capture_output=True, text=True
                        )

                        if sacct_result.returncode == 0 and sacct_result.stdout.strip():
                            sacct_lines = sacct_result.stdout.strip().split('\n')
                            for line in sacct_lines:
                                if line and not line.endswith('.batch') and not line.endswith('.extern'):
                                    parts = line.split('|')
                                    if len(parts) >= 2:
                                        state = parts[0]
                                        exit_code = parts[1]

                                        if state in ["FAILED", "CANCELLED", "TIMEOUT", "NODE_FAIL"] or exit_code != "0:0":
                                            job_failed = True
                                            failure_reason = f"SLURM state: {state}, Exit code: {exit_code}"
                                        break

                    # Check SLURM output files for immediate errors
                    if not job_still_running:
                        logs_dir = os.path.join(output_dir, "logs")
                        error_file = os.path.join(logs_dir, f"run_{accession}.err")
                        output_file = os.path.join(logs_dir, f"run_{accession}.out")

                        # Check error file for obvious failures
                        error_content = ""
                        if os.path.exists(error_file):
                            try:
                                with open(error_file, 'r') as f:
                                    error_content = f.read()

                                # Look for common error patterns
                                error_patterns = [
                                    "IndentationError", "SyntaxError", "ImportError",
                                    "ModuleNotFoundError", "FileNotFoundError",
                                    "command not found", "No such file or directory",
                                    "Permission denied", "Traceback (most recent call last)"
                                ]

                                for pattern in error_patterns:
                                    if pattern in error_content:
                                        job_failed = True
                                        if not failure_reason:
                                            failure_reason = f"Error detected in SLURM stderr: {pattern}"
                                        break
                            except Exception:
                                pass

                        # Also check if the job produced any output at all
                        if not job_failed and os.path.exists(output_file):
                            try:
                                with open(output_file, 'r') as f:
                                    output_content = f.read()

                                # If there's very little output and errors, likely failed early
                                if len(output_content.strip()) < 100 and error_content:
                                    job_failed = True
                                    if not failure_reason:
                                        failure_reason = "Job produced minimal output with errors"
                            except Exception:
                                pass

                    if job_still_running:
                        # Job is still running, no update needed
                        continue
                    elif job_failed:
                        # Job failed
                        status['state'] = 'failed'
                        status['completed_time'] = datetime.now().isoformat()
                        status['failure_reason'] = failure_reason or "Unknown failure"
                        status['success'] = False

                        # Still try to extract what results we can
                        analysis_results = extract_analysis_results(accession, output_dir)
                        status.update(analysis_results)

                        with open(status_path, 'w') as f:
                            json.dump(status, f, indent=2)

                        print(f"Job {job_id} ({accession}) FAILED ❌ - {failure_reason}")

                    else:
                        # Job completed (hopefully successfully)
                        status['state'] = 'completed'
                        status['completed_time'] = datetime.now().isoformat()

                        # Extract analysis results
                        analysis_results = extract_analysis_results(accession, output_dir)
                        status.update(analysis_results)

                        with open(status_path, 'w') as f:
                            json.dump(status, f, indent=2)

                        success_indicator = "✅" if status.get('success', False) else "❌"
                        reflection_info = f" ({status.get('reflection_iterations', 0)} reflections)" if status.get('reflection_iterations', 0) > 0 else ""
                        print(f"Job {job_id} ({accession}) completed {success_indicator}{reflection_info}")

                    # Script files are now cleaned up immediately after submission
                    pass

            except (json.JSONDecodeError, FileNotFoundError, subprocess.SubprocessError) as e:
                print(f"Warning: Error updating status for {status_file}: {str(e)}")
                continue

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv_file", required=True)
    ap.add_argument("--output_dir", default="../UORCA_results")
    ap.add_argument("--resource_dir", default="./data/kallisto_indices/")
    ap.add_argument("--max_parallel", type=int, default=10)
    ap.add_argument("--max_storage_gb", type=float, default=None,
                   help="Maximum total storage (in GB) that can be used by running datasets simultaneously")
    ap.add_argument("--cleanup", action="store_true",
                   help="Clean up FASTQ and SRA files after successful completion of each dataset")
    ap.add_argument("--check_interval", type=int, default=300,
                   help="Interval in seconds to check job status and submit new jobs")
    args = ap.parse_args()

    df = pd.read_csv(args.csv_file)

    # Validate CSV structure
    required_columns = ['Accession']
    if not all(col in df.columns for col in required_columns):
        print(f"Error: CSV must contain columns: {required_columns}")
        print(f"Found columns: {list(df.columns)}")
        exit(1)

    # Check for dataset size columns
    if 'DatasetSizeBytes' not in df.columns:
        if args.max_storage_gb is not None:
            print("Warning: DatasetSizeBytes column not found in CSV. Storage limits will be ignored.")
            args.max_storage_gb = None
        df['DatasetSizeBytes'] = 0

    # Double the dataset sizes for safety as requested
    df['SafeStorageBytes'] = df['DatasetSizeBytes'] * 2
    df['SafeStorageGB'] = df['SafeStorageBytes'] / (1024**3)

    if args.max_storage_gb is None:
        print("No storage limit specified. Using traditional array job submission...")
        # Fall back to original array job submission
        logs_dir = os.path.join(args.output_dir, "logs")
        os.makedirs(logs_dir, exist_ok=True)

        env = Environment(loader=FileSystemLoader("main_workflow/run_helpers"))

        script_txt = env.get_template("run_dataset_array.sbatch.j2").render(
            array_size=len(df) - 1,
            max_parallel=args.max_parallel,
            logs_dir=logs_dir,
            csv_file=os.path.abspath(args.csv_file),
            project_root=os.getcwd(),
            base_output_dir=args.output_dir,
            resource_dir=args.resource_dir,
            cleanup=args.cleanup,
        )

        script_path = pathlib.Path("run_datasets_array.sbatch")
        script_path.write_text(script_txt)
        script_path.chmod(0o755)

        subprocess.run(["sbatch", str(script_path)], check=True)
        
        # Clean up script file immediately after submission
        try:
            script_path.unlink()
        except OSError as e:
            print(f"Warning: Could not remove array script file {script_path}: {e}")
        
        return

    print(f"Storage-aware submission enabled. Maximum storage: {args.max_storage_gb:.2f} GB")
    print(f"Will submit {len(df)} datasets with storage management")

    max_storage_bytes = args.max_storage_gb * (1024**3)
    project_root = os.getcwd()

    # Sort datasets by size (largest first) to optimize completion time
    df_sorted = df.sort_values('SafeStorageBytes', ascending=False).reset_index(drop=True)

    # Validate datasets against storage limits
    oversized_datasets = df_sorted[df_sorted['SafeStorageGB'] > args.max_storage_gb]
    if len(oversized_datasets) > 0:
        print(f"\n⚠️  WARNING: {len(oversized_datasets)} dataset(s) exceed storage limit ({args.max_storage_gb:.2f} GB):")
        for _, row in oversized_datasets.iterrows():
            print(f"   {row['Accession']:12s} - {row['SafeStorageGB']:6.2f} GB")

        print(f"\nThese datasets will be skipped. Consider:")
        print(f"  1. Increasing --max_storage_gb to accommodate larger datasets")
        print(f"  2. Processing these datasets separately with higher storage limits")
        print(f"  3. Using the array job submission (remove --max_storage_gb) for unlimited storage")

        # Remove oversized datasets from processing
        df_sorted = df_sorted[df_sorted['SafeStorageGB'] <= args.max_storage_gb].reset_index(drop=True)

        if len(df_sorted) == 0:
            print(f"\n❌ ERROR: No datasets can fit within the storage limit of {args.max_storage_gb:.2f} GB")
            print(f"Please increase --max_storage_gb or remove the storage limit.")
            exit(1)

        print(f"\nContinuing with {len(df_sorted)} datasets that fit within storage limits...\n")

    # Display scheduling order and analysis
    display_scheduling_order(df_sorted, args.max_storage_gb)

    pending_datasets = df_sorted.to_dict('records')

    while pending_datasets:
        # Update job statuses
        update_job_statuses(args.output_dir)

        # Calculate current storage usage and running job count
        current_storage = get_running_jobs_storage(args.output_dir)
        status_dir = os.path.join(args.output_dir, "job_status")
        running_jobs_count = 0
        if os.path.exists(status_dir):
            for status_file in os.listdir(status_dir):
                if status_file.endswith("_status.json"):
                    status_path = os.path.join(status_dir, status_file)
                    try:
                        with open(status_path, 'r') as f:
                            status = json.load(f)
                        if status.get('state') == 'running':
                            running_jobs_count += 1
                    except (json.JSONDecodeError, FileNotFoundError):
                        continue

        # Check if we can submit more jobs
        if running_jobs_count >= args.max_parallel:
            print(f"At maximum parallel jobs ({args.max_parallel}), waiting for completions...")
            time.sleep(args.check_interval)
            continue

        # Try to submit datasets starting with largest that fits
        datasets_submitted_this_round = []
        datasets_to_remove = []

        for i, dataset in enumerate(pending_datasets):
            dataset_storage = dataset['SafeStorageBytes']

            # Check if this dataset fits within storage and job constraints
            if (current_storage + dataset_storage <= max_storage_bytes and
                running_jobs_count < args.max_parallel):

                # Submit this dataset
                job_id = submit_single_dataset(
                    dataset['Accession'],
                    dataset_storage,
                    args.output_dir,
                    args.resource_dir,
                    args.cleanup,
                    project_root
                )

                current_storage += dataset_storage
                running_jobs_count += 1
                datasets_to_remove.append(i)
                datasets_submitted_this_round.append(dataset['Accession'])

                print(f"Storage usage: {current_storage/(1024**3):.2f}/{args.max_storage_gb:.2f} GB")

                # Stop if we've reached max parallel jobs
                if running_jobs_count >= args.max_parallel:
                    print(f"Reached maximum parallel jobs ({args.max_parallel})")
                    break

        # Remove submitted datasets from pending list (reverse order to maintain indices)
        for i in reversed(datasets_to_remove):
            pending_datasets.pop(i)

        # If no datasets were submitted this round, wait for jobs to complete
        if not datasets_submitted_this_round:
            if pending_datasets:
                largest_pending = pending_datasets[0]  # First item is largest due to sorting
                print(f"⏳ No datasets fit current constraints. Largest pending: {largest_pending['Accession']} "
                      f"({largest_pending['SafeStorageGB']:.2f} GB). Waiting for jobs to complete...")
                print(f"📊 Current usage: {current_storage/(1024**3):.2f}/{args.max_storage_gb:.2f} GB, "
                      f"Running jobs: {running_jobs_count}/{args.max_parallel}")
            time.sleep(args.check_interval)
        else:
            print(f"✅ Submitted {len(datasets_submitted_this_round)} datasets: {', '.join(datasets_submitted_this_round)}")
            if pending_datasets:
                print(f"📋 {len(pending_datasets)} datasets remaining in queue")

    print(f"\nSubmission complete! Submitted {len(df)} datasets total.")
    print(f"Monitor progress with: watch 'ls {args.output_dir}/job_status/*.json | wc -l'")

    # Wait for all jobs to complete and generate final summary
    print("\nWaiting for all jobs to complete...")
    while True:
        status_dir = os.path.join(args.output_dir, "job_status")
        if not os.path.exists(status_dir):
            break

        running_jobs = []
        for status_file in os.listdir(status_dir):
            if status_file.endswith("_status.json"):
                status_path = os.path.join(status_dir, status_file)
                try:
                    with open(status_path, 'r') as f:
                        status = json.load(f)
                    if status.get('state') == 'running':
                        running_jobs.append(status.get('accession', 'unknown'))
                except:
                    continue

        if not running_jobs:
            break

        print(f"Still running: {len(running_jobs)} jobs ({', '.join(running_jobs[:5])}{'...' if len(running_jobs) > 5 else ''})")
        time.sleep(args.check_interval)
        update_job_statuses(args.output_dir)

    # Generate final summary report
    generate_summary_report(args.output_dir)

def generate_summary_report(output_dir):
    """Generate a summary report of all dataset analyses."""
    status_dir = os.path.join(output_dir, "job_status")

    if not os.path.exists(status_dir):
        print("No job status directory found.")
        return

    successful = []
    failed = []
    running = []

    for status_file in os.listdir(status_dir):
        if status_file.endswith("_status.json"):
            status_path = os.path.join(status_dir, status_file)
            try:
                with open(status_path, 'r') as f:
                    status = json.load(f)

                accession = status.get('accession', 'unknown')
                reflections = status.get('reflection_iterations', 0)
                job_state = status.get('state', 'unknown')

                if job_state == 'running':
                    running.append({
                        'accession': accession,
                        'storage_gb': status.get('storage_gb', 0),
                        'submitted_time': status.get('submitted_time', 'Unknown')
                    })
                elif status.get('success', False):
                    successful.append({
                        'accession': accession,
                        'reflections': reflections,
                        'storage_gb': status.get('storage_gb', 0),
                        'checkpoints_completed': status.get('checkpoints_completed', 0),
                        'furthest_checkpoint': status.get('furthest_checkpoint', 'unknown')
                    })
                else:
                    # Job failed or completed unsuccessfully
                    error_info = status.get('error_message', 'Unknown error')
                    failure_reason = status.get('failure_reason', None)

                    # Combine error information
                    if failure_reason and failure_reason not in error_info:
                        error_info = f"{failure_reason}; {error_info}" if error_info != 'Unknown error' else failure_reason

                    failed.append({
                        'accession': accession,
                        'reflections': reflections,
                        'storage_gb': status.get('storage_gb', 0),
                        'state': job_state,
                        'failure_reason': failure_reason,
                        'checkpoints_completed': status.get('checkpoints_completed', 0),
                        'furthest_checkpoint': status.get('furthest_checkpoint', 'none')
                    })
            except Exception as e:
                print(f"Warning: Could not read status file {status_file}: {str(e)}")
                continue

    # Write summary to file
    summary_path = os.path.join(output_dir, "batch_analysis_summary.txt")
    with open(summary_path, 'w') as f:
        f.write(f"UORCA Batch Analysis Summary\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"=" * 50 + "\n\n")

        total_jobs = len(successful) + len(failed) + len(running)
        completed_jobs = len(successful) + len(failed)

        f.write(f"Overall Results:\n")
        f.write(f"  Total datasets: {total_jobs}\n")
        f.write(f"  Completed: {completed_jobs}\n")
        f.write(f"  Still running: {len(running)}\n")
        f.write(f"  Successful: {len(successful)} ({len(successful)/completed_jobs*100:.1f}% of completed)\n")
        f.write(f"  Failed: {len(failed)} ({len(failed)/completed_jobs*100:.1f}% of completed)\n")
        f.write("\n")

        if successful:
            f.write(f"Successful Analyses:\n")
            total_size_successful = sum(s['storage_gb'] for s in successful)
            f.write(f"  Total size processed: {total_size_successful:.2f} GB\n")
            for s in sorted(successful, key=lambda x: x['reflections']):
                checkpoints_info = f"({s.get('checkpoints_completed', 0)} checkpoints, reached {s.get('furthest_checkpoint', 'unknown')})"
                f.write(f"  {s['accession']}: {s['reflections']} reflections, {s['storage_gb']:.2f} GB {checkpoints_info}\n")
            f.write("\n")

        if running:
            f.write(f"Still Running:\n")
            for s in running:
                f.write(f"  {s['accession']}: {s['storage_gb']:.2f} GB, submitted at {s['submitted_time']}\n")
            f.write("\n")

        if failed:
            f.write(f"Failed Analyses:\n")
            # Group failures by type
            slurm_failures = [s for s in failed if s.get('failure_reason') and 'SLURM' in s.get('failure_reason', '')]
            analysis_failures = [s for s in failed if s not in slurm_failures]

            if slurm_failures:
                f.write(f"  SLURM/System Failures ({len(slurm_failures)}):\n")
                for s in slurm_failures:
                    f.write(f"    {s['accession']}: {s['storage_gb']:.2f} GB\n")
                    f.write(f"      Reason: {s.get('failure_reason', 'Unknown')}\n")
                f.write("\n")

            if analysis_failures:
                f.write(f"  Analysis Failures ({len(analysis_failures)}):\n")
                for s in sorted(analysis_failures, key=lambda x: x['reflections'], reverse=True):
                    f.write(f"    {s['accession']}: {s['reflections']} reflections, {s['storage_gb']:.2f} GB\n")
                    f.write(f"      Furthest checkpoint: {s.get('furthest_checkpoint', 'unknown')} ({s.get('checkpoints_completed', 0)}/7 completed)\n")

    print(f"\n" + "=" * 50)
    print(f"BATCH ANALYSIS SUMMARY")
    print(f"=" * 50)
    total_jobs = len(successful) + len(failed) + len(running)
    completed_jobs = len(successful) + len(failed)
    print(f"Total datasets: {total_jobs}")
    print(f"Completed: {completed_jobs}")
    print(f"Still running: {len(running)}")
    if completed_jobs > 0:
        print(f"Successful: {len(successful)} ({len(successful)/completed_jobs*100:.1f}% of completed)")
        print(f"Failed: {len(failed)} ({len(failed)/completed_jobs*100:.1f}% of completed)")
    else:
        print("No jobs have completed yet.")
    # Write detailed CSV for easy processing
    csv_path = os.path.join(output_dir, "batch_analysis_results.csv")
    import csv

    with open(csv_path, 'w', newline='') as csvfile:
        fieldnames = ['accession', 'success', 'reflection_iterations', 'storage_gb', 'checkpoints_completed', 'furthest_checkpoint']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        # Write successful datasets
        for s in successful:
            writer.writerow({
                'accession': s['accession'],
                'success': True,
                'reflection_iterations': s['reflections'],
                'storage_gb': s['storage_gb'],
                'checkpoints_completed': s.get('checkpoints_completed', 0),
                'furthest_checkpoint': s.get('furthest_checkpoint', 'unknown')
            })

        # Write still running datasets
        for s in running:
            writer.writerow({
                'accession': s['accession'],
                'success': None,
                'reflection_iterations': 0,
                'storage_gb': s['storage_gb'],
                'checkpoints_completed': 0,
                'furthest_checkpoint': 'none'
            })

        # Write failed datasets
        for s in failed:
            writer.writerow({
                'accession': s['accession'],
                'success': False,
                'reflection_iterations': s['reflections'],
                'storage_gb': s['storage_gb'],
                'checkpoints_completed': s.get('checkpoints_completed', 0),
                'furthest_checkpoint': s.get('furthest_checkpoint', 'none')
            })

    print(f"\nDetailed summary saved to: {summary_path}")
    print(f"CSV results saved to: {csv_path}")

    if running:
        print(f"\nStill running datasets:")
        for s in running[:5]:  # Show first 5 running
            print(f"  {s['accession']}: {s['storage_gb']:.2f} GB")
        if len(running) > 5:
            print(f"  ... and {len(running)-5} more")

    if failed:
        print(f"\nFailed datasets:")
        slurm_failures = [s for s in failed if s.get('failure_reason') and 'SLURM' in s.get('failure_reason', '')]
        analysis_failures = [s for s in failed if s not in slurm_failures]

        if slurm_failures:
            print(f"  SLURM/System failures: {len(slurm_failures)}")
            for s in slurm_failures[:3]:
                print(f"    {s['accession']}: {s.get('failure_reason', 'Unknown')}")

        if analysis_failures:
            print(f"  Analysis failures: {len(analysis_failures)}")
            for s in analysis_failures[:3]:
                furthest = s.get('furthest_checkpoint', 'none')
                completed = s.get('checkpoints_completed', 0)
                print(f"    {s['accession']}: {s['reflections']} reflections, reached {furthest} ({completed}/7)")

        if len(failed) > 6:
            print(f"  ... and {len(failed)-6} more (see {summary_path} for details)")

def test_analysis_extraction():
    """Enhanced test function to demonstrate the fix for checkpoint analysis"""
    print("=== TESTING ANALYSIS EXTRACTION ===")
    
    # Test with the sample file
    test_output_dir = "UORCA_results/2025-06-03_UpdatedEvaluation"
    test_accession = "GSE193767"
    
    print(f"Testing with output_dir: {test_output_dir}")
    print(f"Testing with accession: {test_accession}")
    
    # First, let's show the raw JSON content
    analysis_info_paths = [
        os.path.join(test_output_dir, test_accession, "metadata", "analysis_info.json"),
        os.path.join(test_output_dir, test_accession, "analysis_info.json")
    ]
    
    print(f"\n=== RAW JSON ANALYSIS ===")
    for path in analysis_info_paths:
        if os.path.exists(path):
            print(f"Found analysis_info.json at: {path}")
            try:
                with open(path, 'r') as f:
                    analysis_info = json.load(f)
                
                print(f"Raw JSON keys: {list(analysis_info.keys())}")
                print(f"analysis_success: {analysis_info.get('analysis_success', 'MISSING')}")
                print(f"checkpoints present: {'checkpoints' in analysis_info}")
                
                if 'checkpoints' in analysis_info:
                    checkpoints = analysis_info['checkpoints']
                    print(f"Checkpoint keys: {list(checkpoints.keys())}")
                    for name, cp in checkpoints.items():
                        status = cp.get('status', 'unknown')
                        error = cp.get('error_message', None)
                        print(f"  {name}: {status}" + (f" (ERROR: {error})" if error else ""))
                    
                    # Demonstrate the checkpoint summary function
                    completed, failed, furthest = _checkpoint_summary(checkpoints)
                    print(f"Checkpoint summary: {completed} completed, {failed} failed, furthest: {furthest}")
                    
                    # Show old logic vs new logic
                    print(f"\n--- OLD LOGIC (broken) ---")
                    print("Would break on first non-completed checkpoint and override analysis_success")
                    
                    print(f"\n--- NEW LOGIC (fixed) ---")
                    if 'analysis_success' in analysis_info:
                        print(f"Trusts analysis_success={analysis_info['analysis_success']} as authoritative")
                        if failed > 0:
                            print(f"WARNING: Would flag contradiction - success=True but {failed} failures")
                    else:
                        derived_success = failed == 0 and completed > 0
                        print(f"No analysis_success field, would derive: {derived_success}")
                
                break
            except Exception as e:
                print(f"Error reading {path}: {e}")
        else:
            print(f"File not found: {path}")
    
    print(f"\n=== EXTRACTION RESULTS ===")
    results = extract_analysis_results(test_accession, test_output_dir)
    
    print(f"\nFINAL RESULTS:")
    print(f"  Success: {results['success']}")
    print(f"  Reflection iterations: {results['reflection_iterations']}")
    print(f"  Checkpoints completed: {results.get('checkpoints_completed', 'N/A')}")
    print(f"  Checkpoints failed: {results.get('checkpoints_failed', 'N/A')}")
    print(f"  Furthest checkpoint: {results.get('furthest_checkpoint', 'N/A')}")
    
    # Test multiple scenarios
    print(f"\n=== TESTING MULTIPLE SCENARIOS ===")
    
    # Scenario 1: Missing analysis_success but good checkpoints
    print("Scenario 1: Missing analysis_success, should derive from checkpoints")
    test_checkpoints = {
        'metadata_extraction': {'status': 'completed'},
        'fastq_extraction': {'status': 'completed'},
        'rnaseq_analysis': {'status': 'completed'}
    }
    comp, fail, furt = _checkpoint_summary(test_checkpoints)
    derived = fail == 0 and comp > 0
    print(f"  {comp} completed, {fail} failed → success: {derived}")
    
    # Scenario 2: Has analysis_success=True with some checkpoints skipped
    print("Scenario 2: analysis_success=True with skipped checkpoints (should trust True)")
    test_checkpoints = {
        'metadata_extraction': {'status': 'completed'},
        'fastq_extraction': {'status': 'not_started'},  # Skipped
        'rnaseq_analysis': {'status': 'completed'}
    }
    comp, fail, furt = _checkpoint_summary(test_checkpoints)
    print(f"  {comp} completed, {fail} failed, would trust analysis_success=True")
    
    # Scenario 3: Contradiction case
    print("Scenario 3: analysis_success=True but checkpoint failures (should flag warning)")
    test_checkpoints = {
        'metadata_extraction': {'status': 'completed'},
        'fastq_extraction': {'status': 'failed', 'error_message': 'Download failed'},
        'rnaseq_analysis': {'status': 'completed'}
    }
    comp, fail, furt = _checkpoint_summary(test_checkpoints)
    print(f"  {comp} completed, {fail} failed, would flag contradiction with analysis_success=True")
    
    print("=== END TEST ===")
    
    return results

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "test":
        test_analysis_extraction()
    else:
        main()

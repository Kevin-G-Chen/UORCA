#!/usr/bin/env python
import argparse, os, pathlib, subprocess, pandas as pd, time, json
from jinja2 import Environment, FileSystemLoader
from datetime import datetime

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
    
    # Track job status
    status_info = {
        'job_id': job_id,
        'accession': accession,
        'storage_bytes': dataset_size_bytes,
        'storage_gb': dataset_size_bytes / (1024**3),
        'state': 'running',
        'submitted_time': datetime.now().isoformat(),
        'script_path': str(script_path)
    }
    
    status_file = os.path.join(status_dir, f"{accession}_status.json")
    with open(status_file, 'w') as f:
        json.dump(status_info, f, indent=2)
    
    print(f"Submitted {accession} (Job ID: {job_id}, Size: {dataset_size_bytes/(1024**3):.2f} GB)")
    return job_id

def extract_analysis_results(accession, output_dir):
    """Extract analysis results from log files and analysis_info.json."""
    results = {
        'success': False,
        'reflection_iterations': 0,
        'error_message': None
    }
    
    # Check analysis_info.json for success status
    analysis_info_paths = [
        os.path.join(output_dir, accession, "metadata", "analysis_info.json"),
        os.path.join(output_dir, accession, "analysis_info.json")
    ]
    
    for info_path in analysis_info_paths:
        if os.path.exists(info_path):
            try:
                with open(info_path, 'r') as f:
                    analysis_info = json.load(f)
                results['success'] = analysis_info.get('analysis_success', False)
                break
            except (json.JSONDecodeError, FileNotFoundError):
                continue
    
    # Check log files for reflection iterations
    log_dir = os.path.join(output_dir, "logs")
    log_file = os.path.join(log_dir, f"run_{accession}.out")
    
    if os.path.exists(log_file):
        try:
            with open(log_file, 'r') as f:
                log_content = f.read()
            
            # Count reflection iterations
            reflection_count = log_content.count("❌ Analysis attempt")
            results['reflection_iterations'] = max(0, reflection_count)
            
            # Extract error messages if analysis failed
            if not results['success']:
                error_lines = []
                for line in log_content.split('\n'):
                    if 'ERROR' in line or 'FAILED' in line or 'Exception' in line:
                        error_lines.append(line.strip())
                if error_lines:
                    results['error_message'] = '; '.join(error_lines[-3:])  # Last 3 error lines
        except Exception as e:
            results['error_message'] = f"Could not read log file: {str(e)}"
    
    return results

def update_job_statuses(output_dir):
    """Update status of all tracked jobs."""
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
                    
                    # Check job status using squeue
                    result = subprocess.run(
                        ["squeue", "-j", job_id, "-h", "-o", "%T"],
                        capture_output=True, text=True
                    )
                    
                    if result.returncode != 0:
                        # Job not found in queue, likely completed
                        status['state'] = 'completed'
                        status['completed_time'] = datetime.now().isoformat()
                        
                        # Extract analysis results
                        analysis_results = extract_analysis_results(accession, output_dir)
                        status.update(analysis_results)
                        
                        with open(status_path, 'w') as f:
                            json.dump(status, f, indent=2)
                        
                        success_indicator = "✅" if status['success'] else "❌"
                        reflection_info = f" ({status['reflection_iterations']} reflections)" if status['reflection_iterations'] > 0 else ""
                        print(f"Job {job_id} ({accession}) completed {success_indicator}{reflection_info}")
                        
                        # Clean up script file
                        script_path = status.get('script_path')
                        if script_path and os.path.exists(script_path):
                            os.remove(script_path)
            
            except (json.JSONDecodeError, FileNotFoundError, subprocess.SubprocessError):
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
        return
    
    print(f"Storage-aware submission enabled. Maximum storage: {args.max_storage_gb:.2f} GB")
    print(f"Will submit {len(df)} datasets with storage management")
    
    max_storage_bytes = args.max_storage_gb * (1024**3)
    project_root = os.getcwd()
    
    # Sort datasets by size (smallest first) to maximize utilization
    df_sorted = df.sort_values('SafeStorageBytes').reset_index(drop=True)
    
    submitted_count = 0
    pending_datasets = df_sorted.to_dict('records')
    
    print(f"\nDataset submission summary:")
    for dataset in pending_datasets:
        print(f"  {dataset['Accession']}: {dataset['SafeStorageGB']:.2f} GB (doubled for safety)")
    print()
    
    while pending_datasets or submitted_count > 0:
        # Update job statuses
        update_job_statuses(args.output_dir)
        
        # Calculate current storage usage
        current_storage = get_running_jobs_storage(args.output_dir)
        
        # Try to submit new datasets
        datasets_to_remove = []
        for i, dataset in enumerate(pending_datasets):
            dataset_storage = dataset['SafeStorageBytes']
            
            if current_storage + dataset_storage <= max_storage_bytes:
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
                submitted_count += 1
                datasets_to_remove.append(i)
                
                print(f"Storage usage: {current_storage/(1024**3):.2f}/{args.max_storage_gb:.2f} GB")
                
                # Check if we've hit max parallel jobs
                running_jobs = len([f for f in os.listdir(os.path.join(args.output_dir, "job_status", "")) 
                                  if f.endswith("_status.json")])
                if running_jobs >= args.max_parallel:
                    print(f"Reached maximum parallel jobs ({args.max_parallel})")
                    break
        
        # Remove submitted datasets from pending list (reverse order to maintain indices)
        for i in reversed(datasets_to_remove):
            pending_datasets.pop(i)
        
        if pending_datasets:
            print(f"Waiting for jobs to complete... ({len(pending_datasets)} datasets pending)")
            time.sleep(args.check_interval)
        elif submitted_count == 0:
            print("All datasets submitted!")
            break
    
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
    total_reflections = 0
    
    for status_file in os.listdir(status_dir):
        if status_file.endswith("_status.json"):
            status_path = os.path.join(status_dir, status_file)
            try:
                with open(status_path, 'r') as f:
                    status = json.load(f)
                
                accession = status.get('accession', 'unknown')
                reflections = status.get('reflection_iterations', 0)
                total_reflections += reflections
                
                if status.get('success', False):
                    successful.append({
                        'accession': accession,
                        'reflections': reflections,
                        'storage_gb': status.get('storage_gb', 0)
                    })
                else:
                    failed.append({
                        'accession': accession,
                        'reflections': reflections,
                        'error': status.get('error_message', 'Unknown error'),
                        'storage_gb': status.get('storage_gb', 0)
                    })
            except:
                continue
    
    # Write summary to file
    summary_path = os.path.join(output_dir, "batch_analysis_summary.txt")
    with open(summary_path, 'w') as f:
        f.write(f"UORCA Batch Analysis Summary\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"=" * 50 + "\n\n")
        
        f.write(f"Overall Results:\n")
        f.write(f"  Total datasets: {len(successful) + len(failed)}\n")
        f.write(f"  Successful: {len(successful)} ({len(successful)/(len(successful)+len(failed))*100:.1f}%)\n")
        f.write(f"  Failed: {len(failed)} ({len(failed)/(len(successful)+len(failed))*100:.1f}%)\n")
        f.write(f"  Total reflection iterations: {total_reflections}\n")
        f.write(f"  Average reflections per dataset: {total_reflections/(len(successful)+len(failed)):.1f}\n\n")
        
        if successful:
            f.write(f"Successful Analyses:\n")
            total_size_successful = sum(s['storage_gb'] for s in successful)
            f.write(f"  Total size processed: {total_size_successful:.2f} GB\n")
            for s in sorted(successful, key=lambda x: x['reflections']):
                f.write(f"  {s['accession']}: {s['reflections']} reflections, {s['storage_gb']:.2f} GB\n")
            f.write("\n")
        
        if failed:
            f.write(f"Failed Analyses:\n")
            for s in sorted(failed, key=lambda x: x['reflections'], reverse=True):
                f.write(f"  {s['accession']}: {s['reflections']} reflections, {s['storage_gb']:.2f} GB\n")
                f.write(f"    Error: {s['error'][:100]}{'...' if len(s['error']) > 100 else ''}\n")
    
    print(f"\n" + "=" * 50)
    print(f"BATCH ANALYSIS SUMMARY")
    print(f"=" * 50)
    print(f"Total datasets: {len(successful) + len(failed)}")
    print(f"Successful: {len(successful)} ({len(successful)/(len(successful)+len(failed))*100:.1f}%)")
    print(f"Failed: {len(failed)} ({len(failed)/(len(successful)+len(failed))*100:.1f}%)")
    print(f"Total reflection iterations: {total_reflections}")
    print(f"Average reflections per dataset: {total_reflections/(len(successful)+len(failed)):.1f}")
    # Write detailed CSV for easy processing
    csv_path = os.path.join(output_dir, "batch_analysis_results.csv")
    import csv
    
    with open(csv_path, 'w', newline='') as csvfile:
        fieldnames = ['accession', 'success', 'reflection_iterations', 'storage_gb', 'error_message']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        
        # Write successful datasets
        for s in successful:
            writer.writerow({
                'accession': s['accession'],
                'success': True,
                'reflection_iterations': s['reflections'],
                'storage_gb': s['storage_gb'],
                'error_message': ''
            })
        
        # Write failed datasets
        for s in failed:
            writer.writerow({
                'accession': s['accession'],
                'success': False,
                'reflection_iterations': s['reflections'],
                'storage_gb': s['storage_gb'],
                'error_message': s['error']
            })
    
    print(f"\nDetailed summary saved to: {summary_path}")
    print(f"CSV results saved to: {csv_path}")
    
    if failed:
        print(f"\nFailed datasets:")
        for s in failed[:5]:  # Show first 5 failures
            print(f"  {s['accession']}: {s['reflections']} reflections")
        if len(failed) > 5:
            print(f"  ... and {len(failed)-5} more (see {summary_path} for details)")

if __name__ == "__main__":
    main()

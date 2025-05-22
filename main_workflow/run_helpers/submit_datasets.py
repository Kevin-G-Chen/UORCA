#!/usr/bin/env python
import argparse, os, pathlib, subprocess, pandas as pd
from jinja2 import Environment, FileSystemLoader

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv_file", required=True)
    ap.add_argument("--output_dir", default="../UORCA_results")
    ap.add_argument("--resource_dir", default="./data/kallisto_indices/")
    ap.add_argument("--max_parallel", type=int, default=10)
    ap.add_argument("--cleanup", action="store_true",
                   help="Clean up FASTQ and SRA files after successful completion of each dataset")
    args = ap.parse_args()

    logs_dir = os.path.join(args.output_dir, "logs")
    os.makedirs(logs_dir, exist_ok=True)

    env = Environment(loader=FileSystemLoader("main_workflow/run_helpers"))

    df = pd.read_csv(args.csv_file)

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

if __name__ == "__main__":
    main()

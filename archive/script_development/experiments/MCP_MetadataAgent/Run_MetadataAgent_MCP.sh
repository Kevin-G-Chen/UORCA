#!/usr/bin/env bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --partition=tki_agpdev
#SBATCH --constraint="clx"
#SBATCH -t 0-1:00
#SBATCH -o ./sbatch_output/run_metadata_mcp_server.out
#SBATCH -e ./sbatch_output/run_metadata_mcp_server.out.err

source ~/.bashrc

# Run the MCP server
uv run client.py \
--metadata_path ../../../../data/BenchmarkDatasets/GSE123821/GSE123821_metadata.csv

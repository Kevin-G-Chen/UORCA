#!/usr/bin/env bash

# Script to run UORCA Explorer using Apptainer container
# Usage: ./run_uorca_explorer.sh [results_directory] [port]
# Default results directory: ../UORCA_results
# Default port: 8501

set -e

# Function to display usage
show_usage() {
    echo "Usage: $0 [results_directory] [port]"
    echo ""
    echo "Arguments:"
    echo "  results_directory    Path to UORCA results directory (default: ../UORCA_results)"
    echo "  port                Port number for the web app (default: 8501)"
    echo ""
    echo "Examples:"
    echo "  $0                              # Use defaults"
    echo "  $0 /path/to/results             # Custom results directory"
    echo "  $0 /path/to/results 8502        # Custom directory and port"
    echo ""
    echo "After starting, access via SSH tunnel:"
    echo "  ssh -L 8000:127.0.0.1:PORT your_username@hpc_server"
    echo "  Then open http://127.0.0.1:8000 in your browser"
}

# Parse arguments
if [[ "$1" == "-h" ]] || [[ "$1" == "--help" ]]; then
    show_usage
    exit 0
fi

RESULTS_DIR=${1:-"../UORCA_results"}
PORT=${2:-8501}

# Validate port number
if ! [[ "$PORT" =~ ^[0-9]+$ ]] || [ "$PORT" -lt 1024 ] || [ "$PORT" -gt 65535 ]; then
    echo "Error: Invalid port number '$PORT'. Please use a number between 1024-65535."
    exit 1
fi

echo "=========================================="
echo "UORCA Explorer Container Launcher"
echo "=========================================="
echo "Results directory: $RESULTS_DIR"
echo "Port: $PORT"
echo "=========================================="

# Check if running in SLURM environment
if [[ -n "$SLURM_JOB_ID" ]]; then
    echo "Warning: Running inside a SLURM job. This may not work as expected for interactive web apps."
    echo "Consider running this script outside of SLURM for better interactivity."
fi

# Load apptainer module
echo "Loading Apptainer module..."
if ! module load apptainer 2>/dev/null; then
    echo "Error: Failed to load Apptainer module. Please ensure it's available on your system."
    exit 1
fi

# Create temporary directory
TEMP_DIR=$(pwd)/tmp_apptainer_explorer
mkdir -p ${TEMP_DIR}

# Path to your .sif
SIF=/data/tki_agpdev/kevin/phd/aim1/UORCA/scratch/container_testing/uorca_0.1.0.sif

# Verify SIF exists
if [ ! -f "$SIF" ]; then
    echo "Error: Container file not found at $SIF"
    echo "Please ensure the container is built and accessible."
    echo "You may need to update the SIF path in this script to match your container location."
    exit 1
fi

# Verify results directory exists
if [ ! -d "$RESULTS_DIR" ]; then
    echo "Error: Results directory not found at $RESULTS_DIR"
    echo "Please provide a valid results directory path."
    echo ""
    show_usage
    exit 1
fi

# Check if results directory contains expected UORCA output structure
if [ ! -d "$RESULTS_DIR" ] || [ -z "$(find "$RESULTS_DIR" -name "*.csv" -o -name "analysis_info.json" 2>/dev/null)" ]; then
    echo "Warning: Results directory exists but may not contain UORCA analysis results."
    echo "Make sure you're pointing to a directory with completed UORCA analyses."
fi

# Convert to absolute path for bind mounting
RESULTS_DIR=$(realpath "$RESULTS_DIR")
PROJECT_ROOT=$(pwd)

echo "Using temporary directory: ${TEMP_DIR}"
echo "Container: $SIF"
echo "Project root: $PROJECT_ROOT"
echo "Results directory: $RESULTS_DIR"

# Bind-mounts:
#  - your repo → /workspace
#  - the host's results directory → /UORCA_results in the container
BIND="-B ${PROJECT_ROOT}:/workspace \
      -B ${RESULTS_DIR}:/UORCA_results \
      -B ${TEMP_DIR}:/tmp"

# Check if port is already in use
if ss -tuln 2>/dev/null | grep -q ":${PORT} "; then
    echo "Warning: Port $PORT appears to be in use. The app may fail to start."
    echo "Consider using a different port or stopping other services on this port."
fi

echo ""
echo "Starting UORCA Explorer..."
echo "Container will start momentarily. Please wait for startup messages..."
echo ""
echo "To access the app:"
echo "1. On your laptop, run: ssh -L 8000:127.0.0.1:${PORT} $(whoami)@$(hostname)"
echo "2. Open http://127.0.0.1:8000 in your browser"
echo ""
echo "Press Ctrl+C to stop the application"
echo "=========================================="

# Cleanup function
cleanup() {
    echo ""
    echo "Shutting down UORCA Explorer..."
    echo "Cleaning up temporary directory..."
    rm -rf ${TEMP_DIR} 2>/dev/null || true
    echo "UORCA Explorer stopped. Thank you for using UORCA!"
}

# Set trap to cleanup on exit and common signals
trap cleanup EXIT INT TERM

# Run the container with Streamlit
echo "Launching Streamlit app in container..."
echo "This may take a moment to initialize..."

if ! apptainer exec \
  $BIND \
  --tmpdir=${TEMP_DIR} \
  --cleanenv \
  $SIF \
  bash -lc "\
    cd /workspace && \
    echo 'Container started successfully. Initializing Streamlit...' && \
    uv run streamlit run main_workflow/reporting/uorca_explorer.py --server.port ${PORT} --server.address 0.0.0.0 --server.headless true
  "; then
    echo ""
    echo "Error: Failed to start UORCA Explorer."
    echo "Check the error messages above for troubleshooting information."
    exit 1
fi
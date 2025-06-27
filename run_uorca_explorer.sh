#!/usr/bin/env bash

# Script to run UORCA Explorer using Apptainer or Docker container
# Usage: ./run_uorca_explorer.sh [OPTIONS] [results_directory] [port]
# Default results directory: ../UORCA_results
# Default port: 8501

set -e

# Flag to prevent double cleanup
CLEANUP_DONE=0

# Default values
ENGINE="auto"
DOCKER_IMAGE="kevingchen/uorca:0.1.0"
SIF_PATH="/data/tki_agpdev/kevin/phd/aim1/UORCA/scratch/container_testing/uorca_0.1.0.sif"
ENV_FILE=""

# Helper functions for tool detection
have_apptainer() { command -v apptainer &>/dev/null; }
have_docker() { command -v docker &>/dev/null; }

# Function to display usage
show_usage() {
    echo "Usage: $0 [OPTIONS] [results_directory] [port]"
    echo ""
    echo "Options:"
    echo "  -e, --engine ENGINE     Container engine to use: auto|apptainer|docker (default: auto)"
    echo "  --image IMAGE          Docker image to use (default: kevingchen/uorca:0.1.0)"
    echo "  --sif PATH             Apptainer .sif file path (default: built-in path)"
    echo "  --env-file FILE        Path to .env file for Docker (default: auto-detect)"
    echo "  -h, --help             Show this help message"
    echo ""
    echo "Arguments:"
    echo "  results_directory      Path to UORCA results directory (default: ../UORCA_results)"
    echo "  port                   Port number for the web app (default: 8501)"
    echo ""
    echo "Engine selection:"
    echo "  auto                   Prefer Apptainer if available, else Docker"
    echo "  apptainer             Force use of Apptainer"
    echo "  docker                Force use of Docker"
    echo ""
    echo "Examples:"
    echo "  $0                                    # Use defaults (auto-detect engine)"
    echo "  $0 -e docker                         # Force Docker engine"
    echo "  $0 -e apptainer /path/to/results     # Force Apptainer with custom results"
    echo "  $0 /path/to/results 8502             # Custom directory and port"
    echo "  $0 --image myuser/uorca:latest       # Custom Docker image"
    echo "  $0 --env-file /path/to/.env          # Custom .env file location"
    echo ""
    echo "Environment Variables:"
    echo "  For Docker: OPENAI_API_KEY and other environment variables can be set via:"
    echo "  - .env file in project root (auto-detected)"
    echo "  - Custom .env file with --env-file option"
    echo "  - System environment variables (automatically passed)"
    echo ""
    echo "After starting, access via SSH tunnel (if on remote server):"
    echo "  ssh -L 8000:127.0.0.1:PORT your_username@hpc_server"
    echo "  Then open http://127.0.0.1:8000 in your browser"
    echo ""
    echo "For local usage, simply open http://127.0.0.1:PORT in your browser"
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)
            show_usage
            exit 0
            ;;
        -e|--engine)
            if [[ -z "$2" ]] || [[ "$2" =~ ^- ]]; then
                echo "Error: --engine requires a value (auto|apptainer|docker)"
                exit 1
            fi
            ENGINE="$2"
            shift 2
            ;;
        --image)
            if [[ -z "$2" ]] || [[ "$2" =~ ^- ]]; then
                echo "Error: --image requires a value"
                exit 1
            fi
            DOCKER_IMAGE="$2"
            shift 2
            ;;
        --sif)
            if [[ -z "$2" ]] || [[ "$2" =~ ^- ]]; then
                echo "Error: --sif requires a value"
                exit 1
            fi
            SIF_PATH="$2"
            shift 2
            ;;
        --env-file)
            if [[ -z "$2" ]] || [[ "$2" =~ ^- ]]; then
                echo "Error: --env-file requires a value"
                exit 1
            fi
            ENV_FILE="$2"
            shift 2
            ;;
        -*)
            echo "Error: Unknown option $1"
            show_usage
            exit 1
            ;;
        *)
            # First positional argument is results directory
            break
            ;;
    esac
done

# Parse positional arguments
RESULTS_DIR=${1:-"../UORCA_results"}
PORT=${2:-8501}

# Validate engine selection
if [[ "$ENGINE" != "auto" && "$ENGINE" != "apptainer" && "$ENGINE" != "docker" ]]; then
    echo "Error: Invalid engine '$ENGINE'. Must be one of: auto, apptainer, docker"
    exit 1
fi

# Validate port number
if ! [[ "$PORT" =~ ^[0-9]+$ ]] || [ "$PORT" -lt 1024 ] || [ "$PORT" -gt 65535 ]; then
    echo "Error: Invalid port number '$PORT'. Please use a number between 1024-65535."
    exit 1
fi

# Resolve container engine
if [[ "$ENGINE" == "auto" ]]; then
    if have_apptainer; then
        ENGINE="apptainer"
        echo "Auto-detected: Using Apptainer"
    elif have_docker; then
        ENGINE="docker"
        echo "Auto-detected: Using Docker"
    else
        echo "Error: Neither Apptainer nor Docker is available on this system."
        echo "Please install one of these container engines or specify a different engine."
        exit 1
    fi
else
    # User explicitly chose an engine - validate it exists
    if [[ "$ENGINE" == "apptainer" ]] && ! have_apptainer; then
        echo "Error: Apptainer not found but was explicitly requested."
        echo "Please install Apptainer or use --engine docker"
        exit 1
    fi
    if [[ "$ENGINE" == "docker" ]] && ! have_docker; then
        echo "Error: Docker not found but was explicitly requested."
        echo "Please install Docker or use --engine apptainer"
        exit 1
    fi
fi

echo "=========================================="
echo "UORCA Explorer Container Launcher"
echo "=========================================="
echo "Container engine: $ENGINE"
echo "Results directory: $RESULTS_DIR"
echo "Port: $PORT"
if [[ "$ENGINE" == "docker" ]]; then
    echo "Docker image: $DOCKER_IMAGE"
else
    echo "Apptainer SIF: $SIF_PATH"
fi
echo "=========================================="

# Check if running in SLURM environment
if [[ -n "$SLURM_JOB_ID" ]]; then
    echo "Warning: Running inside a SLURM job. This may not work as expected for interactive web apps."
    echo "Consider running this script outside of SLURM for better interactivity."
fi

# Engine-specific setup
if [[ "$ENGINE" == "apptainer" ]]; then
    # Load apptainer module
    echo "Loading Apptainer module..."
    if ! module load apptainer 2>/dev/null; then
        echo "Error: Failed to load Apptainer module. Please ensure it's available on your system."
        exit 1
    fi

    # Verify SIF exists
    if [ ! -f "$SIF_PATH" ]; then
        echo "Error: Container file not found at $SIF_PATH"
        echo "Please ensure the container is built and accessible."
        echo "You may need to update the SIF path using --sif option."
        exit 1
    fi
else
    # Docker setup
    echo "Checking Docker image availability..."
    if ! docker image inspect "$DOCKER_IMAGE" &>/dev/null; then
        echo "Docker image $DOCKER_IMAGE not found locally. Attempting to pull..."
        if ! docker pull "$DOCKER_IMAGE"; then
            echo "Error: Failed to pull Docker image $DOCKER_IMAGE"
            echo "Please check the image name and your internet connection."
            exit 1
        fi
    fi
fi

# Handle environment variables for Docker
if [[ "$ENGINE" == "docker" ]]; then
    # Auto-detect .env file if not specified
    if [[ -z "$ENV_FILE" ]]; then
        if [[ -f "$(pwd)/.env" ]]; then
            ENV_FILE="$(pwd)/.env"
            echo "Auto-detected .env file: $ENV_FILE"
        fi
    fi

    # Validate custom .env file if specified
    if [[ -n "$ENV_FILE" ]] && [[ ! -f "$ENV_FILE" ]]; then
        echo "Error: Specified .env file not found: $ENV_FILE"
        exit 1
    fi

    # Check for required environment variables
    if [[ -z "$OPENAI_API_KEY" ]] && [[ -n "$ENV_FILE" ]]; then
        # Try to source the .env file to check for OPENAI_API_KEY
        if grep -q "OPENAI_API_KEY" "$ENV_FILE" 2>/dev/null; then
            echo "Found OPENAI_API_KEY in .env file"
        else
            echo "Warning: OPENAI_API_KEY not found in .env file or environment"
            echo "The AI features may not work without this API key"
        fi
    elif [[ -z "$OPENAI_API_KEY" ]]; then
        echo "Warning: OPENAI_API_KEY not set in environment"
        echo "The AI features may not work without this API key"
        echo "Consider setting it in a .env file or as an environment variable"
    fi
fi

# Create temporary directory
TEMP_DIR=$(pwd)/tmp_container_explorer
mkdir -p ${TEMP_DIR}

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
echo "Project root: $PROJECT_ROOT"
echo "Results directory: $RESULTS_DIR"
echo ""
echo "PATH MAPPING:"
echo "  Host path: ${RESULTS_DIR}"
echo "  → Container path: /UORCA_results"
echo ""
echo "NOTE: When using the app, enter '/UORCA_results' in the results directory field."

# Check if port is already in use
if ss -tuln 2>/dev/null | grep -q ":${PORT} " || netstat -tuln 2>/dev/null | grep -q ":${PORT} "; then
    echo "Warning: Port $PORT appears to be in use. The app may fail to start."
    echo "Consider using a different port or stopping other services on this port."
fi

echo ""
echo "=========================================="
echo "Starting UORCA Explorer..."
echo "Container will start momentarily."
echo ""
if [[ -n "$SSH_CONNECTION" ]] || [[ -n "$SSH_CLIENT" ]]; then
    echo "ACCESS INSTRUCTIONS (SSH/Remote):"
    echo "1. On your laptop, run: ssh -L 8000:127.0.0.1:${PORT} $(whoami)@$(hostname)"
    echo "2. Open http://127.0.0.1:8000 in your browser"
else
    echo "ACCESS INSTRUCTIONS (Local):"
    echo "Open http://127.0.0.1:${PORT} in your browser"
fi
echo ""
echo "IMPORTANT NOTES:"
if [[ "$ENGINE" == "docker" ]]; then
    echo "- Docker will map container port 8501 to host port ${PORT}"
    echo "- The app URL shown by Streamlit (http://0.0.0.0:8501) maps to http://127.0.0.1:${PORT}"
else
    echo "- The app will show a URL like 'http://0.0.0.0:${PORT}' - use your access URL instead"
fi
echo "- Enter '/UORCA_results' as the results directory path in the app"
echo ""
echo "Press Ctrl+C to stop the application"
echo "=========================================="

# Cleanup function
cleanup() {
    # Prevent double execution
    if [ "$CLEANUP_DONE" -eq 1 ]; then
        return
    fi
    CLEANUP_DONE=1

    echo ""
    echo "Shutting down UORCA Explorer..."
    echo "Cleaning up temporary directory..."
    rm -rf ${TEMP_DIR} 2>/dev/null || true
    echo "UORCA Explorer stopped. Thank you for using UORCA!"
}

# Set trap to cleanup on exit and common signals
trap cleanup EXIT INT TERM

# Launch container based on engine
echo "Launching Streamlit app in container..."
echo "This may take a moment to initialize..."

if [[ "$ENGINE" == "apptainer" ]]; then
    # Apptainer execution (existing logic)
    if ! apptainer exec \
        -B ${PROJECT_ROOT}:/workspace \
        -B ${RESULTS_DIR}:/UORCA_results \
        -B ${TEMP_DIR}:/tmp \
        --tmpdir=${TEMP_DIR} \
        --cleanenv \
        $SIF_PATH \
        bash -lc "\
            cd /workspace && \
            echo 'Container started successfully. Initializing Streamlit...' && \
            export STREAMLIT_BROWSER_GATHER_USAGE_STATS=false && \
            export STREAMLIT_SERVER_HEADLESS=true && \
            uv run streamlit run main_workflow/reporting/uorca_explorer.py --server.port ${PORT} --server.address 0.0.0.0 --server.headless true
        "; then
        echo ""
        echo "Error: Failed to start UORCA Explorer with Apptainer."
        echo "Check the error messages above for troubleshooting information."
        exit 1
    fi
else
    # Docker execution
    # Build environment variable arguments
    ENV_ARGS=""

    # Add .env file if available
    if [[ -n "$ENV_FILE" ]]; then
        ENV_ARGS="--env-file $ENV_FILE"
    fi

    # Pass through important environment variables if they exist
    for var in OPENAI_API_KEY ANTHROPIC_API_KEY HUGGINGFACE_API_TOKEN; do
        if [[ -n "${!var}" ]]; then
            ENV_ARGS="$ENV_ARGS -e $var=${!var}"
        fi
    done

    if ! docker run --rm -it --init \
        -p ${PORT}:8501 \
        -v ${PROJECT_ROOT}:/workspace \
        -v ${RESULTS_DIR}:/UORCA_results:ro \
        -v ${TEMP_DIR}:/tmp \
        $ENV_ARGS \
        -e STREAMLIT_BROWSER_GATHER_USAGE_STATS=false \
        -e STREAMLIT_SERVER_HEADLESS=true \
        -e STREAMLIT_SERVER_FILE_WATCHER_TYPE=poll \
        -e PYTHONUNBUFFERED=1 \
        -w /workspace \
        "${DOCKER_IMAGE}" \
        bash -c "\
            echo 'Container started successfully. Initializing Streamlit...' && \
            uv run streamlit run main_workflow/reporting/uorca_explorer.py --server.port 8501 --server.address 0.0.0.0 --server.headless true --server.fileWatcherType poll
        "; then
        echo ""
        echo "Error: Failed to start UORCA Explorer with Docker."
        echo "Check the error messages above for troubleshooting information."
        exit 1
    fi
fi

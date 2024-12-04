#!/bin/bash

# process_geo.sh
# Description: Master script to download GEO metadata and corresponding FASTQ files.

set -euo pipefail

# Function to display usage information
usage() {
    echo "Usage: $0 -g <GEO_ACCESSION> [-n <NUM_SPOTS>] [-o <OUTPUT_DIR>] [-f]"
    echo "  -g: GEO accession number (required)"
    echo "  -n: Number of spots to download (optional, default is all)"
    echo "  -o: Output directory (optional, default is current directory)"
    echo "  -f: Overwrite existing files if they exist (optional)"
    exit 1
}

# Function to log messages
log_message() {
    local message="$1"
    echo "$(date +"%Y-%m-%d %H:%M:%S") - $message" | tee -a "$MASTER_LOG"
}

# Function to log and execute commands
log_and_execute() {
    local cmd="$1"
    log_message "Executing: $cmd"
    eval "$cmd"
}

# Initialize variables
GEO_ACCESSION=""
NUM_SPOTS=""
OUTPUT_DIR="."
FORCE=false

# Parse command-line options using getopts
while getopts "g:n:o:f" opt; do
    case "$opt" in
        g)
            GEO_ACCESSION="$OPTARG"
            ;;
        n)
            NUM_SPOTS="$OPTARG"
            ;;
        o)
            OUTPUT_DIR="$OPTARG"
            ;;
        f)
            FORCE=true
            ;;
        *)
            usage
            ;;
    esac
done

# Shift off the options and optional --
shift $((OPTIND -1))

# Validate required options
if [ -z "$GEO_ACCESSION" ]; then
    echo "Error: GEO accession is required." >&2
    usage
fi

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

# Initialize master log file
MASTER_LOG="$OUTPUT_DIR/master.log"
echo "Processing started at $(date)" > "$MASTER_LOG"
# Error logging function
log_error() {
    local message="$1"
    echo "$(date +"%Y-%m-%d %H:%M:%S") - ERROR: $message" | tee -a "$MASTER_LOG" >&2
}

# Redirect stderr to log file while preserving terminal output
exec 2> >(tee -a "$MASTER_LOG" >&2)

# Function to log messages
log_message() {
    local message="$1"
    echo "$(date +"%Y-%m-%d %H:%M:%S") - $message" | tee -a "$MASTER_LOG"
}

# Function to check if scripts exist and are executable
check_scripts() {
    # Get directory where this script is located
    SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
    # Debug output
    log_message "Debug: Script directory is $SCRIPT_DIR"
    log_message "Debug: Current working directory is $(pwd)"
    log_message "Debug: PATH is $PATH"
    local scripts=("download_metadata.R" "download_FASTQs.sh")
    for script in "${scripts[@]}"; do
        script_path="$SCRIPT_DIR/$script"
        echo "Debug: Checking for script at $script_path"
        if [ ! -f "$script_path" ]; then
            echo "Error: Required script '$script' not found in $SCRIPT_DIR" >&2
            echo "Debug: Directory contents of $SCRIPT_DIR:"
            ls -la "$SCRIPT_DIR"
            exit 1
        fi
        if [ "${script##*.}" = "sh" ] && [ ! -x "$script_path" ]; then
            echo "Debug: Making $script_path executable"
            chmod +x "$script_path"
        fi
        echo "Debug: Found and validated $script_path"
    done
}

# Function to execute the R script for metadata
execute_metadata_script() {
    log_message "Starting metadata download using download_metadata.R"
    cmd="Rscript $SCRIPT_DIR/download_metadata.R --geo_accession '$GEO_ACCESSION' --output_dir '$OUTPUT_DIR'"
    echo "Debug: Executing metadata command: $cmd"
    echo "Debug: R executable location: $(which Rscript)"
    echo "Debug: R version: $(Rscript --version 2>&1)"
    log_and_execute "$cmd" 2>> "$MASTER_LOG"
    echo "Debug: Metadata script exit code: $?"
    log_message "Metadata download completed."
}

# Function to execute the Bash script for FASTQ downloads
execute_fastqs_script() {
    log_message "Starting FASTQ download using download_fastqs.sh"
    echo "Debug: FASTQ script permissions:"
    ls -l "$SCRIPT_DIR/download_FASTQs.sh"
    cmd="$SCRIPT_DIR/download_FASTQs.sh -g '$GEO_ACCESSION' -o '$OUTPUT_DIR'"
    echo "Debug: FASTQ command to execute: $cmd"

    if [ -n "$NUM_SPOTS" ]; then
        cmd="$cmd -n $NUM_SPOTS"
    fi

    if [ "$FORCE" = true ]; then
        cmd="$cmd -f"
    fi

    log_and_execute "$cmd" 2>> "$MASTER_LOG"
    log_message "FASTQ download completed."
}

# Add error trap
trap 'echo "Debug: Error on line $LINENO. Exit code: $?"' ERR

# Main execution
main() {
    check_scripts

    log_message "Processing GEO accession: $GEO_ACCESSION"
    log_message "Output directory: $OUTPUT_DIR"
    if [ "$FORCE" = true ]; then
        log_message "Force overwrite is enabled."
    fi
    if [ -n "$NUM_SPOTS" ]; then
        log_message "Number of spots to download: $NUM_SPOTS"
    else
        log_message "Number of spots to download: All"
    fi

    # Step 1: Download Metadata
    execute_metadata_script

    # Step 2: Download FASTQ Files
    execute_fastqs_script

    log_message "All processing completed successfully for GEO accession: $GEO_ACCESSION"
}

# Run the main function
main

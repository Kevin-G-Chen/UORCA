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

# Initialize master log
MASTER_LOG="$OUTPUT_DIR/processing_master.log"
echo "Master processing started at $(date)" > "$MASTER_LOG"

# Function to check if scripts exist and are executable
check_scripts() {
    local scripts=("download_metadata.R" "download_fastqs.sh")
    for script in "${scripts[@]}"; do
        if [ ! -f "$script" ]; then
            echo "Error: Required script '$script' not found in the current directory." >&2
            exit 1
        fi
        if [ "${script##*.}" = "sh" ] && [ ! -x "$script" ]; then
            chmod +x "$script"
        fi
    done
}

# Function to execute the R script for metadata
execute_metadata_script() {
    log_message "Starting metadata download using download_metadata.R"
    cmd="Rscript download_metadata.R --geo_accession '$GEO_ACCESSION' --output_dir '$OUTPUT_DIR'"
    log_and_execute "$cmd" 2>> "$MASTER_LOG"
    log_message "Metadata download completed."
}

# Function to execute the Bash script for FASTQ downloads
execute_fastqs_script() {
    log_message "Starting FASTQ download using download_fastqs.sh"
    cmd="./download_fastqs.sh -g '$GEO_ACCESSION' -o '$OUTPUT_DIR'"

    if [ -n "$NUM_SPOTS" ]; then
        cmd="$cmd -n $NUM_SPOTS"
    fi

    if [ "$FORCE" = true ]; then
        cmd="$cmd -f"
    fi

    log_and_execute "$cmd" 2>> "$MASTER_LOG"
    log_message "FASTQ download completed."
}

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

#!/bin/bash

# download_fastqs.sh
# Description: Fetches SRA IDs for a given GEO accession and downloads corresponding FASTQ files.

set -euo pipefail

# Default value for FORCE flag
FORCE=false
OUTPUT_DIR="."
ERROR_LOG=""  # Will be initialized after OUTPUT_DIR is finalized
# Error logging function
log_error() {
    local message="$1"
    echo "$(date +"%Y-%m-%d %H:%M:%S") - ERROR: $message" | tee -a "$ERROR_LOG" >&2
}

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
    echo "$(date +"%Y-%m-%d %H:%M:%S") - $message" | tee -a "$ERROR_LOG"
}

# Function to log and execute commands
log_and_execute() {
    local cmd="$1"
    local log_output=${2:-true}  # Optional parameter to control output logging

    log_message "Executing: $cmd"
    if [ "$log_output" = true ]; then
        eval "$cmd"
    else
        # Execute without capturing output in logs
        local output
        output=$(eval "$cmd" 2>> "$ERROR_LOG")
        echo "$output"
    fi
}

# Function to check dependencies
check_dependencies() {
    local dependencies=(esearch efetch xtract fastq-dump gzip)
    for cmd in "${dependencies[@]}"; do
        if ! command -v "$cmd" &> /dev/null; then
            log_message "Error: $cmd is not installed or not in PATH."
            exit 1
        fi
    done
}

# Function to get SRA IDs for a GEO accession
get_sra_ids() {
    local GEO_ACCESSION=$1
    local OUTPUT_DIR=$2
    local SRA_ID_FILE="$OUTPUT_DIR/sra_ids.txt"

    # Fetch related GSM sample IDs
    log_message "Fetching SRA IDs for GEO accession: $GEO_ACCESSION"
    cmd="esearch -db gds -query '${GEO_ACCESSION}[ACCN] AND GSM[ETYP]' | efetch -format docsum | xtract -pattern DocumentSummary -element Accession"
    # Execute command directly instead of through log_and_execute to avoid parsing the log message
    samples=$(eval "$cmd" 2>> "$ERROR_LOG")

    if [ -z "$samples" ]; then
        log_message "No GSM samples found for GEO accession: $GEO_ACCESSION"
        exit 1
    fi

    # Initialize the SRA ID file
    echo -e "sample_ID\texperiment\tSRA_ID" > "$SRA_ID_FILE"

    # Process each GSM sample ID
        for sample in $samples; do
            # Validate sample ID format (should be GSM followed by numbers)
            if [[ ! "$sample" =~ ^GSM[0-9]+$ ]]; then
                log_message "Warning: Invalid sample ID format: $sample, skipping"
                continue
            fi

            log_message "Processing sample: $sample"

            # Get corresponding Experiment and Run IDs with explicit error handling
            cmd="esearch -db sra -query '$sample' | efetch -format docsum | xtract -pattern DocumentSummary -element Experiment@acc,Run@acc"
            ids=$(eval "$cmd" 2>> "$ERROR_LOG" | tr '\n' ',' | sed 's/,$//')
            if [ $? -ne 0 ]; then
                log_message "Error retrieving IDs for sample $sample"
                continue
            fi

        # If no IDs found, set to "No IDs found"
        if [ -z "$ids" ]; then
            ids="No_IDs_found\tNo_IDs_found"
        else
            # Replace comma with tab
            ids=$(echo "$ids" | tr ',' '\t')
        fi

        # Append to SRA ID file
        echo -e "${sample}\t${ids}" >> "$SRA_ID_FILE"
    done

    log_message "SRA IDs have been saved to $SRA_ID_FILE"
    # Convert to long format if needed
    convert_to_long_format "$SRA_ID_FILE"
}

# Function to download and process FASTQ files
download_fastqs() {
    local INPUT_FILE=$1
    local NUM_READS=$2
    local OUTPUT_DIR=$3

    # Read the file line by line, skipping the header
    tail -n +2 "$INPUT_FILE" | while IFS=$'\t' read -r sample_ID experiment SRA_ID; do
        if [[ "$SRA_ID" == "No_IDs_found" ]]; then
        log_message "Skipping sample $sample_ID due to missing SRA_ID."
            continue
        fi

        log_message "Processing $SRA_ID..."

        # Define FASTQ file patterns
        fastq_files=("${OUTPUT_DIR}/${SRA_ID}"_*.fastq)
        fastq_exists=false
        for fq in "${fastq_files[@]}"; do
            if [ -f "$fq" ]; then
                fastq_exists=true
                break
            fi
        done

        if [ "$fastq_exists" = true ] && [ "$FORCE" = false ]; then
        log_message "FASTQ files for $SRA_ID already exist. Skipping download. Use -f to overwrite."
            continue
        fi

        # Step 1: Use fast-dump to extract reads with error logging
        if [ -n "$NUM_READS" ]; then
            cmd="fastq-dump --split-files -X $NUM_READS -O $OUTPUT_DIR $SRA_ID"
            if ! log_and_execute "$cmd" 2>> "$ERROR_LOG"; then
                log_message "Error: fastq-dump failed for $SRA_ID"
                exit 1
            fi
        else
            if ! fastq-dump --split-files -O "$OUTPUT_DIR" "$SRA_ID" 2>> "$ERROR_LOG"; then
            log_message "Error: fastq-dump failed for $SRA_ID. Check $ERROR_LOG for details."
                exit 1
            fi
        fi

        # Step 2: Sort the FASTQ files (Use with caution)
        for fastq_file in "${OUTPUT_DIR}/${SRA_ID}"_*.fastq; do
            if [ -f "$fastq_file" ]; then
                sorted_file="${fastq_file}.sorted"
                if [ -f "$sorted_file" ] && [ "$FORCE" = false ]; then
                log_message "Sorted file $sorted_file already exists. Skipping sorting. Use -f to overwrite."
                    continue
                fi

                log_message "Sorting $fastq_file..."
                paste - - - - < "$fastq_file" | sort -k1,1 -t " " | tr "\t" "\n" > "${sorted_file}"
                mv "${sorted_file}" "$fastq_file"
            fi
        done

        # Step 3: Gzip the FASTQ files
        for fastq_file in "${OUTPUT_DIR}/${SRA_ID}"_*.fastq; do
            if [ -f "$fastq_file" ]; then
                gz_file="${fastq_file}.gz"
                if [ -f "$gz_file" ] && [ "$FORCE" = false ]; then
                log_message "Gzipped file $gz_file already exists. Skipping gzipping. Use -f to overwrite."
                    continue
                fi

                log_message "Gzipping $fastq_file..."
                gzip -f "$fastq_file"
            fi
        done
    done
}


# Function to convert sra_ids file to long format
convert_to_long_format() {
    local INPUT_FILE=$1
    local TEMP_FILE="${INPUT_FILE}.temp"
    local HEADER="sample_ID\texperiment\tSRA_ID"

    # Write header to temp file
    echo -e "$HEADER" > "$TEMP_FILE"

    # Process each line (skipping header)
    tail -n +2 "$INPUT_FILE" | while IFS=$'\t' read -r sample_id experiment srr_ids; do
        # Split SRR IDs on whitespace and create a row for each
        for srr_id in $srr_ids; do
            echo -e "${sample_id}\t${experiment}\t${srr_id}" >> "$TEMP_FILE"
        done
    done

    # Replace original file with temp file
    mv "$TEMP_FILE" "$INPUT_FILE"
    log_message "Converted SRA IDs file to long format"
}

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

# Check dependencies
check_dependencies

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"
# Initialize single log file
ERROR_LOG="$OUTPUT_DIR/processing.log"
echo "Processing started at $(date)" > "$ERROR_LOG"


# Trap to clean up on errors
trap 'log_error "An error occurred on line $LINENO: $BASH_COMMAND"; exit 1' ERR

# Redirect stderr to log file while preserving terminal output
exec 2> >(tee -a "$ERROR_LOG" >&2)

# Step 1: Get SRA IDs
get_sra_ids "$GEO_ACCESSION" "$OUTPUT_DIR"

# Step 2: Download and process FASTQ files
echo "Downloading and processing FASTQ files..." >&2
download_fastqs "$OUTPUT_DIR/sra_ids.txt" "$NUM_SPOTS" "$OUTPUT_DIR"

log_message "FASTQ download and processing completed for $GEO_ACCESSION"

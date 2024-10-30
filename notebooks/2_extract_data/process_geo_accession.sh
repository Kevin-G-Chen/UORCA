#!/bin/bash

# Initialize error log will be moved later after OUTPUT_DIR is set

set -euo pipefail

# Function to check dependencies
check_dependencies() {
    local dependencies=(esearch efetch xtract fasterq-dump gzip)
    for cmd in "${dependencies[@]}"; do
        if ! command -v "$cmd" &> /dev/null; then
            echo "Error: $cmd is not installed or not in PATH." >&2
            exit 1
        fi
    done
}

# Function to get SRA IDs for a GEO accession
get_sra_ids() {
    local GEO_ACCESSION=$1
    local OUTPUT_DIR=$2

    local OUTPUT_FILE="$OUTPUT_DIR/results.txt"

    # Initialize the output file with headers
    echo -e "sample_ID\texperiment\tSRA_ID" > "$OUTPUT_FILE"

    # Fetch related GSM sample IDs
    samples=$(esearch -db gds -query "${GEO_ACCESSION}[ACCN] AND GSM[ETYP]" | \
              efetch -format docsum | \
              xtract -pattern DocumentSummary -element Accession)

    # Process each GSM sample ID
    for sample in $samples; do
        echo "Processing sample: $sample" >&2

        # Get corresponding Experiment and Run IDs
        ids=$(esearch -db sra -query "$sample" | \
              efetch -format docsum | \
              xtract -pattern DocumentSummary -element Experiment@acc,Run@acc | \
              tr '\n' ',' | sed 's/,$//')

        # If no IDs found, set to "No IDs found"
        if [ -z "$ids" ]; then
            ids="No IDs found\tNo IDs found"
        else
            # Replace comma with tab
            ids=$(echo "$ids" | tr ',' '\t')
        fi

        # Append to output file
        echo -e "${sample}\t${ids}" >> "$OUTPUT_FILE"
    done
}

# Function to download and process FASTQ files
download_fastqs() {
    local INPUT_FILE=$1
    local NUM_READS=$2
    local OUTPUT_DIR=$3

    # Read the file line by line, skipping the header
    tail -n +2 "$INPUT_FILE" | while IFS=$'\t' read -r sample_ID experiment SRA_ID; do
        if [[ "$SRA_ID" == "No IDs found" ]]; then
            echo "Skipping sample $sample_ID due to missing SRA_ID." >&2
            continue
        fi

        # Validate SRA ID format
        if [[ ! "$SRA_ID" =~ ^SRR[0-9]{8,}$ ]]; then
            echo "Warning: Invalid SRA ID format '$SRA_ID'. Skipping." >&2
            continue
        fi

        echo "Processing $SRA_ID..." >&2

        # Step 1: Use fastq-dump to extract reads with error logging
        if [ -n "$NUM_READS" ]; then
            if ! fastq-dump --split-files -X "$NUM_READS" -O "$OUTPUT_DIR" "$SRA_ID" 2>> "$ERROR_LOG"; then
                echo "Error: fastq-dump failed for $SRA_ID. Check $ERROR_LOG for details." >&2
                exit 1
            fi
        else
            if ! fastq-dump --split-files -O "$OUTPUT_DIR" "$SRA_ID" 2>> "$ERROR_LOG"; then
                echo "Error: fastq-dump failed for $SRA_ID. Check $ERROR_LOG for details." >&2
                exit 1
            fi
        fi

        # Step 2: Sort the FASTQ files (Use with caution)
        for fastq_file in "${OUTPUT_DIR}/${SRA_ID}"_*.fastq; do
            if [ -f "$fastq_file" ]; then
                echo "Sorting $fastq_file..." >&2
                paste - - - - < "$fastq_file" | sort -k1,1 -t " " | tr "\t" "\n" > "${fastq_file}.sorted"
                mv "${fastq_file}.sorted" "$fastq_file"
            fi
        done

        # Step 3: Gzip the FASTQ files
        gzip "${OUTPUT_DIR}/${SRA_ID}"_*.fastq
    done
}

# Function to display usage information
usage() {
    echo "Usage: $0 -g <GEO_ACCESSION> [-n <NUM_SPOTS>] [-o <OUTPUT_DIR>]"
    echo "  -g: GEO accession number (required)"
    echo "  -n: Number of spots to download (optional, default is all)"
    echo "  -o: Output directory (optional, default is current directory)"
    exit 1
}

# Main script execution starts here

# Check dependencies
check_dependencies

# Parse command-line options
while getopts ":g:n:o:" opt; do
    case $opt in
        g) GEO_ACCESSION="$OPTARG" ;;
        n) NUM_SPOTS="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        \?) echo "Invalid option -$OPTARG" >&2; usage ;;
        :) echo "Option -$OPTARG requires an argument." >&2; usage ;;
    esac
done

# Check if GEO accession is provided
if [ -z "${GEO_ACCESSION:-}" ]; then
    echo "Error: GEO accession is required." >&2
    usage
fi

# Set default values if not provided
OUTPUT_DIR=${OUTPUT_DIR:-.}
NUM_SPOTS=${NUM_SPOTS:-}

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

# Initialize error log
ERROR_LOG="$OUTPUT_DIR/processing.log"
echo "Processing started at $(date)" > "$ERROR_LOG"

export TERM=${TERM:-dumb}

# Trap to clean up on errors
trap 'echo "An error occurred. Check the logs for details." >&2; exit 1' ERR

# Step 1: Get SRA IDs
echo "Getting SRA IDs for $GEO_ACCESSION..." >&2
get_sra_ids "$GEO_ACCESSION" "$OUTPUT_DIR"

# Step 2: Download and process FASTQ files
echo "Downloading and processing FASTQ files..." >&2
download_fastqs "$OUTPUT_DIR/results.txt" "$NUM_SPOTS" "$OUTPUT_DIR"

echo "Processing completed for $GEO_ACCESSION" >&2
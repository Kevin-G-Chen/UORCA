#!/bin/bash

# download_fastqs.sh
# Description: Fetches SRA IDs for a given GEO accession and downloads corresponding FASTQ files.

set -euo pipefail

# Default value for FORCE flag
FORCE=false

# Function to display usage information
# Function to log and execute commands
log_and_execute() {
    local cmd="$1"
    echo "$(date +'%Y-%m-%d %H:%M:%S') - Executing: $cmd" >> "$ERROR_LOG"
    eval "$cmd"
}
usage() {
    echo "Usage: $0 -g <GEO_ACCESSION> [-n <NUM_SPOTS>] [-o <OUTPUT_DIR>] [--force]"
    echo "  -g, --geo_accession: GEO accession number (required)"
    echo "  -n, --num_spots: Number of spots to download (optional, default is all)"
    echo "  -o, --output_dir: Output directory (optional, default is current directory)"
    echo "  -f, --force: Overwrite existing files if they exist (optional)"
    exit 1
}

# Function to check dependencies
check_dependencies() {
    local dependencies=(esearch efetch xtract fastq-dump gzip)
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
    local SRA_ID_FILE="$OUTPUT_DIR/sra_ids.txt"

    # Fetch related GSM sample IDs
    echo "Fetching SRA IDs for GEO accession: $GEO_ACCESSION" >&2
    cmd="esearch -db gds -query '${GEO_ACCESSION}[ACCN] AND GSM[ETYP]' | efetch -format docsum | xtract -pattern DocumentSummary -element Accession"
    samples=$(log_and_execute "$cmd")

    if [ -z "$samples" ]; then
        echo "No GSM samples found for GEO accession: $GEO_ACCESSION" >&2
        exit 1
    fi

    # Initialize the SRA ID file
    echo -e "sample_ID\texperiment\tSRA_ID" > "$SRA_ID_FILE"

    # Process each GSM sample ID
    for sample in $samples; do
        echo "Processing sample: $sample" >&2

        # Get corresponding Experiment and Run IDs
        cmd="esearch -db sra -query '$sample' | efetch -format docsum | xtract -pattern DocumentSummary -element Experiment@acc,Run@acc | tr '\n' ',' | sed 's/,\$//'"
        ids=$(log_and_execute "$cmd")

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

    echo "SRA IDs have been saved to $SRA_ID_FILE" >&2
}

# Function to download and process FASTQ files
download_fastqs() {
    local INPUT_FILE=$1
    local NUM_READS=$2
    local OUTPUT_DIR=$3
    local ERROR_LOG="$OUTPUT_DIR/processing.log"

    # Read the file line by line, skipping the header
    tail -n +2 "$INPUT_FILE" | while IFS=$'\t' read -r sample_ID experiment SRA_ID; do
        if [[ "$SRA_ID" == "No_IDs_found" ]]; then
            echo "Skipping sample $sample_ID due to missing SRA_ID." >&2
            continue
        fi

        # Validate SRA ID format
        if [[ ! "$SRA_ID" =~ ^SRR[0-9]{8,}$ ]]; then
            echo "Warning: Invalid SRA ID format '$SRA_ID'. Skipping." >&2
            continue
        fi

        echo "Processing $SRA_ID..." >&2

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
            echo "FASTQ files for $SRA_ID already exist. Skipping download. Use --force to overwrite." >&2
            continue
        fi

        # Step 1: Use fast-dump to extract reads with error logging
        if [ -n "$NUM_READS" ]; then
        cmd="fastq-dump --split-files -X $NUM_READS -O $OUTPUT_DIR $SRA_ID"
        if ! log_and_execute "$cmd" 2>> "$ERROR_LOG"; then
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
                sorted_file="${fastq_file}.sorted"
                if [ -f "$sorted_file" ] && [ "$FORCE" = false ]; then
                    echo "Sorted file $sorted_file already exists. Skipping sorting. Use --force to overwrite." >&2
                    continue
                fi

                echo "Sorting $fastq_file..." >&2
                paste - - - - < "$fastq_file" | sort -k1,1 -t " " | tr "\t" "\n" > "${sorted_file}"
                mv "${sorted_file}" "$fastq_file"
            fi
        done

        # Step 3: Gzip the FASTQ files
        for fastq_file in "${OUTPUT_DIR}/${SRA_ID}"_*.fastq; do
            if [ -f "$fastq_file" ]; then
                gz_file="${fastq_file}.gz"
                if [ -f "$gz_file" ] && [ "$FORCE" = false ]; then
                    echo "Gzipped file $gz_file already exists. Skipping gzipping. Use --force to overwrite." >&2
                    continue
                fi

                echo "Gzipping $fastq_file..." >&2
                gzip -f "$fastq_file"
            fi
        done
    done
}

# Parse command-line options using getopt
TEMP=$(getopt -o g:n:o:f --long geo_accession:,num_spots:,output_dir:,force -n 'download_fastqs.sh' -- "$@")
if [ $? != 0 ]; then
    echo "Error: Failed to parse options." >&2
    usage
fi

eval set -- "$TEMP"

# Initialize variables
GEO_ACCESSION=""
NUM_SPOTS=""
OUTPUT_DIR="."

# Extract options
while true; do
    case "$1" in
        -g|--geo_accession)
            GEO_ACCESSION="$2"
            shift 2
            ;;
        -n|--num_spots)
            NUM_SPOTS="$2"
            shift 2
            ;;
        -o|--output_dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -f|--force)
            FORCE=true
            shift
            ;;
        --)
            shift
            break
            ;;
        *)
            echo "Error: Invalid option $1" >&2
            usage
            ;;
    esac
done

# Validate required options
if [ -z "$GEO_ACCESSION" ]; then
    echo "Error: GEO accession is required." >&2
    usage
fi

# Check dependencies
check_dependencies

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

# Initialize error log
ERROR_LOG="$OUTPUT_DIR/processing.log"
echo "Processing started at $(date)" > "$ERROR_LOG"

# Trap to clean up on errors
trap 'echo "An error occurred. Check $ERROR_LOG for details." >&2; exit 1' ERR

# Step 1: Get SRA IDs
get_sra_ids "$GEO_ACCESSION" "$OUTPUT_DIR"

# Step 2: Download and process FASTQ files
echo "Downloading and processing FASTQ files..." >&2
download_fastqs "$OUTPUT_DIR/sra_ids.txt" "$NUM_SPOTS" "$OUTPUT_DIR"

echo "FASTQ download and processing completed for $GEO_ACCESSION" >&2
#!/bin/bash

# Function to display usage information
usage() {
    echo "Usage: $0 -g <GEO_ACCESSION> [-n <NUM_SPOTS>] [-o <OUTPUT_DIR>]"
    echo "  -g: GEO accession number (required)"
    echo "  -n: Number of spots to download (optional, default is all)"
    echo "  -o: Output directory (optional, default is current directory)"
    exit 1
}

# Parse command-line options
while getopts ":g:n:o:" opt; do
    case $opt in
        g) GEO_ACCESSION="$OPTARG" ;;
        n) NUM_SPOTS="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        \?) echo "Invalid option -$OPTARG" >&2; usage ;;
    esac
done

# Check if GEO accession is provided
if [ -z "$GEO_ACCESSION" ]; then
    echo "Error: GEO accession is required."
    usage
fi

# Set default values if not provided
OUTPUT_DIR=${OUTPUT_DIR:-.}
NUM_SPOTS=${NUM_SPOTS:-0}  # 0 means download all spots

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

# Step 1: Get SRA IDs
echo "Getting SRA IDs for $GEO_ACCESSION..."
./get_sra_ids.sh "$GEO_ACCESSION" "$OUTPUT_DIR"

# Step 2: Download and process FASTQ files
echo "Downloading and processing FASTQ files..."
if [ $NUM_SPOTS -eq 0 ]; then
    ./download_FASTQs.sh -i "$OUTPUT_DIR/results.txt" -o "$OUTPUT_DIR"
else
    ./download_FASTQs.sh -i "$OUTPUT_DIR/results.txt" -n "$NUM_SPOTS" -o "$OUTPUT_DIR"
fi

echo "Processing completed for $GEO_ACCESSION"

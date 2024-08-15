#!/bin/bash

# Usage function
usage() {
    echo "Usage: $0 -i <input_file> -n <num_reads> -o <output_dir>"
    echo "  -i: Path to the input file (e.g., your_file.txt)"
    echo "  -n: Number of reads to include from each FASTQ file"
    echo "  -o: Output directory to save the processed FASTQ files"
    exit 1
}

# Parse command-line arguments
while getopts ":i:n:o:" opt; do
    case "${opt}" in
        i) input_file="${OPTARG}" ;;
        n) num_reads="${OPTARG}" ;;
        o) output_dir="${OPTARG}" ;;
        *) usage ;;
    esac
done

# Check if all arguments are provided
if [ -z "${input_file}" ] || [ -z "${num_reads}" ] || [ -z "${output_dir}" ]; then
    usage
fi

# Ensure the output directory exists
mkdir -p "${output_dir}"

# Read the file line by line, skipping the header
tail -n +2 "$input_file" | while IFS=$'\t' read -r sample_ID experiment SRA_ID; do
    echo "Processing $SRA_ID..."

    # Step 1: Use fastq-dump to extract a limited number of reads and save to the output directory
    fastq-dump --split-files -X "${num_reads}" -O "${output_dir}" "$SRA_ID"
    
    # Step 2: Sort the FASTQ files
    for fastq_file in "${output_dir}/${SRA_ID}_1.fastq" "${output_dir}/${SRA_ID}_2.fastq"; do
        if [ -f "$fastq_file" ]; then
            paste - - - - < "$fastq_file" | sort -k1,1 -t " " | tr "\t" "\n" > "${fastq_file}.sorted"
            mv "${fastq_file}.sorted" "$fastq_file"
        fi
    done
    
    # Step 3: Gzip the FASTQ files
    gzip "${output_dir}/${SRA_ID}_1.fastq" "${output_dir}/${SRA_ID}_2.fastq"
    
done

echo "Processing completed."
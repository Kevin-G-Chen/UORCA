#!/bin/bash

# Usage function
usage() {
    echo "Usage: $0 -i <input_file> -n <num_lines> -o <output_dir>"
    echo "  -i: Path to the input file (e.g., your_file.txt)"
    echo "  -n: Number of lines to include from each FASTQ file"
    echo "  -o: Output directory to save the processed FASTQ files"
    exit 1
}

# Parse command-line arguments
while getopts ":i:n:o:" opt; do
    case "${opt}" in
        i) input_file="${OPTARG}" ;;
        n) num_lines="${OPTARG}" ;;
        o) output_dir="${OPTARG}" ;;
        *) usage ;;
    esac
done

# Check if all arguments are provided
if [ -z "${input_file}" ] || [ -z "${num_lines}" ] || [ -z "${output_dir}" ]; then
    usage
fi

# Ensure the output directory exists
mkdir -p "${output_dir}"

# Read the file line by line, skipping the header
tail -n +2 "$input_file" | while IFS=$t read -r sample_ID experiment SRA_ID; do
    echo "Processing $SRA_ID..."

    # Step 1: Prefetch the SRA data
    prefetch "$SRA_ID"
    
    # Step 2: Convert SRA to FASTQ
    fasterq-dump --split-files "$SRA_ID" -O "${output_dir}"
    
    # Step 3: Subset the FASTQ files
    for fastq_file in "${output_dir}/${SRA_ID}_1.fastq" "${output_dir}/${SRA_ID}_2.fastq"; do
        if [ -f "$fastq_file" ]; then
            head -n "${num_lines}" "$fastq_file" > "${fastq_file}.tmp" && mv "${fastq_file}.tmp" "$fastq_file"
        fi
    done
    
    # Step 4: Sort the FASTQ files
    for fastq_file in "${output_dir}/${SRA_ID}_1.fastq" "${output_dir}/${SRA_ID}_2.fastq"; do
        if [ -f "$fastq_file" ]; then
            paste - - - - < "$fastq_file" | sort -k1,1 -t " " | tr "\t" "\n" > "${fastq_file}.sorted"
            mv "${fastq_file}.sorted" "$fastq_file"
        fi
    done
    
    # Step 5: Gzip the FASTQ files
    gzip "${output_dir}/${SRA_ID}_1.fastq" "${output_dir}/${SRA_ID}_2.fastq"
    
    # Step 6: Clean up intermediate files (if any)
    rm -f "$SRA_ID".sra
done

echo "Processing completed."


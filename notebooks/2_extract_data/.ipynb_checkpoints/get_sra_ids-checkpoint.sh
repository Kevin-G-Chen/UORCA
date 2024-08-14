#!/bin/bash

# Ensure that the necessary tools are in the PATH
export PATH=$PATH:/home/myuser/edirect

# Check if a GEO accession was provided, otherwise display an error and exit
if [ -z "$1" ]; then
  echo "Error: No GEO accession provided."
  echo "Usage: $0 <GEO_ACCESSION>"
  exit 1
fi

# Set the GEO accession to the provided argument
GEO_ACCESSION=$1

# Initialize the output file with headers
echo -e "sample_ID\texperiment\tSRA_ID" > results.txt

# Print a message indicating the GEO accession being processed
echo "Processing GEO accession: $GEO_ACCESSION"

# Fetch related GSM sample IDs
samples=$(esearch -db gds -query "$GEO_ACCESSION[ACCN] AND GSM[ETYP]" | \
           efetch -format docsum | \
           xtract -pattern DocumentSummary -element Accession)

# Process each GSM sample ID
for sample in $samples; do
  # Print a message indicating the current sample being processed
  echo "Processing sample: $sample" 2>/dev/null

  # Run the esearch and efetch commands to get the corresponding Experiment and Run IDs
  ids=$(esearch -db sra -query "$sample" | efetch -format docsum | xtract -pattern DocumentSummary -element Experiment@acc,Run@acc)

  # If no IDs are found, set ids to "No IDs found"
  if [ -z "$ids" ]; then
    ids="No IDs found\tNo IDs found"
  fi

  # Append the sample ID, Experiment ID, and Run ID to the results.txt file
  echo -e "$sample\t$ids" >> results.txt
done
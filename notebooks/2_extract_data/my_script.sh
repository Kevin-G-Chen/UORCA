#!/bin/bash

# Initialize the output file with headers
echo -e "experiment\tSRA_ID" > results.txt

# Use a for loop to iterate over each line in output.txt
for experiment in $(cat output.txt); do
  # Print a message indicating the current experiment being processed
  echo "Processing experiment: $experiment"

  # Run the esearch and efetch commands to get the corresponding SRA ID
  sra_id=$(esearch -db sra -query "$experiment" | efetch -format docsum | xtract -pattern DocumentSummary -element Run@acc)

  # If no SRA ID is found, set sra_id to "No SRA ID found"
  if [ -z "$sra_id" ]; then
    sra_id="No SRA ID found"
  fi

  # Append the experiment and SRA ID to the results.txt file
  echo -e "$experiment\t$sra_id" >> results.txt
done


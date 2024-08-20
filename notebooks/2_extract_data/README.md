The main scripts here:
- get_sra_ids.sh is a bash script where, given a GEO accession, it will identify the SRA IDs corresponding to the samples nested in the dataset accession. The results are associated with results.txt (I will probably modify this script to produce a slightly more descriptive result...)
- download_FASTQs.sh is a bash script that will download a subset of FASTQ files associated with SRA IDs. It takes as input the result from the above script.

The other files:
- InteractiveBashTests contains some testing I do, exploring other options for FASTQ downloads
- the R scripts contain my metadata download tests. The execution of one of these scripts is included in ExtractNCBIGEOData.

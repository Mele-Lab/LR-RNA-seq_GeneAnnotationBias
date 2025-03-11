#!/bin/bash
# This script extracts fasta sequences for unique genomic regions 
# that exceed a specified size threshold using BEDTools.

# Set the threshold size for unique regions
thr_size=1000

# Load the BEDTools module
module load BEDTools/2.30.0-GCC-11.3.0

# Iterate over each file containing unique regions exceeding the threshold size
for file in $(ls */unique_regions_genome_sup${thr_size}.bed); do
    # Get the directory path of the current file
    path="$(dirname $file)"
    echo "Processing file: $file"

    # Extract the fasta name from the directory path
    fasta_name=$(echo $path | cut -d"_" -f1)

    # Define the path to the corresponding fasta file
    fasta_path="[PATH]/data/genome/${fasta_name}/${fasta_name}.fa"
    echo "Using fasta: $fasta_path"

    # Use BEDTools to extract fasta sequences for the unique regions
    bedtools getfasta \
    -fi $fasta_path \
    -bed $file \
    -fo "${path}/unique_regions_genome_sup${thr_size}.fasta"
done


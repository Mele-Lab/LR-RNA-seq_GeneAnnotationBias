#!/bin/bash

# Check if the correct number of arguments are provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <gff_file.gz>"
    exit 1
fi

# Assign variables
gff_file=$1
contig_map_file=data/ref_contigs_correspondence.tsv
output_file="${gff_file%.*}_updated.gff.gz"

# Create a temporary sed script file
sed_script=$(mktemp)

# Generate the sed script from the contig map
while IFS=$'\t' read -r original_contig new_contig; do
    echo "s/\b${original_contig}\b/${new_contig}/g" >> "$sed_script"
done < "$contig_map_file"

# Process the gzipped input file: filter out lines containing "SIRV" or "ERCC", substitute contig names, and compress the output
zcat "$gff_file" | grep -v -e "SIRV" -e "ERCC" -e "chrEBV" | sed -f "$sed_script" | gzip > "$output_file"

# Clean up
rm "$sed_script"

echo "Contig names substituted and filtered successfully. Output file: $output_file"

#!/bin/bash

# add column with sample anem


# Define input directory
INPUTDIR="data/poder"  # Replace with your actual input directory

# Loop through matching files
for file in "$INPUTDIR"/*/*.transcript_counts.tsv; do
    # Skip if no files match
    [ -e "$file" ] || continue

    # Extract the base filename portion before .transcript_counts.tsv
    base=$(basename "$file" .transcript_counts.tsv)

    # Add the new column to the file
    # Assuming the file is tab-delimited, use awk to append the column
    awk -v col="$base" 'BEGIN {OFS="\t"} {print $0, col}' "$file" > $INPUTDIR/${base}_stringent_samplecol.tsv


done

# concatenate
awk 'FNR==1 && NR!=1 { next } { print }' $INPUTDIR/*_stringent_samplecol.tsv > $INPUTDIR/concat_isoquant_counts.tsv

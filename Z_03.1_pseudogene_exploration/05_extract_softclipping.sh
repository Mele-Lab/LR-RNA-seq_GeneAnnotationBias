#!/bin/bash

BED=$1
SAMPLE=$(echo $INPUTFASTQ | sed 's-.*/--' | sed 's/_preprocessed\.fastq\.gz//')

# Use samtools to convert BAM to SAM and parse CIGAR strings
cat $BED| awk '{
    cigar = $7;
    # Initialize S count
    S_count = 0;
    # Iterate through the CIGAR string
    while (match(cigar, /[0-9]+S/)) {
        # Extract the number of S bases
        S_count += substr(cigar, RSTART, RLENGTH - 1);
        # Remove the processed part of the CIGAR string
        cigar = substr(cigar, RSTART + RLENGTH);
    }
    # Print the total S count for this line
    print S_count;
}' > $TMPDIR/$SAMPLE.softclipping

paste $BED $TMPDIR/$SAMPLE.softclipping > data_unmasked_unfiltered/${SAMPLE}_softclipping.bed
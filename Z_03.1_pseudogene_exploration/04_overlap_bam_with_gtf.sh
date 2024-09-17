#!/bin/bash

module load bedtools samtools
PATH=/gpfs/projects/bsc83/utils/bedops:$PATH

GTF=/gpfs/projects/bsc83/Data/gene_annotations/gencode/v47/modified/gencode.v47.primary_assembly.sirvset4.all_genes.bed
BAM=$1

SAMPLE=$(echo $BAM | sed 's-.*/--' | sed 's/\.bam//')



# samtools view $BAM | awk '
# {
#     cigar = $6;
#     sum = 0;
#     while (match(cigar, /[0-9]+/)) {
#         sum += substr(cigar, RSTART, RLENGTH);
#         cigar = substr(cigar, RSTART + RLENGTH);
#     }
#     printf "%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $3, $4-1, $4+sum-1, $1, $2,$5,$6,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22;

# }' > data_unmasked_unfiltered/$SAMPLE.bed



echo "BED intersect\n\n"

# # intersect
bedtools intersect -a data_unmasked_unfiltered/$SAMPLE.bed -b $GTF -wa -wb -f 0.7 -F 0.9 -e > $TMPDIR/${SAMPLE}_overlapped.bed

# echo "parse S\n\n"

# retrieve S and paste them 
cat $TMPDIR/${SAMPLE}_overlapped.bed | awk '{
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

echo "Paste\n\n"

paste $TMPDIR/$SAMPLE.softclipping $TMPDIR/${SAMPLE}_overlapped.bed > data_unmasked_unfiltered/${SAMPLE}_softclipping.bed
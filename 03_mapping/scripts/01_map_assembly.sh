#!/bin/bash

####### README
# This script maps fastqs of Q10 without prioritizing annotated junctions 
#   and keeps only primary alignments to genome contigs
#   This are the most similar settings to minimap for lyric


module load minimap2
module load samtools

INPUTFASTQ=$1
TYPE=assembly_mapping
INDEX=ref/GRCh38.primary_assembly.sirvset4.genome.mmi
# automatically parse previous info for the outputs
SAMPLE=$(echo $INPUTFASTQ | sed 's-.*/--' | sed 's/_preprocessed_Q10\.fastq\.gz//')
OUTSAM=data/$TYPE/$SAMPLE.sam
OUTBAM=data/$TYPE/genomic/$SAMPLE.bam
OUTSIRVBAM=data/$TYPE/sirv/${SAMPLE}_SIRV.bam


# Create dirs
mkdir data/$TYPE
mkdir data/$TYPE/genomic
mkdir data/$TYPE/sirv

# Map
minimap2 \
    -ax splice \
    -t 112 \
    --MD \
    --secondary=no \
    -L \
    -o $OUTSAM \
    -a $INDEX \
    $INPUTFASTQ

# SAM to BAM
# Keep only reads mapped to genome contigs
# Filter out non-primary alignments, unmapped reads and supplementary alignments
# Sort BAM
samtools view \
    -@112 \
    -b \
    -F 256 \
    -F 4 \
    -F 2048 \
    -L ref/GRCh38_primary_assembly_contigs.bed $OUTSAM |\
    samtools sort -@112 -o $OUTBAM

# Index BAM
samtools index -@112 -b $OUTBAM

# Now the same for SIRV mapping reads
samtools view \
    -@112 \
    -b \
    -F 256 \
    -F 4 \
    -F 2048 \
    -L ref/sirv_erc_contigs.bed $OUTSAM |\
    samtools sort -@112 -o $OUTSIRVBAM

samtools index -@112 -b $OUTSIRVBAM

# Remove SAM
rm $OUTSAM


#!/bin/bash

####### README
# This script maps fastqs of Q7 prioritizing annotated junctions and filters by MAPQ10
module load minimap2
module load samtools

INPUTFASTQ=$1
TYPE=unmerged
JUNCBED=/gpfs/projects/bsc83/Data/gene_annotations/gencode/v47/modified/gencode.v47.primary_assembly.sirvset4.annotation.junc.bed
INDEX=ref/GRCh38.primary_assembly.sirvset4.genome.mmi
# automatically parse previous info for the outputs
SAMPLE=$(echo $INPUTFASTQ | sed 's-.*/--' | sed 's/_preprocessed_Q7\.fastq\.gz//')
OUTSAM=data/$TYPE/$SAMPLE.sam
OUTBAM=data/$TYPE/genomic/$SAMPLE.bam
OUTSIRVBAM=data/$TYPE/sirv/${SAMPLE}_SIRV.bam

# create dirs
mkdir data/$TYPE
mkdir data/$TYPE/genomic
mkdir data/$TYPE/sirv

# map
minimap2 \
    -ax splice \
    --junc-bed $JUNCBED \
    -t 112 \
    --MD \
    -o $OUTSAM \
    -a $INDEX \
    $INPUTFASTQ


# SAM to BAM
# Keep only reads mapped to genome contigs
# Filter out alignments with MAPQ<10
# Sort BAM
samtools view \
    -@112 \
    -b \
    --min-MQ 10 \
    -L ref/GRCh38_primary_assembly_contigs.bed $OUTSAM |\
    samtools sort -@112 -o $OUTBAM

# Index BAM
samtools index -@112 -b $OUTBAM

# Now the same for SIRV mapping reads
samtools view \
    -@112 \
    -b \
    --min-MQ 10 \
    -L ref/sirv_erc_contigs.bed $OUTSAM |\
    samtools sort -@112 -o $OUTSIRVBAM

samtools index -@112 -b $OUTSIRVBAM

# Remove SAM
rm $OUTSAM
#!/bin/bash

module load samtools

ITERATIONS=10
NUM_READS=5000000
SAMPLE=$1
SAMPLENAME=$(echo $1 | sed 's-.*/--' | sed 's/\.bam//')

mkdir data/bam/
mkdir data/bam/temp

# keep header of the bam file
samtools view -@112 -H $SAMPLE > data/bam/temp/${SAMPLENAME}_header.sam

# keep beheaded bam file in sam format
samtools view -@112 $SAMPLE > data/bam/temp/${SAMPLENAME}_beheaded.sam

for iteration in $(seq 1 $ITERATIONS)
do
    # Select X random lines of the sam
    shuf -n $NUM_READS data/bam/temp/${SAMPLENAME}_beheaded.sam \
        > data/bam/temp/${SAMPLENAME}_beheaded_downsampled_iter${iteration}_${NUM_READS}.sam

    # Rehead the random sam
    cat data/bam/temp/${SAMPLENAME}_header.sam data/bam/temp/${SAMPLENAME}_beheaded_downsampled_iter${iteration}_${NUM_READS}.sam \
        > data/bam/temp/${SAMPLENAME}_downsampled_iter${iteration}_${NUM_READS}.sam
    
    # SAM to BAM
    samtools view -@ 112 -b data/bam/temp/${SAMPLENAME}_downsampled_iter${iteration}_${NUM_READS}.sam \
        -o data/bam/${SAMPLENAME}_downsampled_iter${iteration}_${NUM_READS}reads.bam
done
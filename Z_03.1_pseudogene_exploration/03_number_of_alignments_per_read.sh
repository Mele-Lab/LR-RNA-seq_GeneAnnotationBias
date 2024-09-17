#!/bin/bash

module load samtools

SAMPLE=$1
SAMPLENAME=$(echo $1 | sed 's-.*/--' | sed 's/\.bam//')


samtools view -@112 $SAMPLE | cut -f1 | sort | uniq -c > data_counts/${SAMPLENAME}_unmasked.tsv

cat data_counts/${SAMPLENAME}_unmasked.tsv | awk 'BEGIN{OFS="\t"}{print $1,$2,"unmasked"}' > data_counts/${SAMPLENAME}_unmasked.tsv_mod
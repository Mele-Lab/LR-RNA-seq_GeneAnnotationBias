#!/bin/bash
module load samtools

UNALIGNED_BAM=$1 #01_basecalling/data/20240212_HS_2_PY2_GM10493.bam
OUTPUT_BAM=$2
samtools view $UNALIGNED_BAM | awk '{$1=$1"_"$13; print $0}' > $OUTPUT_BAM

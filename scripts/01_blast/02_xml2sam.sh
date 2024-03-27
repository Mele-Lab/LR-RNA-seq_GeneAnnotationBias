#!/bin/bash

module load samtools
# Input file

XML=$1
REF=$2 # fasta of the reference used in blast
READS=$3 # fasta
NAME=$(echo "$XML" | sed 's/.*\///' | sed 's/\.xml//')
RPATH="01_blast/01_results/"$4

echo $RPATH
# Transform xml format to sam
/gpfs/projects/bsc83/utils/Blast2Bam/bin/blast2bam --minAlignLength 44 $XML $REF $READS > ${RPATH}.sam

echo "transformed to sam"

# filter out unmapped reads
awk '$2!=4 {print}' ${RPATH}.sam > $RPATH.filtered.sam

# sam2bam, sort and index sam
samtools view -b $RPATH.filtered.sam > $RPATH.bam
samtools sort $RPATH.bam > $RPATH.sorted.bam
samtools index $RPATH.sorted.bam > $RPATH.sorted.bam.bai

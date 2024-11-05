#!/bin/bash

module load anaconda bedtools
source activate base
source activate lorals2

SAMPLENAME=$1
KGCODE=$2 # code of samples in the 1000G
TYPE=$3

INPUT=data/$TYPE/04_calc_ase/${SAMPLENAME}_ase_annotated_filtered.tsv
BAM=/gpfs/projects/bsc83/Projects/pantranscriptome/pclavell/03_mapping/data/${TYPE}_general_mapping/genomic/${SAMPLENAME}.bam # genome mapping
TBAM=data/$TYPE/03_transcriptome_mapped/${SAMPLENAME}_preprocessed_Q7.fastq.gz_transcriptome.bam
OUTPUT=data/$TYPE/05_calc_asts/${SAMPLENAME}

mkdir -p data/$TYPE/05_calc_asts
calc_asts \
    -m quant \
    -b $BAM \
    -i $INPUT \
    -o $OUTPUT \
    -x $TBAM

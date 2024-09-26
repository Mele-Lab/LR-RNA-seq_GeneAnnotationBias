#!/bin/bash

SAMPLENAME=$1
KGCODE=$2 # code of samples in the 1000G
INPUT=data/04_calc_ase/${SAMPLENAME}_ase_annotated.tsv
BAM=/gpfs/projects/bsc83/Projects/pantranscriptome/pclavell/03_mapping/data/general_mapping/genomic/${SAMPLENAME}.bam # genome mapping
TBAM=data/03_transcriptome_mapped/${SAMPLENAME}_transcriptome.bam
OUTPUT=data/06_calc_asts/${SAMPLENAME}_asts.tsv

mkdir -p data/05_calc_asts
calc_asts \
    -m quant \
    -b $BAM \
    -i $INPUT \
    -o $OUTPUT \
    -x $TBAM \
    --filter

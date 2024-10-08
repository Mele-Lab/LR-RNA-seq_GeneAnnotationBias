#!/bin/bash

module load anaconda bedtools
source activate base
source activate lorals2

SAMPLENAME=$1
KGCODE=$2 # code of samples in the 1000G

INPUT=data/05_calc_asts/${SAMPLENAME}_asts_quant.tsv
OUTDIR=data/05_calc_asts/${SAMPLENAME}

process_asts \
    -i $INPUT \
    -g /gpfs/projects/bsc83/Data/gene_annotations/gencode/v47/modified/gencode.v47.primary_assembly.annotation.transcripidv_geneidv_match.tsv \
    -o $OUTDIR
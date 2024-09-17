#!/bin/bash

module load anaconda
conda init
source activate base
conda activate /gpfs/projects/bsc83/utils/conda_envs/nanoplot

INPUT=$1
SAMPLENAME=$(echo $1 | sed 's-.*/--' | sed 's/\.bam//')

NanoPlot \
    -t 112 \
    -o data_pseudogene_masked/genomic_nanoplot \
    -p ${SAMPLENAME}_nanoplot \
    --tsv_stats \
    --bam $INPUT
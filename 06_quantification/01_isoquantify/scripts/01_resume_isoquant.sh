#!/bin/bash

module load anaconda minimap2 samtools

conda init
source activate base
conda activate /gpfs/projects/bsc83/utils/conda_envs/isoquant

BAM=$1
SAMPLENAME=$(echo $BAM | sed 's-.*/--' | sed 's/\.bam//')

TYPE=data/240916_merge_with_geneentries

/gpfs/projects/bsc83/utils/conda_envs/isoquant/IsoQuant/isoquant.py \
    --threads 112 \
    --resume \
    --debug \
    --output $TYPE
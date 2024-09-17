#!/bin/bash

module load anaconda minimap2 samtools

conda init
source activate base
conda activate /gpfs/projects/bsc83/utils/conda_envs/isoquant

/gpfs/projects/bsc83/utils/conda_envs/isoquant/IsoQuant/isoquant.py \
    --threads 48 \
    --resume \
    --output 01_quantification_isoquant/isoquant_output
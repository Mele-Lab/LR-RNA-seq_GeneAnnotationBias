#!/bin/bash

module load anaconda
source activate base
conda init bash
conda activate /gpfs/projects/bsc83/utils/conda_envs/nanoplot
/gpfs/projects/bsc83/utils/conda_envs/nanoplot/bin/NanoPlot \
            -t 112 \
            -o testnanoplotdir\
            -p test1 \
            --tsv_stats \
            -f png \
            --ubam /gpfs/projects/bsc83/Projects/pantranscriptome/pclavell/01_basecalling/data/modifications/20240212_HS_16_CH4_GM18772_10percent_subsampling.bam
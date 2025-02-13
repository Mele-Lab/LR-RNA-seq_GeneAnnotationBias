#!/bin/bash

module load anaconda
module load R/4.3.2
conda init
source activate base
conda activate /gpfs/projects/bsc83/utils/conda_envs/nanoplot

FASTQ=$1
SAMPLE=$(echo $FASTQ | sed 's-.*/--')


#python scripts/snakemake/fastq_to_tsv.py --fastq $FASTQ -o nanoplotstats/$SAMPLE

echo $SAMPLE
Rscript scripts/correlate_Q_length.R nanoplotstats/$SAMPLE $SAMPLE
#rm nanoplotstats/$SAMPLE
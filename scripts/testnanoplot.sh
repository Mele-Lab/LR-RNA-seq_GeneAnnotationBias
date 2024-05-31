#!/bin/bash

module load miniconda
source activate base
conda init
conda activate /gpfs/projects/bsc83/utils/conda_envs/fastqfilter
fastq-filter -q 7 data/tenpercentbam/tenpercentbam.fastq.gz -o filtered_tenpercen.fastq.gz
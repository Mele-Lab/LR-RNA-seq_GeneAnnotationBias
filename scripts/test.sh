#!/bin/bash

# module load dorado

# dorado trim \
# 	--threads 80 \
# 	--emit-fastq \
# 	-v \
#        data/tenpercentbam/tenpercentbam.fastq > outputtrim.fastq	

# module load fastqc
# fastqc \
# 	-o testoutdir \
# 	-t 48 \
# 	-a snakemake/ref/adapters.tsv \
# 	-f fastq \
# 	outputtrim.fastq
module load anaconda
source activate base
conda init
conda activate /gpfs/projects/bsc83/utils/conda_envs/porechop

porechop \
	-i data/tenpercentbam/tenpercentbam.fastq.gz \
	-o porechoptest.fastq.gz \
	--threads 112 \
	--verbosity 2
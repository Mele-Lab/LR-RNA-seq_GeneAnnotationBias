#!/bin/bash

# module load dorado

# dorado trim \
# 	--threads 80 \
# 	--emit-fastq \
# 	-v \
#        data/tenpercentbam/tenpercentbam.fastq > outputtrim.fastq	

# module load fastqc
# fastqc \
# 	-o testoutdir4 \
# 	-t 112 \
# 	-a snakemake/ref/adapters.tsv \
# 	-f fastq \
# 	porechoptest.fastq.gz 




# IT RUNS
module load anaconda
source activate base
conda init
conda activate /gpfs/projects/bsc83/utils/conda_envs/porechop

porechop \
	-i data/tenpercentbam/tenpercentbam.fastq.gz \
	-o porechoptest.fastq.gz \
	--threads 112 \
	--verbosity 2	

module load fastqc
fastqc \
	-o testoutdir4 \
	-t 112 \
	-a snakemake/ref/adapters.tsv \
	-f fastq \
	porechoptest.fastq.gz 

# Change this in adapters.py
# ADAPTERS = [Adapter('SQK-LSK114',
#              start_sequence=('HEAD_ADAPTER', 'AATGTACTTCGTTCAGTTACGTATTGCT'),
#              end_sequence=('TAIL_ADAPTER', 'GCAATACGTAACTGAACGAAGT')),
#             Adapter('CapTrapFull',
#              start_sequence=('FULL_FIVE_LINKER', 'TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGNNNNNNNNNNNNNNNNGTGGTATCAACGCAGAGTAC'),
#              end_sequence=('FULL_THREE_LINKER', 'GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGCGATGCTTTTTTTTTTTTTTTT')),
#             Adapter('CapTrapPrimers',
#              start_sequence=('FIVE_LINKER', 'TCGTCGGCAGCGTC'),
#              end_sequence=('THREE_LINKER', 'GTCTCGTGGGCTCGG'))]
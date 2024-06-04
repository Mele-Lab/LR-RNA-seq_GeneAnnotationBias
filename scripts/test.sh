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




# # IT RUNS
module load anaconda
source activate base
conda init
conda activate /gpfs/projects/bsc83/utils/conda_envs/porechop

porechop \
	-i data/tenpercentbam/tenpercentbam.fastq.gz \
	-o porechoptest.fastq.gz \
	--threads 112 \
	--verbosity 3 \
	--min_split_read_size 200 \
	--extra_middle_trim_good_side 2 \
	--extra_middle_trim_bad_side 2



# module load minimap2
# minimap2 \
# 	-ax splice \
# 	--junc-bed /gpfs/projects/bsc83/MN4/bsc83/Data/gene_annotation/gencode/release_44/modified/gencode.v44.chr_patch_hapl_scaff.annotation_transcriptlevel.bed \
# 	--MD \
# 	-t 112 \
# 	-o minimap_output_postporechop_split.bam \
# 	/gpfs/projects/bsc83/Data/gene_annotations/gencode/v44/modified/gencode.v44.transcripts_simpleID.fa \
# 	porechoptest.fastq.gz



module load fastqc
fastqc \
	-o testoutdir_split_newparams \
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
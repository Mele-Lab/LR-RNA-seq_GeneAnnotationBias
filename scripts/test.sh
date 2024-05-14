#!/bin/bash

module load bedtools
bedtools bamtofastq -i /gpfs/projects/bsc83/Projects/pantranscriptome/pclavell/01_basecalling/data/modifications/4kbam_entries.bam  -fq test/convertedto.fastq


# module load anaconda
#         conda init
#         source activate base
#         conda activate /gpfs/projects/bsc83/utils/conda_envs/duplextools_env
#         /gpfs/projects/bsc83/utils/conda_envs/duplextools_env/bin/duplex_tools split_on_adapter \
#             --threads 112 \
#             --adapter_type ONT_sequencing_adapter+CapTrapSeqJoint \
#             /gpfs/projects/bsc83/Projects/pantranscriptome/pclavell/01_basecalling/data/modifications/ \
#             test/test \
#             PCR
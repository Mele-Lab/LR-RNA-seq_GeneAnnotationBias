#!/bin/bash

module load anaconda
conda activate /gpfs/projects/bsc83/utils/conda_envs/duplextools_env
/gpfs/projects/bsc83/utils/conda_envs/duplextools_env/bin/duplex_tools split_on_adapter \
    --threads 112 \
    --adapter_type 'ONT_sequencing_adapter+CapTrapSeqJoint' \
    20240212_HS_16_CH4_GM18772_10percent_subsampling.fastq \
    test \
    PCR
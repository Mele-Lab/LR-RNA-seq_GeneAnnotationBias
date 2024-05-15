#!/bin/bash

module load anaconda
conda init
source activate base
conda activate /gpfs/projects/bsc83/utils/conda_envs/duplextools_env
/gpfs/projects/bsc83/utils/conda_envs/duplextools_env/bin/duplex_tools split_on_adapter \
    --threads 112 \
    --adapter_type 'ONT_sequencing_adapter+CapTrapSeqJoint' \
    testdata \
    test \
    PCR
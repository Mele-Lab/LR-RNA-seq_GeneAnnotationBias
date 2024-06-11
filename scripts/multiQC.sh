#!/bin/bash

module load hdf5 python/3.12.1

#export PATH="/gpfs/projects/bsc83/utils/python_libs/multiqc/:$PATH"

multiqc scripts/data/*/qc/fastqc/
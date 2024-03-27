#!/bin/bash


module load intel mkl impi gcc/7.2.0
module load miniconda3/py39_4.10.3 && source activate snakemake
snakemake --unlock
snakemake \
-s Snakefile \
-j 100 \
--latency-wait 120 \
--cluster "sbatch \
  -q bsc_ls \
  -c {resources.threads}  \
  --time=2:00:00"

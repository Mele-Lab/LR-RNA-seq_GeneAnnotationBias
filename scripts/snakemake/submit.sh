#!/bin/bash


module load miniconda && source activate sqanti3-snakemake
#snakemake --unlock
snakemake \
-s snakemake/Snakefile \
-j 100 \
--latency-wait 120 \
--cluster "sbatch \
  -q bsc_ls \
  -c {resources.threads}  \
  --time=2:00:00"

#!/bin/bash


module load miniconda && source activate sqanti3-snakemake
snakemake \
-s Snakefile \
-j 100 \
--latency-wait 120 \
--cluster "sbatch \
  -q debug \
  -c {resources.threads}  \
  --time=2:00:00" \
  -n

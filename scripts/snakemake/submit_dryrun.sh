#!/bin/bash


#module load intel mkl impi gcc/7.2.0
module load miniconda && source activate sqanti3-snakemake
snakemake \
-s ONT_preprocessing/scripts/snakemake/Snakefile \
-j 100 \
--latency-wait 120 \
--cluster "sbatch \
  -q debug \
  -c {resources.threads}  \
  --time=2:00:00" \
  -n

#!/bin/bash


module load miniconda && source activate sqanti3-snakemake
#snakemake --unlock
snakemake \
-s snakemake/Snakefile \
-j 100 \
--latency-wait 120 \
--cluster "sbatch \
  -q gp_bscls \
  -c 1  \
  -A bsc83 \
  -o slurm_out/%a.out \
  --time=2:00:00"

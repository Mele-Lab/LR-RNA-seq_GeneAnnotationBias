#!/bin/bash


module load miniconda && source activate sqanti3-snakemake
#snakemake --unlock
snakemake \
  -s snakemake/Snakefile \
  -j 1 \
  --latency-wait 120 \
  --cluster "sbatch \
  -q gp_bscls \
  -c 1 \
  -A bsc83 \
  -o slurm_out/%x_%j.out \
  --time=2:00:00"

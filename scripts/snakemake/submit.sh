#!/bin/bash


module load miniconda && source activate sqanti3-snakemake
#snakemake --unlock
snakemake \
  -s $PWD/scripts/snakemake/Snakefile \
  -j 100 \
  --configfile $PWD/scripts/snakemake/config.yml \
  --directory $PWD/scripts \
  --latency-wait 120 \
  --cluster "sbatch \
  -q gp_bscls \
  -c {resources.threads} \
  -A bsc83 \
  -o smk_out/%j_%x.out \
  --time=2:00:00"

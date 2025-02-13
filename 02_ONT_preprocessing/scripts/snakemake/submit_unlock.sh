#!/bin/bash


module load miniconda && source activate sqanti3-snakemake
snakemake --unlock \
  -s $PWD/scripts/snakemake/Snakefile \
  -j 100 \
  --configfile $PWD/scripts/snakemake/config.yml \
  --directory $PWD/scripts \
  --keep-going \
  --latency-wait 120 \
  --rerun-incomplete \
  --cluster "sbatch \
  -q gp_bscls \
  -c {resources.threads} \
  -A bsc83 \
  -o smk_out/{wildcards.sample}/%j_%x.out \
  -t {resources.runtime}"

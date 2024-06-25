#!/bin/bash


# module load miniconda && source activate sqanti3-snakemake
# snakemake \
#   -s snakemake/Snakefile \
#   --dag \
#   -j 1 \
#   --latency-wait 120 \
#   --cluster "sbatch \
#   -q gp_debug \
#   -c {resources.threads}  \
#   --time=2:00:00" \
#   -n \
#   | dot -Tpng > "workflow_dag.png"

module load miniconda && source activate sqanti3-snakemake
#snakemake --unlock
snakemake \
  -s $PWD/scripts/snakemake/Snakefile \
  -j 100 \
  --dag \
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
  -t {resources.runtime}" \
  -n \
  | dot -Tpng > "workflow_dag.png"
#!/bin/bash


module load miniconda && source activate sqanti3-snakemake
snakemake \
  -s snakemake/Snakefile \
  --dag \
  -j 1 \
  --latency-wait 120 \
  --cluster "sbatch \
  -q gp_debug \
  -c {resources.threads}  \
  --time=2:00:00" \
  -n \
  | dot -Tpng > "workflow_dag.png"


#!/bin/bash

#SBATCH --job-name=count
#SBATCH --output=/gpfs/projects/bsc83/Projects/gencode_diversity/deduplication/logs/slurm-%j-count.out
#SBATCH --error=/gpfs/projects/bsc83/Projects/gencode_diversity/deduplication/logs/slurm-%j-count.err
#SBATCH --chdir=
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --tasks-per-node=1
#SBATCH --qos=debug
#SBATCH --time=2:00:00
#SBATCH --constraint=

echo '----------------------'
date
echo '----------------------'
echo ' '

bash 04_dedup/02.1_count_duplicates.sh 04_dedup/01_dedup_results/FAX00275_juncbed_transcriptome.sorted_0ED_percontig.tsv      

echo ' '
echo '----------------------'
date
echo '----------------------'
echo ' '


#!/bin/bash

#SBATCH --job-name=filter
#SBATCH --output=/gpfs/projects/bsc83/Projects/gencode_diversity/deduplication/logs/slurm-%j-filter.out
#SBATCH --error=/gpfs/projects/bsc83/Projects/gencode_diversity/deduplication/logs/slurm-%j-filter.err
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

bash 03_genome_mapping/02_filter_and_transform_results.sh 03_genome_mapping/01_alignment_results/FAX00275_juncbed_transcriptome.sam      

echo ' '
echo '----------------------'
date
echo '----------------------'
echo ' '


#!/bin/bash

#SBATCH --job-name=testdedup
#SBATCH --output=/gpfs/projects/bsc83/Projects/pantranscriptome/Pau/01_deduplication/logs/slurm-%j-testdedup.out
#SBATCH --error=/gpfs/projects/bsc83/Projects/pantranscriptome/Pau/01_deduplication/logs/slurm-%j-testdedup.err
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

bash 05_deduplication/01_deduplication.sh 03_genome_mapping/01_alignment_results/FAX00275_juncbed_transcriptome.sorted.bam 05_deduplication/01_deduplicated_bams      

echo ' '
echo '----------------------'
date
echo '----------------------'
echo ' '


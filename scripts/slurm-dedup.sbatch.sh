#!/bin/bash

#SBATCH --job-name=dedup
#SBATCH --output=/gpfs/projects/bsc83/Projects/gencode_diversity/deduplication/logs/slurm-%j-dedup.out
#SBATCH --error=/gpfs/projects/bsc83/Projects/gencode_diversity/deduplication/logs/slurm-%j-dedup.err
#SBATCH --chdir=
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --tasks-per-node=1
#SBATCH --qos=bsc_ls
#SBATCH --time=10:00:00
#SBATCH --constraint=

echo '----------------------'
date
echo '----------------------'
echo ' '

bash 04_dedup/01_dedup.sh 03_genome_mapping/01_alignment_results/FAX00275_juncbed_transcriptome.sorted.bam 4      

echo ' '
echo '----------------------'
date
echo '----------------------'
echo ' '


#!/bin/bash

#SBATCH --job-name=align
#SBATCH --output=/gpfs/projects/bsc83/Projects/gencode_diversity/deduplication/logs/slurm-%j-align.out
#SBATCH --error=/gpfs/projects/bsc83/Projects/gencode_diversity/deduplication/logs/slurm-%j-align.err
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

bash 03_genome_mapping/01_minimap.sh 02_extract_UMI/01_extracted_UMI/FAX00275_with_extracted_UMI.fasta      

echo ' '
echo '----------------------'
date
echo '----------------------'
echo ' '


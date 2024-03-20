#!/bin/bash

#SBATCH --job-name=extractumi
#SBATCH --output=/gpfs/projects/bsc83/Projects/gencode_diversity/deduplication/logs/slurm-%j-extractumi.out
#SBATCH --error=/gpfs/projects/bsc83/Projects/gencode_diversity/deduplication/logs/slurm-%j-extractumi.err
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

bash 02_extract_UMI/01_extract_UMI.sh 01_blast/01_results/FAX00275.filtered.sam      

echo ' '
echo '----------------------'
date
echo '----------------------'
echo ' '


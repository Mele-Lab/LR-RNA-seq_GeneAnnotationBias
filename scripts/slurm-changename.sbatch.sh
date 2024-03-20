#!/bin/bash

#SBATCH --job-name=changename
#SBATCH --output=/gpfs/projects/bsc83/Projects/gencode_diversity/deduplication/logs/slurm-%j-changename.out
#SBATCH --error=/gpfs/projects/bsc83/Projects/gencode_diversity/deduplication/logs/slurm-%j-changename.err
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

mv 01_blast/01_results/00.xml 01_blast/01_results/FAX00275.xml      

echo ' '
echo '----------------------'
date
echo '----------------------'
echo ' '


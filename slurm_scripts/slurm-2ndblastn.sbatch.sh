#!/bin/bash

#SBATCH --job-name=2ndblastn
#SBATCH --output=/gpfs/projects/bsc83/Projects/gencode_diversity/deduplication/logs/slurm-%j-2ndblastn.out
#SBATCH --error=/gpfs/projects/bsc83/Projects/gencode_diversity/deduplication/logs/slurm-%j-2ndblastn.err
#SBATCH --chdir=
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --tasks-per-node=1
#SBATCH --qos=bsc_ls
#SBATCH --time=5:00:00
#SBATCH --constraint=

echo '----------------------'
date
echo '----------------------'
echo ' '

bash 01_blast/01_blast.sh 00_data/FAX06702_pass_d8defc63_883fb55f_CONCAT.fasta 01_blast/00_reference/fiveLinkDb FAX06702      

echo ' '
echo '----------------------'
date
echo '----------------------'
echo ' '


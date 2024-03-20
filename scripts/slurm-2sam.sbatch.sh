#!/bin/bash

#SBATCH --job-name=2sam
#SBATCH --output=/gpfs/projects/bsc83/Projects/gencode_diversity/deduplication/logs/slurm-%j-2sam.out
#SBATCH --error=/gpfs/projects/bsc83/Projects/gencode_diversity/deduplication/logs/slurm-%j-2sam.err
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

bash 01_blast/02_xml2sam.sh 01_blast/01_results/FAX06702.xml 01_blast/00_reference/fiveLinkDb.fasta 00_data/FAX06702_pass_d8defc63_883fb55f_CONCAT.fasta FAX06702      

echo ' '
echo '----------------------'
date
echo '----------------------'
echo ' '


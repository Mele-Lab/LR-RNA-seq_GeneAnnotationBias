#!/bin/bash

#SBATCH --job-name=blast2
#SBATCH --output=/gpfs/projects/bsc83/Projects/gencode_diversity/deduplication/logs/slurm-%j-blast2.out
#SBATCH --error=/gpfs/projects/bsc83/Projects/gencode_diversity/deduplication/logs/slurm-%j-blast2.err
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

bash 01_blast/01_blast.sh 00_data/FAX00275_pass_0184915a_b0ce8917_CONCAT_head100000.fasta 01_blast/00_reference/fiveLinkDb 01_blast/01_results/top100res.xml       

echo ' '
echo '----------------------'
date
echo '----------------------'
echo ' '


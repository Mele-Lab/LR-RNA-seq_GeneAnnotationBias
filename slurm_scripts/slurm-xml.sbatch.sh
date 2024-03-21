#!/bin/bash

#SBATCH --job-name=xml
#SBATCH --output=/gpfs/projects/bsc83/Projects/gencode_diversity/deduplication/logs/slurm-%j-xml.out
#SBATCH --error=/gpfs/projects/bsc83/Projects/gencode_diversity/deduplication/logs/slurm-%j-xml.err
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

bash 01_blast/02_xml2sam.sh 01_blast/01_results/FAX00275.xml 01_blast/00_reference/fiveLinkDb.fasta 00_data/FAX00275_pass_0184915a_b0ce8917_CONCAT.fasta      

echo ' '
echo '----------------------'
date
echo '----------------------'
echo ' '


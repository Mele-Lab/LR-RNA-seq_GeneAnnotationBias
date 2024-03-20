#!/bin/bash

#SBATCH --job-name=index
#SBATCH --output=/gpfs/projects/bsc83/Projects/gencode_diversity/deduplication/logs/slurm-%j-index.out
#SBATCH --error=/gpfs/projects/bsc83/Projects/gencode_diversity/deduplication/logs/slurm-%j-index.err
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

module load minimap2; minimap2 -d 03_genome_mapping/gencodev44_transcriptome.mmi /gpfs/projects/bsc83/Data/gene_annotation/gencode/release_44/gencode.v44.transcripts_simple_transcript_entry.fa      

echo ' '
echo '----------------------'
date
echo '----------------------'
echo ' '


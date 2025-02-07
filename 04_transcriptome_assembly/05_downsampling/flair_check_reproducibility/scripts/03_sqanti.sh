#!/bin/bash

MYGTF=$1
MYPREFIX=$(echo $1 | sed 's-.*/--' | sed 's/\.gtf//')

GENOME=/gpfs/projects/bsc83/Data/assemblies/GRCh38/GRCh38.primary_assembly.genome.fa
GENCODE=/gpfs/projects/bsc83/Data/gene_annotations/gencode/v47/gencode.v47.primary_assembly.annotation.gtf

module load anaconda
conda init
source activate base
conda activate /gpfs/projects/bsc83/utils/conda_envs/SQANTI3-5.2.1

mkdir data/sqanti
mkdir data/sqanti/$MYPREFIX
python /gpfs/projects/bsc83/utils/conda_envs/SQANTI3-5.2.1/sqanti3_qc.py \
    $MYGTF \
    $GENCODE \
    $GENOME \
    --cpus 112 \
    --dir data/sqanti/$MYPREFIX \
    --output $MYPREFIX \
    --report pdf \
    --isoAnnotLite \
    
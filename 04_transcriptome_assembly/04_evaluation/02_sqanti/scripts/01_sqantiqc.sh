#!/bin/bash

MYGTF=$1
MYPREFIX=$2
GENOME=/gpfs/projects/bsc83/Data/assemblies/GRCh38/modified/GRCh38.primary_assembly.sirvset4.genome.fa 
GENCODE=/gpfs/projects/bsc83/Data/gene_annotations/gencode/v47/gencode.v47.primary_assembly.annotation.gtf
MYGTF=$GENCODE
module load anaconda
conda init
source activate base
conda activate /gpfs/projects/bsc83/utils/conda_envs/SQANTI3-5.2.1
head $MYGTF
mkdir -p data/$MYPREFIX
python /gpfs/projects/bsc83/utils/conda_envs/SQANTI3-5.2.1/sqanti3_qc.py \
    $MYGTF \
    $GENCODE \
    $GENOME \
    --cpus 112 \
    --dir data/$MYPREFIX \
    --output $MYPREFIX
    
#!/bin/bash

module load anaconda minimap2 samtools

conda init
source activate base
conda activate /gpfs/projects/bsc83/utils/conda_envs/isoquant

GENOMEFASTA=/gpfs/projects/bsc83/Data/assemblies/GRCh38/GRCh38.primary_assembly.genome.fa
#ANNOTATIONGTF=/gpfs/projects/bsc83/Data/gene_annotations/gencode/v47/gencode.v47.primary_assembly.annotation.gtf
#ANNOTATIONGTF=/gpfs/projects/bsc83/Projects/pantranscriptome/pclavell/06_quantification/01_isoquantify/ref/gencode.v47.primary_assembly.annotation.db
#ANNOTATIONGTF=/gpfs/projects/bsc83/Projects/pantranscriptome/novelannotations/merged/240909_merged_withoutchrEBV_placeholdergene.gtf
#ANNOTATIONGTF=/gpfs/projects/bsc83/Projects/pantranscriptome/pclavell/06_quantification/01_isoquantify/ref/240909_merged_withoutchrEBV_placeholdergene.db
ANNOTATIONGTF=/gpfs/projects/bsc83/Projects/pantranscriptome/novelannotations/merged/poder_v1.sorted.gtf
# ANNOTATIONGTF=/gpfs/projects/bsc83/Projects/pantranscriptome/novelannotations/merged/240910_uma.gtf


BAM=$1
TYPE=data/poder
#TYPE=data/gencodev47
SAMPLENAME=$(echo $1 | sed 's-.*/--' | sed 's/\.bam//')
mkdir -p $TYPE
# Quantification
/gpfs/projects/bsc83/utils/conda_envs/isoquant/IsoQuant/isoquant.py \
    --data_type ont \
    --reference $GENOMEFASTA \
    --genedb $ANNOTATIONGTF \
    --complete_genedb \
    --bam $BAM \
    --prefix $SAMPLENAME \
    --threads 112 \
    --output $TYPE \
    --no_model_construction \
    --transcript_quantification unique_only \
    --gene_quantification unique_splicing_consistent

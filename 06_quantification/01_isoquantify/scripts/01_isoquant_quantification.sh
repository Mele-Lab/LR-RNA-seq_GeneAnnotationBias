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
ANNOTATIONGTF=/gpfs/projects/bsc83/Projects/pantranscriptome/novelannotations/merged/240917_merge_associatedgene2isoform_noambigousISM_FSM_genic_add_gene_entries.gtf


BAM=$1
TYPE=data/240917_merge_with_geneentries_nocomplete
#TYPE=data/gencodev47
SAMPLENAME=$(echo $1 | sed 's-.*/--' | sed 's/\.bam//')
mkdir -p $TYPE
# Discovery and quantification
/gpfs/projects/bsc83/utils/conda_envs/isoquant/IsoQuant/isoquant.py \
    --data_type ont \
    --reference $GENOMEFASTA \
    --genedb $ANNOTATIONGTF \
    --bam $BAM \
    --prefix $SAMPLENAME \
    --threads 112 \
    --output $TYPE

#!/bin/bash

module load anaconda minimap2 samtools

conda init
source activate base
conda activate /gpfs/projects/bsc83/utils/conda_envs/isoquant

GENOMEFASTA=/gpfs/projects/bsc83/Data/assemblies/GRCh38/GRCh38.primary_assembly.genome.fa
#ANNOTATIONGTF=/gpfs/projects/bsc83/Data/gene_annotations/gencode/v47/gencode.v47.primary_assembly.annotation.gtf
ANNOTATIONGTF=/gpfs/projects/bsc83/Projects/pantranscriptome/pclavell/04_transcriptome_assembly/02_isoquant/ref/gencode.v47.primary_assembly.annotation.db

BAM=$1
SAMPLENAME=$(echo $1 | sed 's-.*/--' | sed 's/\.bam//')

# Discovery and quantification
/gpfs/projects/bsc83/utils/conda_envs/isoquant/IsoQuant/isoquant.py \
    --data_type ont \
    --reference $GENOMEFASTA \
    --genedb $ANNOTATIONGTF \
    --complete_genedb \
    --bam $BAM \
    --report_novel_unspliced true \
    --prefix $SAMPLENAME \
    --threads 112 \
    --sqanti_output \
    --model_construction_strategy sensitive_ont \
    --splice_correction_strategy assembly \
    --delta 8 \
    --output data

# for unguided remove genedb

# --genedb $ANNOTATIONGTF \


# ########### TAMARA'S SETTINGS
# isoquant.py \
#     --fl_data \
#     --check_canonical \
#     --report_novel_unspliced true \
#     -o isoquant_output_ont_custom \
#     --delta 8 \

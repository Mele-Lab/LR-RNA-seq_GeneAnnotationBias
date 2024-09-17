#!/bin/bash

module load anaconda
source activate base
conda init
conda activate /gpfs/projects/bsc83/utils/conda_envs/flair/

GENOMEFASTA=/gpfs/projects/bsc83/Data/assemblies/GRCh38/GRCh38.primary_assembly.genome.fa
ANNOTATIONGTF=/gpfs/projects/bsc83/Data/gene_annotations/gencode/v47/gencode.v47.primary_assembly.annotation.gtf

BAM=$1
SAMPLENAME=$(echo $1 | sed 's-.*/--' | sed 's/\.bam//')

# convert alignment bam to bed so flair can correct the splicing junctions
mkdir data/00_beds
/gpfs/projects/bsc83/utils/conda_envs/flair/bin/bam2Bed12 -i $BAM > data/00_beds/$SAMPLENAME.bed12

# FLAIR CORRECT
echo "FLAIR CORRECT"

mkdir data/01_correct
/gpfs/projects/bsc83/utils/conda_envs/flair/bin/flair correct \
    -q data/00_beds/$SAMPLENAME.bed12 \
    -f $ANNOTATIONGTF \
    -g $GENOMEFASTA \
    --ss_window 8 \
    --threads 112 \
    --output data/01_correct/$SAMPLENAME

# FLAIR COLLAPSE
echo "FLAIR COLLAPSE"

mkdir data/02_collapse
/gpfs/projects/bsc83/utils/conda_envs/flair/bin/flair collapse \
    -g $GENOMEFASTA \
    --gtf $ANNOTATIONGTF \
    -q data/01_correct/${SAMPLENAME}_all_corrected.bed \
    -r /gpfs/projects/bsc83/Projects/pantranscriptome/pclavell/02_ONT_preprocessing/data/q10/${SAMPLENAME}_preprocessed_Q10.fastq.gz \
    --stringent \
    --check_splice \
    --generate_map \
    --annotation_reliant generate \
    --threads 112 \
    --output data/02_collapse/$SAMPLENAME



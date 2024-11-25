#!/bin/bash

GENOMEFASTA=/gpfs/projects/bsc83/Data/assemblies/GRCh38/GRCh38.primary_assembly.genome.fa
ANNOTATIONGTF=/gpfs/projects/bsc83/Data/gene_annotations/gencode/v47/gencode.v47.primary_assembly.annotation.gtf

# The following file must be in the same directory so just filename not path
FILEPATH=$1
SAMPLENAME=$2
mkdir -p data/espresso_s/
mkdir -p data/espresso_s/$SAMPLENAME
echo -e "$FILEPATH\t$SAMPLENAME" > data/espresso_s/$SAMPLENAME/samples.tsv

module load anaconda
source activate espresso


perl /apps/GPP/ANACONDA/2023.07/envs/espresso/bin/ESPRESSO_S.pl \
    -L data/espresso_s/$SAMPLENAME/samples.tsv \
    -F $GENOMEFASTA  \
    -A $ANNOTATIONGTF \
    -O data/espresso_s/$SAMPLENAME \
    -T 112

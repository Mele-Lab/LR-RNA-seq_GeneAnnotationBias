#!/bin/bash

GENOMEFASTA=/gpfs/projects/bsc83/Data/assemblies/GRCh38/GRCh38.primary_assembly.genome.fa
ANNOTATIONGTF=/gpfs/projects/bsc83/Data/gene_annotations/refseq/v110/modified/GCF_000001405.40_GRCh38.p14_genomic.UCSC_contignames.gtf

# The following file must be in the same directory so just filename not path
SAMPLENAME=$2

module load anaconda
source activate espresso

mkdir -p data/espresso_q
mkdir -p data/espresso_q/$SAMPLENAME
echo "0---THIRD STEP---0"
python /gpfs/projects/bsc83/utils/conda_envs/espresso/snakemake/scripts/combine_espresso_c_output_for_q.py \
    --c-base-work-dir data/espresso_c/$SAMPLENAME \
    --new-base-dir data/espresso_q/$SAMPLENAME

perl /apps/GPP/ANACONDA/2023.07/envs/espresso/bin/ESPRESSO_Q.pl \
    -L data/espresso_q/$SAMPLENAME/samples.tsv.updated \
    -A $ANNOTATIONGTF \
    -T 112

#!/bin/bash

module load anaconda bedtools bedops
source activate base
source activate lorals2

SAMPLENAME=$1
KGCODE=$2 # code of samples in the 1000G
INPUT=data/04_calc_ase/${SAMPLENAME}_ase.tsv # output from calc_ase
GTF=/gpfs/projects/bsc83/Data/gene_annotations/gencode/v47/gencode.v47.primary_assembly.annotation.gtf
VCF=data/01_samples_vcf/${KGCODE}_het.vcf.gz
OUTPUT=data/04_calc_ase/${SAMPLENAME}_ase_annotated.tsv
ANNOT=$(echo $(basename $GTF) | sed 's/\.gtf//')

#gtf2bed < $GTF > ref/${ANNOT}.bed


annotate_ase \
    -i $INPUT \
    -b ref/${ANNOT}.bed \
    -f $VCF \
    -o $OUTPUT \
    --blacklist ref/hg38-blacklist.v2.bed \
    --mapping ref/wgEncodeCrgMapabilityAlign100mer.hg38.mappingover5.bed.gz
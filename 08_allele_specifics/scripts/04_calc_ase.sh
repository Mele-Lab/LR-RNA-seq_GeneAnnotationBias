#!/bin/bash

module load anaconda bedtools
source activate base
source activate lorals2
#used input: 6_NI1_GM18486 NA18486

SAMPLENAME=$1
KGCODE=$2 # code of samples in the 1000G
TYPE=$3
#BAM=/gpfs/projects/bsc83/Projects/pantranscriptome/pclavell/03_mapping/data/general_mapping/genomic/${SAMPLENAME}.bam # genome mapping
BAM=/gpfs/projects/bsc83/Projects/pantranscriptome/pclavell/08_allele_specifics/data/$TYPE/02_haplotype_alignments/${SAMPLENAME}_preprocessed_Q7_reads_aln_sorted.merged.bam
VCF=data/01_samples_vcf/${KGCODE}_het.vcf.gz
OUTPUT=data/$TYPE/04_calc_ase/${SAMPLENAME}_ase.tsv

mkdir -p data/$TYPE/04_calc_ase

which calc_ase

calc_ase \
    -b $BAM \
    -f $VCF \
    -o $OUTPUT




#!/bin/bash

FASTQ=/gpfs/projects/bsc83/Projects/pantranscriptome/pclavell/02_ONT_preprocessing/data/q7/${1}_preprocessed_Q7.fastq.gz
KGCODE=$2

# Check that this is the GTF of interest
JUNCBED=/gpfs/projects/bsc83/Data/gene_annotations/gencode/v47/modified/gencode.v47.primary_assembly.sirvset4.annotation.junc.bed


HAPLOTYPES=/gpfs/projects/bsc83/Projects/pantranscriptome/pclavell/08_allele_specifics/data/01_samples_vcf/$KGCODE
OUTDIR=/gpfs/projects/bsc83/Projects/pantranscriptome/pclavell/08_allele_specifics/data/02_haplotype_alignments

echo $FASTQ
echo $HAPLOTYPES

mkdir -p $OUTDIR
bash scripts/02sub_hap_aligner.sh --fastq $FASTQ \
    --reference $HAPLOTYPES \
    --threads 112 \
    --outdir $OUTDIR \
    --juncbed $JUNCBED
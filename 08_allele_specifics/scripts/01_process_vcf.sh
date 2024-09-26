#!/bin/bash

VCF=/gpfs/projects/bsc83/Data/1000G/x30_vcf/original_data/1kGP_high_coverage_Illumina.all_chr.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
SAMPLES=../00_metadata/array_samples_1000G30x

GENOME=/gpfs/projects/bsc83/Data/assemblies/GRCh38/GRCh38.primary_assembly.genome.fa
OUTDIR=data/01_samples_vcf



bash scripts/01sub_process_vcf.sh \
    --vcf $VCF \
    --fasta $GENOME \
    --outdir $OUTDIR \
    --samples $SAMPLES
#!/bin/bash

module load bcftools

bcftools concat \
    --threads 112 \
    -O z \
    -o /gpfs/projects/bsc83/Data/1000G/x30_vcf/original_data/1kGP_high_coverage_Illumina.all_chr.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
    /gpfs/projects/bsc83/Data/1000G/x30_vcf/original_data/1kGP_high_coverage_Illumina.chr*.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
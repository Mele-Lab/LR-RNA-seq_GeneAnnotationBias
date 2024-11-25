#!/bin/bash

zcat /gpfs/projects/bsc83/Data/1000G/x30_vcf/original_data/1kGP_high_coverage_Illumina.chr6.filtered.SNV_INDEL_SV_phased_panel.vcf.gz | awk '$2 >= 33070173 && $2 <= 33070317 {print $0}' > data/example_HLA_variants.vcf 
cat header_1000G_vcf.tsv data/example_HLA_variants.vcf > data/example_HLA_variants_v2.vcf 



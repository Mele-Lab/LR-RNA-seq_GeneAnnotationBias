#!/bin/bash

GENCODE=/gpfs/projects/bsc83/Data/gene_annotations/gencode/v47/gencode.v47.primary_assembly.annotation.gtf
PODER=/gpfs/projects/bsc83/Projects/pantranscriptome/novelannotations/merged/poder_v1.gtf
VARIANTS=/gpfs/projects/bsc83/Data/1000G/x30_vcf/original_data/1kGP_high_coverage_Illumina.all_chr.filtered.SNV_INDEL_SV_phased_panel.vcf.gz

awk -F"\t" '$3=="gene" {print $1,$4,$5}' $GENCODE > data/genes
awk -F"\t" '$3=="gene" {print $1,$4,$5}' $PODER >> data/genes

sort data/genes | uniq | awk 'BEGIN{FS=" ";OFS="\t"}{$1=$1;print $0}'> data/unique_genes.tsv

module load bedtools

#zcat $VARIANTS | cut -f1,2,3,4,5,8 > data/variants.vcf
bedtools intersect -a $VARIANTS -b data/unique_genes.tsv | gzip > 1000G_variants_in_genes.vcf.gz
#!/bin/bash

module load plink

cat scripts/data/merged/final/samples_1000G_chr_merged_biallelic.vcf |\
 grep "#" > finalparsed.vcf

cat 05_variant_calling/scripts/data/merged/final/samples_1000G_chr_merged_biallelic.vcf |\
 grep -v "#" |\
 awk 'BEGIN{OFS="\t"}{for(i=10; i<=NF; i++){sub(/:.*/, "", $i)}; for(i=10; i<=NF; i++){gsub(/\//, "|", $i)}; print}' >> finalparsed.vcf


plink2 \
    --vcf finalparsed.vcf \
    --allow-extra-chr \
    --set-missing-var-ids @:# \
    --make-bed \
    --pca \
    --threads 112\
    --out data/pca
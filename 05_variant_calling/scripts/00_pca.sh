#!/bin/bash

module load plink
#scripts/data/merged/final/samples_1000G_chr_merged_biallelic_parsed.vcf
zcat scripts/data/old/merged/fortyfiveincompletesamplestestmerge.vcf.gz |\
 grep "#" > data/finalparsed.vcf

zcat scripts/data/old/merged/fortyfiveincompletesamplestestmerge.vcf.gz |\
 grep -v "#" |\
 awk 'BEGIN{OFS="\t"}{for(i=10; i<=NF; i++){sub(/:.*/, "", $i)}; for(i=10; i<=NF; i++){gsub(/\//, "|", $i)}; print}' >> data/finalparsed.vcf


# plink2 \
#     --vcf finalparsed.vcf \
#     --allow-extra-chr \
#     --set-missing-var-ids @:# \
#     --make-bed \
#     --pca \
#     --threads 112\
#     --out data/pca
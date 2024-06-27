#!/bin/bash

module load plink

# prune and create pca
plink \
    --vcf $VCF \
    --allow-extra-chr \
    --set-missing-var-ids @:# \
    --extract pantranscriptome_pruned.prune.in \
    --make-bed \
    --pca \
    --out pantranscriptome_pruned
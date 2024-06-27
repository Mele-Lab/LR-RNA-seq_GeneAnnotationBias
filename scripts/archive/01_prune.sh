#!/bin/bash

module load plink

# perform linkage pruning - i.e. identify prune sites
plink \
    --vcf $VCF \
    --allow-extra-chr \
    --set-missing-var-ids @:# \
    --indep-pairwise 50 10 0.1 \
    --out pantranscriptome_pruned


    #--double-id
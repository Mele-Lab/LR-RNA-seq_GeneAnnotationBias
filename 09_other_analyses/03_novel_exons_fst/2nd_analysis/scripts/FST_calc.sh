#!/usr/bin/sh

## Keep only biallelic
vcftools \
--gzvcf /nfs/users/rg/fdegalez/data/variants/1000G_variants_in_genes.withHeader.vcf.gz \
--recode \
--recode-INFO-all \
--min-alleles 2 \
--max-alleles 2 \
--remove-indels \
--stdout |\
gzip > /nfs/users/rg/fdegalez/data/variants/1000G_variants_in_genes.withHeader.onlyBiall.vcf.gz


## Extract per population

for POP in PEL CEU YRI LWK ITU
do
echo $POP
vcftools \
    --gzvcf /nfs/users/rg/fdegalez/data/variants/1000G_variants_in_genes.withHeader.10000first.onlyBiall.vcf.gz \
    --keep /nfs/users/rg/fdegalez/data/variants/metadata/${POP}_samples.txt \
    --recode \
    --recode-INFO-all \
    --stdout |\
    gzip > only_${POP}.biall.vcf.gz
done


## FST between populations

CEU="/nfs/users/rg/fdegalez/data/variants/metadata/CEU_samples.txt"
ITU="/nfs/users/rg/fdegalez/data/variants/metadata/ITU_samples.txt"

vcftools \
--gzvcf /nfs/users/rg/fdegalez/data/variants/1000G_variants_in_genes.withHeader.10000first.onlyBiall.vcf.gz \
--weir-fst-pop ${CEU} \
--weir-fst-pop ${ITU} \
--stdout |\
sort | uniq > CEU_vs_ITU.txt



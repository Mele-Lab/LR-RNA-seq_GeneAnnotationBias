#!/bin/bash
#$ -q short-centos79,long-centos79,rg-el7
#$ -pe smp 8
#$ -l virtual_free=20G
#$ -l h_rt=12:00:00
#$ -cwd
#$ -N FST
#$ -j y
#$ -t 1-10:1
#$ -o logs/FST-$JOB_ID-$TASK_ID.out.log
#$ -e logs/FST-$JOB_ID$TASK_ID.err.log


#################
# start message #
#################
start_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] starting on $(hostname)

# make bash behave more robustly
set -e
set -u
set -o pipefail

###################
# set environment #
###################
# Loading module

###############################################
# Submit array according to a list/file       #
###############################################
file="/nfs/users/rg/fdegalez/data/variants/FST/comb.list"
line=`sed "$((SGE_TASK_ID))q;d" ${file}`

###############
# run command #
###############
## Variables (fixed)
echo "Start..."

conda activate vcftools

#### Functions

VCF="/nfs/users/rg/fdegalez/data/variants/1000G_variants_in_genes.withHeader.vcf.gz"
#VCF="/nfs/users/rg/fdegalez/data/variants/1000G_variants_in_genes.withHeader.onlyBiall.vcf.gz"
list_pop_1=$(echo $line | awk '{print $1}')
list_pop_2=$(echo $line | awk '{print $2}')


vcftools \
--gzvcf ${VCF} \
--weir-fst-pop /nfs/users/rg/fdegalez/data/variants/metadata/${list_pop_1}_samples.txt \
--weir-fst-pop /nfs/users/rg/fdegalez/data/variants/metadata/${list_pop_2}_samples.txt \
--stdout |\
sort | uniq > ./${list_pop_1}_${list_pop_2}.FST.txt


###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds

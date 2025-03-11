#!/bin/bash
#$ -q short-centos79,long-centos79,rg-el7
#$ -pe smp 8
#$ -l virtual_free=1G
#$ -l h_rt=00:05:00
#$ -cwd
#$ -N FeatureUnique
#$ -j y
#$ -t 1-22:1
#$ -o logs/FeatureUniqueRefTable-$JOB_ID-$TASK_ID.out.log
#$ -e logs/FeatureUniqueRefTable-$JOB_ID$TASK_ID.err.log


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
file="[PATH]/mapping_minimap/RNAseqLR_mapping/featureIsoUniqueRegion.list"
line="[PATH]/mapping_minimap/RNAseqLR_mapping/HG002.fq_HG002.fa/isoquant/OUT"
line=`sed "$((SGE_TASK_ID))q;d" ${file}`

FASTA_NAME=$(echo $line | cut -d "/" -f 8 | cut -d "_" -f 2 | cut -d "." -f 1)
echo $FASTA_NAME

###############
# run command #
###############
## Variables (fixed)

module load BEDTools/2.30.0-GCC-11.3.0
cd $line
#### Functions
rm -f OUT.countFeatures_specRegions.txt


for feature in exon transcript gene
do
echo $feature

for ref in hg38 T2T
do

echo $ref 

inter=$(cat OUT.${feature}_models_intersectUniqueRegions_${ref}.formatted.bed | wc -l )
included=$(cat OUT.${feature}_models_includedUniqueRegions_${ref}.formatted.bed | wc -l )

echo $feature $ref $inter "inter" >> OUT.countFeatures_specRegions.txt
echo $feature $ref $included "included" >> OUT.countFeatures_specRegions.txt
done
done
###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds

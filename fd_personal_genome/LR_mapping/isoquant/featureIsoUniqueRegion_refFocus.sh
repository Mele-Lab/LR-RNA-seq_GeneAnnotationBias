#!/bin/bash
#$ -q short-centos79,long-centos79,rg-el7
#$ -pe smp 8
#$ -l virtual_free=5G
#$ -l h_rt=00:15:00
#$ -cwd
#$ -N FeatureUnique
#$ -j y
#$ -t 1-12:1
#$ -o logs/FeatureUniqueRefFocus-$JOB_ID-$TASK_ID.out.log
#$ -e logs/FeatureUniqueRefFocus-$JOB_ID$TASK_ID.err.log


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
line="[PATH]/mapping_minimap/RNAseqLR_mapping/HG002.fq_hg38.fa/isoquant/OUT"
line=`sed "$((SGE_TASK_ID))q;d" ${file}`

FASTA_NAME=$(echo $line | cut -d "/" -f 8 | cut -d "_" -f 2 | cut -d "." -f 1)
echo $FASTA_NAME
ref=$(echo $line | cut -d "/" -f 8 | cut -d "_" -f 1 | cut -d "." -f 1)
echo $ref

###############
# run command #
###############
## Variables (fixed)

module load BEDTools/2.30.0-GCC-11.3.0
cd $line
#### Functions

GTF="OUT.transcript_models.gtf"

grep -v "^#" $GTF | cut -f3 | sort | uniq -c > OUT.transcript_models.countFeatures.txt

for feature in exon transcript gene
do
echo $feature

###hg38
unique_region_ref="[PATH]/mapping_minimap/genomic_comparison/${FASTA_NAME}_${ref}/unique_regions_sup1000.sorted.bed"

## intersect
bedtools intersect \
-a $unique_region_ref \
-b OUT.${feature}_models.bed \
-wa \
-wb > OUT.${feature}_models_intersectUniqueRegions_${ref}.bed

cat OUT.${feature}_models_intersectUniqueRegions_${ref}.bed | awk -v OFS="\t" '{
print $5 ,$6,$7,$8 }' | sort | uniq > OUT.${feature}_models_intersectUniqueRegions_${ref}.formatted.bed
## included
bedtools intersect \
-a $unique_region_ref \
-b OUT.${feature}_models.bed \
-wa \
-wb \
-F 0.99 > OUT.${feature}_models_includedUniqueRegions_${ref}.bed
cat OUT.${feature}_models_includedUniqueRegions_${ref}.bed | awk -v OFS="\t" '{
print $5 ,$6,$7,$8 }' | sort | uniq > OUT.${feature}_models_includedUniqueRegions_${ref}.formatted.bed

done
###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds

#!/bin/bash

##################
# slurm settings #
##################

# For array
#SBATCH --array=1-6

# where to put stdout / stderr
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# time limit in minutes
#SBATCH --time=00:40:00

# queue
#SBATCH --qos=shorter

# memory (MB)
#SBATCH --mem=20G

# job name
#SBATCH --job-name ReadsFollowed

# cpu slots
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

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


###############################################
# Submit array according to a list/file       #
# slurm arrays are 0 based, sed line counting #
# starts from 1                               #
###############################################
genome=`sed "$((SLURM_ARRAY_TASK_ID))q;d" [PATH]/mapping_minimap/readFollowing/listOfGenome.list`



###############
# run command #
###############
## Links to file

## For extraction reads
reads_in_hg38="[PATH]/mapping_minimap/RNAseqLR_mapping/${genome}.fq_${genome}.fa/isoquant/OUT/OUT.transcript_models_intersectUniqueRegions_hg38.readsOnly.tsv"
reads_in_T2T="[PATH]/mapping_minimap/RNAseqLR_mapping/${genome}.fq_${genome}.fa/isoquant/OUT/OUT.transcript_models_intersectUniqueRegions_T2T.readsOnly.tsv"

## For reads info
mapping_PG="[PATH]/mapping_minimap/RNAseqLR_mapping/${genome}.fq_${genome}.fa/alignment.primary.sorted.bed"
mapping_hg38="[PATH]/mapping_minimap/RNAseqLR_mapping/${genome}.fq_hg38.fa/alignment.primary.sorted.bed"
mapping_T2T="[PATH]/mapping_minimap/RNAseqLR_mapping/${genome}.fq_T2T.fa/alignment.primary.sorted.bed"

## For Tr info
tr_PG="[PATH]/mapping_minimap/RNAseqLR_mapping/${genome}.fq_${genome}.fa/isoquant/OUT/OUT.transcript_model_reads.tsv.gz"
tr_hg38="[PATH]/mapping_minimap/RNAseqLR_mapping/${genome}.fq_hg38.fa/isoquant/OUT/OUT.transcript_model_reads.tsv.gz"
tr_T2T="[PATH]/mapping_minimap/RNAseqLR_mapping/${genome}.fq_T2T.fa/isoquant/OUT/OUT.transcript_model_reads.tsv.gz"


##### Personnal genome
# Extrat exon and intron lengths (cum) from reads mapping to PG
cat $mapping_PG | awk '{
  split($11, arr, ",")
  sum_ex=0
  for (i in arr) sum_ex+=arr[i]
  split($12, arr, ",")
  sum_in=0
  for (i in arr) sum_in+=arr[i]
  print $4, $10, sum_ex, sum_in
}' >  readsMappingIn_${genome}_${genome}_info.txt
## Reads from regions of interest
grep -f $reads_in_hg38 readsMappingIn_${genome}_${genome}_info.txt > readsMappingIn_${genome}_${genome}_inSpecRegion_hg38_info.txt
grep -f $reads_in_T2T readsMappingIn_${genome}_${genome}_info.txt > readsMappingIn_${genome}_${genome}_inSpecRegion_T2T_info.txt

zcat ${tr_PG} | grep -f $reads_in_hg38 > transcript_model_reads_${genome}_${genome}_intersectUniqueRegionsOfhg38.tsv
zcat ${tr_PG} | grep -f $reads_in_T2T > transcript_model_reads_${genome}_${genome}_intersectUniqueRegionsOfT2T.tsv

######

###### Reference : hg38
# Extrat exon and intron lengths (cum) from reads mapping to ref
cat $mapping_hg38 | awk '{
  split($11, arr, ",")
  sum_ex=0
  for (i in arr) sum_ex+=arr[i]
  split($12, arr, ",")
  sum_in=0
  for (i in arr) sum_in+=arr[i]
  print $4, $10, sum_ex, sum_in
}' > readsMappingIn_${genome}_hg38_info.txt
## Reads from regions of interest
grep -f $reads_in_hg38 readsMappingIn_${genome}_hg38_info.txt > readsMappingIn_${genome}_inSpecRegion_hg38_info.txt
zcat ${tr_hg38} | grep -f $reads_in_hg38 > transcript_model_reads_${genome}_hg38_intersectUniqueRegionsOfhg38.tsv
######


###### Refrenece : T2T
# Extrat exon and intron lengths (cum) from reads mapping to ref
cat $mapping_T2T | awk '{
  split($11, arr, ",")
  sum_ex=0
  for (i in arr) sum_ex+=arr[i]
  split($12, arr, ",")
  sum_in=0
  for (i in arr) sum_in+=arr[i]
  print $4, $10, sum_ex, sum_in
}' > readsMappingIn_${genome}_T2T_info.txt
## Reads from regions of interest
grep -f $reads_in_T2T readsMappingIn_${genome}_T2T_info.txt > readsMappingIn_${genome}_inSpecRegion_T2T_info.txt
zcat ${tr_T2T} | grep -f $reads_in_T2T > transcript_model_reads_${genome}_T2T_intersectUniqueRegionsOfT2T.tsv
######


###############
# end message #
###############
cgroup_dir=$(awk -F: '{print $NF}' /proc/self/cgroup)
peak_mem=`cat /sys/fs/cgroup$cgroup_dir/memory.peak`
echo [$(date +"%Y-%m-%d %H:%M:%S")] peak memory is $(echo $peak_mem | numfmt --to=iec) \($peak_mem\) bytes
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds

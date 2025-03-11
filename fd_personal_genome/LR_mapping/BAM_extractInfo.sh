#!/bin/bash
#$ -q short-centos79,long-centos79,rg-el7
#$ -pe smp 8
#$ -l virtual_free=2G
#$ -l h_rt=12:00:00
#$ -cwd
#$ -N BAM_extractInfo
#$ -j y
#$ -t 1-22:1
#$ -o logs/BAM_extractInfo-$JOB_ID-$TASK_ID.out.log
#$ -e logs/BAM_extractInfo-$JOB_ID$TASK_ID.err.log


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
file="[PATH]/mapping_minimap/RNAseqLR_mapping/BAM_extractInfo.queries.tsv"
line=`sed "$((SGE_TASK_ID))q;d" ${file}`
bam=$(echo $line | awk '{print $1}')
output_path=$(dirname $bam)

###############
# run command #
###############
## Variables (fixed)
echo $bam

#### Functions

# Align the genome
samtools view $bam | awk '
  {
    nm = ms = as = tp = s1 = s2 = de = rl = "NA"
    for (i=12; i<=NF; i++) {
      if ($i ~ /^NM:i:/) nm = substr($i, 6)
      else if ($i ~ /^ms:i:/) ms = substr($i, 6)
      else if ($i ~ /^AS:i:/) as = substr($i, 6)
      else if ($i ~ /^tp:A:/) tp = substr($i, 6)
      else if ($i ~ /^s1:i:/) s1 = substr($i, 6)
      else if ($i ~ /^s2:i:/) s2 = substr($i, 6)
      else if ($i ~ /^de:f:/) de = substr($i, 6)
    }
    print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" nm "\t" ms "\t" as "\t" tp "\t" s1 "\t" s2 "\t" de
  }
' > $output_path/alignment.sorted.custom.bed

###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds

#!/bin/bash
#$ -q short-centos79,long-centos79,rg-el7
#$ -pe smp 8
#$ -l virtual_free=10G
#$ -l h_rt=00:30:00
#$ -cwd
#$ -N extractFlag
#$ -j y
#$ -t 1-50:1
#$ -o logs/extractFlag-$JOB_ID-$TASK_ID.out.log
#$ -e logs/extractFlag-$JOB_ID$TASK_ID.err.log


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
file="[PATH]/mapping_minimap/RNAseqLR_mapping/extractFlag.tsv"
line=`sed "$((SGE_TASK_ID))q;d" ${file}`
file=$(echo $line | awk '{print $1}')


###############
# run command #
###############
## Variables (fixed)
echo $file

#### Functions
# extract the path of the file
file_path=$(dirname $file)
file_name=$(basename $file | sed "s/\.bam//g")

samtools view $file | cut -f2 | sort | uniq -c > ${file_path}/${file_name}.extractFlag.txt

###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds

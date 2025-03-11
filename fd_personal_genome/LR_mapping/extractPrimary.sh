#!/bin/bash
#$ -q short-centos79,long-centos79,rg-el7
#$ -pe smp 8
#$ -l virtual_free=15G
#$ -l h_rt=00:15:00
#$ -cwd
#$ -N exPrim
#$ -j y
#$ -t 1-22:1
#$ -o logs/exPrim-$JOB_ID-$TASK_ID.out.log
#$ -e logs/exPrim-$JOB_ID-$TASK_ID.err.log


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
module load BEDTools/2.30.0-GCC-11.3.0

###############################################
# Submit array according to a list/file       #
###############################################
file="[PATH]/mapping_minimap/RNAseqLR_mapping/extractPrimary.list"
line=`sed "$((SGE_TASK_ID))q;d" ${file}`

###############
# run command #
###############
## Variables (fixed)

PATH_ToWork=$(echo $line)
cd $PATH_ToWork

genome_name=$(echo "$PATH_ToWork" | awk -F'_' '{print $(NF)}' | awk -F'.' '{print $1}')


BAM="alignment.primary.sorted.bam"
bedtools bamtobed -i alignment.primary.sorted.bam -bed12 > alignment.primary.sorted.bed


###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds

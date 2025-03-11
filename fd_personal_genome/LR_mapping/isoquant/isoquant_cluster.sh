#!/bin/bash
#$ -q short-centos79,long-centos79,rg-el7
#$ -pe smp 8
#$ -l virtual_free=40G
#$ -l h_rt=12:00:00
#$ -cwd
#$ -N Isoquant
#$ -j y
#$ -t 1-22:1
#$ -o logs/Isoquant-$JOB_ID-$TASK_ID.out.log
#$ -e logs/Isoquant-$JOB_ID$TASK_ID.err.log


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
file="[PATH]/mapping_minimap/RNAseqLR_mapping/isoquant_list.txt"
line=`sed "$((SGE_TASK_ID))q;d" ${file}`

###############
# run command #
###############
## Variables (fixed)

PATH_ToWork=$(echo $line)
cd $PATH_ToWork

genome_name=$(echo "$PATH_ToWork" | awk -F'_' '{print $(NF)}' | awk -F'.' '{print $1}')


conda activate isoquant

BAM="alignment.sorted.bam"
FASTA="[PATH]/data/genome/${genome_name}/${genome_name}.fa"

samtools view -H ${BAM} > alignment.primary.bam
samtools view -F 256 -F4 -F 2048 ${BAM} >> alignment.primary.bam
cat alignment.primary.bam | samtools sort  --threads 16  >  alignment.primary.sorted.bam
samtools index alignment.primary.sorted.bam
rm alignment.primary.bam

isoquant.py \
--reference ${FASTA} \
--bam alignment.primary.sorted.bam \
--data_type "ont" \
-o "isoquant"


###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds

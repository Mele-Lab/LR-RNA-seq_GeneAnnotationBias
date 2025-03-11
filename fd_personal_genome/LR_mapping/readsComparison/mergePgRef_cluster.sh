#!/bin/bash

##################
# slurm settings #
##################

# For array
#SBATCH --array=1-12

# where to put stdout / stderr
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# time limit in minutes
#SBATCH --time=02:30:00

# queue
#SBATCH --qos=shorter

# memory (MB)
#SBATCH --mem=80G

# job name
#SBATCH --job-name mergePgRef

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

line=`sed "1q;d" mergePgRef.list`

line=`sed "$((SLURM_ARRAY_TASK_ID))q;d" mergePgRef.list`

genome=$(echo $line | cut -f1 -d" ")
ref=$(echo $line | cut -f2 -d" ")

echo $genome
echo $ref

###############
# run command #
###############


########### Merging bam/bed and reads/tr info (5)
echo "Merging (3)"
python3 [PATH][PATH]/mapping_minimap/readFollowing/mergeReads_step3_fixed.py \
  [PATH]/mapping_minimap/RNAseqLR_mapping/${genome}.fq_${genome}.fa/alignment.primary.sorted.allInfo.tsv \
  [PATH]/mapping_minimap/RNAseqLR_mapping/${genome}.fq_${ref}.fa/alignment.primary.sorted.allInfo.tsv \
  ./mergePgRef/${genome}_${ref}.alignment.primary.sorted.allInfo.tsv

#cat ./mergePgRef/${genome}_${ref}.alignment.primary.sorted.allInfo.tsv | head
wc -l ./mergePgRef/${genome}_${ref}.alignment.primary.sorted.allInfo.tsv
cat ./mergePgRef/${genome}_${ref}.alignment.primary.sorted.allInfo.tsv | cut -f1 | sort | uniq | wc -l


###############
# end message #
###############
echo "##########################################################"
echo "##########################################################"
echo "##########################################################"

cgroup_dir=$(awk -F: '{print $NF}' /proc/self/cgroup)
peak_mem=`cat /sys/fs/cgroup$cgroup_dir/memory.peak`
echo [$(date +"%Y-%m-%d %H:%M:%S")] peak memory is $(echo $peak_mem | numfmt --to=iec) \($peak_mem\) bytes
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds

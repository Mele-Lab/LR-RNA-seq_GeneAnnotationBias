#!/bin/bash

##################
# slurm settings #
##################

# For array
#SBATCH --array=7-12

# where to put stdout / stderr
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# time limit in minutes
#SBATCH --time=24:00:00

# queue
#SBATCH --qos=long

# memory (MB)
#SBATCH --mem=50G

# job name
#SBATCH --job-name lift_isoTr

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

line=`sed "1q;d" liftTrIso_cluster.list`
line=`sed "$((SLURM_ARRAY_TASK_ID))q;d" liftTrIso_cluster.list`

genome=$(echo $line | cut -f1 -d" ")
ref=$(echo $line | cut -f2 -d" ")

echo $genome
echo $ref

reference="[PATH]/data/genome/${genome}/${genome}.fa"
gtf="[PATH]/mapping_minimap/RNAseqLR_mapping/${genome}.fq_${genome}.fa/isoquant/OUT/OUT.transcript_models.gtf"
target="[PATH]/data/genome/${ref}/${ref}.fa"
output="${genome}_${ref}"

mkdir ./results/$output
cd ./results/$output


echo $reference
echo $gtf

source activate
conda activate liftoff

echo $target
echo $output

## Liftover
liftoff \
$target \
$reference \
-g $gtf \
-o $output


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
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) secondse eg	g� 

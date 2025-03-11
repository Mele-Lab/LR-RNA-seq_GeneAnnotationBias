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
#SBATCH --time=1:00:00

# queue
#SBATCH --qos=shorter

# memory (MB)
#SBATCH --mem=10G

# job name
#SBATCH --job-name gffcompare

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
gffcompare="[PATH]/tools/gffcompare/gffcompare"
gffread="[PATH]/tools/gffread/gffread"

###############################################
# Submit array according to a list/file       #
# slurm arrays are 0 based, sed line counting #
# starts from 1                               #
###############################################

line=`sed "$((SLURM_ARRAY_TASK_ID))q;d" "[PATH]/mapping_minimap/isoquant_analysis/mapIsoTrToRef/liftTrIso_cluster.list"`

genome=$(echo $line | cut -f1 -d" ")
ref=$(echo $line | cut -f2 -d" ")

echo $genome
echo $ref

mkdir results/${genome}_${ref}
cd results/${genome}_${ref}


##Test data
GTF2="[PATH]/no_backup_rg/mapping_minimap/isoquant_analysis/mapIsoTrToRef/results/${genome}_${ref}/${genome}_${ref}.gtf"
GTF1="[PATH]/data/annotation/241018_v47_poder_merge.placeholder_gene_name.withCDS.trBiotype.gtf"
#fasta=""
output="${genome}_GencodePoder"


## Prelimianry steps 
grep -v "^#" $GTF1 | awk '{if ($3 == "transcript" || $3 == "exon") print $0}' > GTF1.tmp
grep -v "^#" $GTF2 | awk '{if ($3 == "transcript" || $3 == "exon") print $0}' > GTF2.tmp

## Run gffcompare
$gffcompare GTF2.tmp \
-r GTF1.tmp \
-o $output \
-V



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
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) secondse eg	g? 





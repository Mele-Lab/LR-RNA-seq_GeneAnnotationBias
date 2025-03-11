#!/bin/bash

##################
# slurm settings #
##################

# For array
#SBATCH --array=1-18

# where to put stdout / stderr
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# time limit in minutes
#SBATCH --time=05:00:00

# queue
#SBATCH --qos=normal

# memory (MB)
#SBATCH --mem=50G

# job name
#SBATCH --job-name fullReadsInfo

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
module load SAMtools/1.20-GCC-13.2.0


###############################################
# Submit array according to a list/file       #
# slurm arrays are 0 based, sed line counting #
# starts from 1                               #
###############################################

line=`sed "1q;d" createFullReadsInfo_cluster.list`
line=`sed "$((SLURM_ARRAY_TASK_ID))q;d" createFullReadsInfo_cluster.list`

genome=$(echo $line | cut -f1 -d" ")
ref=$(echo $line | cut -f2 -d" ")

echo $genome
echo $ref

###############
# run command #
###############



########### Extraction of Bam information
echo "Extraction of Bam information"
samtools view "[PATH]/mapping_minimap/RNAseqLR_mapping/${genome}.fq_${ref}.fa/alignment.primary.sorted.bam" | awk '
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
    print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" nm "\t" ms "\t" as "\t" tp "\t" s1 "\t" s2 "\t" de
  }
' > "[PATH]/mapping_minimap/RNAseqLR_mapping/${genome}.fq_${ref}.fa/alignment.primary.sorted.infoFromBam.tsv"

echo "test1"
#cat "[PATH]/mapping_minimap/RNAseqLR_mapping/${genome}.fq_${ref}.fa/alignment.primary.sorted.infoFromBam.tsv" | head
echo "test2"
wc -l "[PATH]/mapping_minimap/RNAseqLR_mapping/${genome}.fq_${ref}.fa/alignment.primary.sorted.infoFromBam.tsv"
cat "[PATH]/mapping_minimap/RNAseqLR_mapping/${genome}.fq_${ref}.fa/alignment.primary.sorted.infoFromBam.tsv" | cut -f1 | sort | uniq | wc -l

########### Extraction of Bed information
echo "Extraction of Bed information"
cat "[PATH]/mapping_minimap/RNAseqLR_mapping/${genome}.fq_${ref}.fa/alignment.primary.sorted.bed" | awk '{
  split($11, arr, ",")
  sum_ex=0
  for (i in arr) sum_ex+=arr[i]
  sum_in=$3-$2+1-sum_ex
  print $4 "\t" $10 "\t" sum_ex "\t" sum_in
}' >  [PATH]/mapping_minimap/RNAseqLR_mapping/${genome}.fq_${ref}.fa/alignment.primary.sorted.formattedExonInfo.bed

#cat [PATH]/mapping_minimap/RNAseqLR_mapping/${genome}.fq_${ref}.fa/alignment.primary.sorted.formattedExonInfo.bed | head
wc -l [PATH]/mapping_minimap/RNAseqLR_mapping/${genome}.fq_${ref}.fa/alignment.primary.sorted.formattedExonInfo.bed
cat [PATH]/mapping_minimap/RNAseqLR_mapping/${genome}.fq_${ref}.fa/alignment.primary.sorted.formattedExonInfo.bed | cut -f1 | sort | uniq | wc -l


########### Merging bam/bed information
echo "Merging (1)"
python3 [PATH]/mapping_minimap/readFollowing/mergeReads_step1.py \
  [PATH]/mapping_minimap/RNAseqLR_mapping/${genome}.fq_${ref}.fa/alignment.primary.sorted.formattedExonInfo.bed \
  [PATH]/mapping_minimap/RNAseqLR_mapping/${genome}.fq_${ref}.fa/alignment.primary.sorted.infoFromBam.tsv \
  [PATH]/mapping_minimap/RNAseqLR_mapping/${genome}.fq_${ref}.fa/alignment.primary.sorted.mappingInfo.tsv

#cat [PATH]/mapping_minimap/RNAseqLR_mapping/${genome}.fq_${ref}.fa/alignment.primary.sorted.mappingInfo.tsv | head
wc -l [PATH]/mapping_minimap/RNAseqLR_mapping/${genome}.fq_${ref}.fa/alignment.primary.sorted.mappingInfo.tsv
cat [PATH]/mapping_minimap/RNAseqLR_mapping/${genome}.fq_${ref}.fa/alignment.primary.sorted.mappingInfo.tsv | cut -f1 | sort | uniq | wc -l


########### Collapsing reads/tr info
echo "Collasping reads/tr"
python3 [PATH]/mapping_minimap/readFollowing/collapseReadsTr.py \
  "[PATH]/mapping_minimap/RNAseqLR_mapping/${genome}.fq_${ref}.fa/isoquant/OUT/OUT.transcript_model_reads.tsv.gz" \
  "[PATH]/mapping_minimap/RNAseqLR_mapping/${genome}.fq_${ref}.fa/isoquant/OUT/OUT.transcript_model_reads.collapsed.tsv.gz"

#zcat [PATH]/mapping_minimap/RNAseqLR_mapping/${genome}.fq_${ref}.fa/isoquant/OUT/OUT.transcript_model_reads.collapsed.tsv.gz | head
zcat [PATH]/mapping_minimap/RNAseqLR_mapping/${genome}.fq_${ref}.fa/isoquant/OUT/OUT.transcript_model_reads.collapsed.tsv.gz | wc -l 
zcat [PATH]/mapping_minimap/RNAseqLR_mapping/${genome}.fq_${ref}.fa/isoquant/OUT/OUT.transcript_model_reads.collapsed.tsv.gz | cut -f1 | sort | uniq | wc -l



########### Merging bam/bed and reads/tr info
echo "Merging (2)"
python3 [PATH]/mapping_minimap/readFollowing/mergeReads_step2.py \
  [PATH]/mapping_minimap/RNAseqLR_mapping/${genome}.fq_${ref}.fa/alignment.primary.sorted.mappingInfo.tsv \
  [PATH]/mapping_minimap/RNAseqLR_mapping/${genome}.fq_${ref}.fa/isoquant/OUT/OUT.transcript_model_reads.collapsed.tsv.gz \
  [PATH]/mapping_minimap/RNAseqLR_mapping/${genome}.fq_${ref}.fa/alignment.primary.sorted.allInfo.tsv

#cat [PATH]/mapping_minimap/RNAseqLR_mapping/${genome}.fq_${ref}.fa/alignment.primary.sorted.allInfo.tsv | head
wc -l [PATH]/mapping_minimap/RNAseqLR_mapping/${genome}.fq_${ref}.fa/alignment.primary.sorted.allInfo.tsv
cat [PATH]/mapping_minimap/RNAseqLR_mapping/${genome}.fq_${ref}.fa/alignment.primary.sorted.allInfo.tsv | cut -f1 | sort | uniq | wc -l


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

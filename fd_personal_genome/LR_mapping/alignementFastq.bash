#!/bin/bash
#$ -q short-centos79,long-centos79,rg-el7
#$ -pe smp 8
#$ -l virtual_free=50G
#$ -l h_rt=12:00:00
#$ -cwd
#$ -N fastqAlignement
#$ -j y
#$ -t 1-20:1
#$ -o logs/fastqAlignement-$JOB_ID-$TASK_ID.out.log
#$ -e logs/fastqAlignement-$JOB_ID$TASK_ID.err.log


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
file="[PATH]/mapping_minimap/RNAseqLR_mapping/mapping_queries.tsv"
line=`sed "$((SGE_TASK_ID+1))q;d" ${file}`
fastq=$(echo $line | awk '{print $1}')
fasta=$(echo $line | awk '{print $2}')
output_dir=$(echo $line | awk '{print $3}')

###############
# run command #
###############
## Variables (fixed)
threads=12
minimap2="[PATH]/tools/minimap2-2.28_x64-linux/minimap2"
echo $fastq
echo $fasta
echo $output_dir

#### Functions
# Directory creation
mkdir -p $output_dir
cd $output_dir

# Align the genome
$minimap2 \
--MD \
-x splice \
-t ${threads} \
--secondary=yes \
-L \
-a ${fasta} \
${fastq} > alignment.sam

samtools view -bS alignment.sam > alignment.bam
samtools sort alignment.bam -o alignment.sorted.bam

# Remove multi-mapping reads
#samtools view -H alignment.sorted.bam > alignmentWoMulti.sorted.bam 
#samtools view -F 256 alignment.sorted.bam  >> alignmentWoMulti.sorted.bam 

bedtools bamtobed -i alignment.sorted.bam -bed12 > alignment.sorted.bed
#bedtools bamtobed -i alignmentWoMulti.sorted.bam -bed12 > alignmentWoMulti.sorted.bed

###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds

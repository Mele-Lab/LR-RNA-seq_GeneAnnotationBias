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
#SBATCH --time=24:00:00

# queue
#SBATCH --qos=long

# memory (MB)
#SBATCH --mem=50G

# job name
#SBATCH --job-name blast

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

source activate
conda activate blast
gffread="[PATH]/tools/gffread/gffread"


###############################################
# Submit array according to a list/file       #
# slurm arrays are 0 based, sed line counting #
# starts from 1                               #
###############################################

line=`sed "1q;d" blastAnalysis_cluster.list`
line=`sed "$((SLURM_ARRAY_TASK_ID))q;d" blastAnalysis_cluster.list`

genome=$(echo $line | cut -f1 -d" ")
ref=$(echo $line | cut -f2 -d" ")

echo $genome
echo $ref

mkdir -p results/${genome}_${ref}
cd results/${genome}_${ref}


## Create FASTA for the transcripts created by isqouant in eahc of the genome
## Ref
FASTA="[PATH]/data/genome/${genome}/${genome}.fa"
GTF="[PATH]/mapping_minimap/RNAseqLR_mapping/${genome}.fq_${genome}.fa/isoquant/OUT/OUT.transcript_models.gtf"

$gffread \
-w transcripts_PG.fasta \
-g ${FASTA} \
${GTF}

## Personal genome 
FASTA="[PATH]/data/genome/${ref}/${ref}.fa"
GTF="[PATH]/mapping_minimap/RNAseqLR_mapping/${genome}.fq_${ref}.fa/isoquant/OUT/OUT.transcript_models.gtf"

$gffread \
-w transcripts_${ref}.fasta \
-g ${FASTA} \
${GTF}

## Extract transcript_id from fasta file
grep ">" transcripts_PG.fasta | cut -f2 -d">"  > transcript_id_PG.txt
grep ">" transcripts_${ref}.fasta | cut -f2 -d">" > transcript_id_${ref}.txt

# ## Blast alignment
# blastn \
# -query transcripts_PG.fasta \
# -subject transcripts_${ref}.fasta \
# -outfmt 6 \
# -out alignment.tsv

# blastn \
# -query transcripts_PG.fasta \
# -subject transcripts_${ref}.fasta \
# -outfmt 6 \
# -max_hsps 1 \
# -max_target_seqs 1 \
# -out alignment_bestHit.tsv

# cut -f1 alignment_bestHit.tsv | sort | uniq > tmp.txt
# grep -v -f tmp.txt transcript_id_PG.txt > transcript_id_PG_notFoundBlast_${ref}.txt
# rm tmp.txt

## Blast alignment (other way)
blastn \
-subject transcripts_PG.fasta \
-query transcripts_${ref}.fasta \
-outfmt 6 \
-out alignment_fromRef.tsv

blastn \
-subject transcripts_PG.fasta \
-query transcripts_${ref}.fasta \
-outfmt 6 \
-max_hsps 1 \
-max_target_seqs 1 \
-out alignment_bestHit_fromRef.tsv

cut -f1 alignment_bestHit_fromRef.tsv | sort | uniq > tmp.txt
grep -v -f tmp.txt transcript_id_${ref}.txt > transcript_id_${ref}_notFoundBlast_PG.txt
rm tmp.txt


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

#!/bin/bash
#$ -q short-centos79,long-centos79,rg-el7
#$ -pe smp 8
#$ -l virtual_free=8G
#$ -l h_rt=0:20:00
#$ -cwd
#$ -N kmerCounts
#$ -j y
#$ -t 1-59:1
#$ -o logs/kmerCounts-randomRegion-$JOB_ID-$TASK_ID.out.log
#$ -e logs/kmerCounts-randomRegion-$JOB_ID$TASK_ID.err.log


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
###############################################
file="[PATH]/mapping_minimap/repetitionQuantification/randomRegion/kmerCount_randomRegion.list"
line=`sed "$((SGE_TASK_ID))q;d" ${file}`
#fasta="/nfs/no_backup/rg/fdegalez/mapping_minimap/repetitionQuantification/randomRegion/shuffledFASTA/HG002_HG00621.shuffled.fasta"
fasta=$(echo $line | awk '{print $1}')


###############
# run command #
###############
## Variables (fixed)
echo $fasta

cd [PATH]/mapping_minimap/repetitionQuantification/randomRegion/kmers

#### Functions
kmer_size=6

name=$(basename $fasta | cut -d "." -f1)


## Create spec directory
mkdir -p $name
cd $name

## Create sub directory to split the contigs
mkdir -p contigs
cd contigs
python3 [PATH]/mapping_minimap/repetitionQuantification/split_fasta.py $fasta

## One line fasta

for contigs in $(ls *.fasta)
do
echo $contigs
grep -v ">" $contigs | sed "s/\n//g" > tmp && mv tmp $contigs
done

cd ..



## Create sub direcotry to extract k-mers frequency
mkdir -p kmer${kmer_size}

for contigs in $(ls contigs/*.fasta)
do
echo $contigs
name=$(basename $contigs | sed "s/.fasta//g")
python3 [PATH]/script/kmer_calc_fa.py $contigs 6 | sort -rnk2,2  > kmer${kmer_size}/$name.txt
done


###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds

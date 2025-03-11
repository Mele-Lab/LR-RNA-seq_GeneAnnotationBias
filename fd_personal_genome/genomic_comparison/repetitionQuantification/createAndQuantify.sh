#!/usr/bin/sh

## Creation of random set of regions
module load BEDTools/2.30.0-GCC-11.3.0

for file in $(cat unique_regions_sup1000.list)
do
duo_name=$(dirname $file | awk -F "/" '{print $NF}') 
echo $duo_name
wc -l ${file}
FASTA_NAME=$(dirname $file | awk -F "/" '{print $NF}' | cut -d "_" -f1) 
FASTA="[PATH]/data/genome/${FASTA_NAME}/${FASTA_NAME}.size"

cat $file | awk '{if ($4<4000000) print $0}' | sort -rnk4,4 > tmp

bedtools shuffle \
-i tmp \
-g ${FASTA} \
-noOverlapping \
-seed 160524 > shuffledBED/${duo_name}.shuffled.bed
done

rm tmp


## Create FASTA from these files
for file in $(ls shuffledBED/*)
do
echo $file
basename=$(basename $file)
output=$(echo $basename | sed "s/.bed/.fasta/g")
fasta_name=$(echo $file | cut -d "/" -f2 | cut -d "_" -f1)
fasta="[PATH]/data/genome/${fasta_name}/${fasta_name}.fa"
bedtools getfasta \
-fi $fasta \
-bed $file \
-fo shuffledFASTA/$output
done
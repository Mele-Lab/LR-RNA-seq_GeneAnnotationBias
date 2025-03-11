#!/bin/bash
#$ -q short-centos79,long-centos79,rg-el7
#$ -pe smp 8
#$ -l virtual_free=50G
#$ -l h_rt=12:00:00
#$ -cwd
#$ -N genomeComparison
#$ -j y

# Set -e to stop on error
set -e
# Set -u to treat unset variables as an error
set -u

## Functions

# To extract the identifier from the file path
extract_id() {
  # Use basename to get the file name and remove the extension
  basename "$1" | cut -d. -f1
}

## Variables (as arguments)
fasta_query=$FASTA_QUERY
fasta_target=$FASTA_TARGET

## Variables (fixed)
threads=12
minimap2="[PATH]/tools/minimap2-2.28_x64-linux/minimap2"

## Output generation
id_query=$(extract_id "$fasta_query")
id_target=$(extract_id "$fasta_target")
output_dir=$(echo "${id_target}_${id_query}")

# Directory creation
mkdir -p $output_dir
cd $output_dir

# Loading module
module load BEDTools/2.30.0-GCC-11.3.0

# Align the genome
# asm5 : profile used by minimap for intraspecieds alignement of FASTA
# supposing less than 5% of divergence.
#
# Please see minimap2 documentation for more informations on the parameters.
$minimap2 -ax asm5 \
$fasta_target \
$fasta_query \
-t ${threads} \
--secondary=yes  > alignment.sam

# Convert SAM to BAM
samtools view -bS alignment.sam > alignment.bam
# Sort the BAM file
samtools sort alignment.bam -o alignment.sorted.bam

# Remove multi-mapping reads
samtools view -H alignment.sorted.bam > alignmentWoMulti.sorted.bam 
samtools view -F 256 -F 272 alignment.sorted.bam  >> alignmentWoMulti.sorted.bam 

# Generate a coverage BED file
bedtools genomecov -ibam alignment.sorted.bam -bga > genome_coverage.bed
bedtools genomecov -ibam alignmentWoMulti.sorted.bam  -bga > genome_coverage_WoMulti.bed

# Filter for zero-coverage regions and add size of it
awk '$4 == 0 {print $1"\t"$2"\t"$3"\t"$3-$2}' genome_coverage.bed > unique_regions.bed
awk '$4 == 0 {print $1"\t"$2"\t"$3"\t"$3-$2}' genome_coverage_WoMulti.bed > unique_regions_WoMulti.bed


#!/bin/bash

module load minimap2 samtools

FASTQ=/gpfs/projects/bsc83/Projects/pantranscriptome/pclavell/02_ONT_preprocessing/data/q7/${1}_preprocessed_Q7.fastq.gz
# CHANGE GTF ACCORDINGLY !!!!!!!!!!!!!!!!!
PATHGTF=/gpfs/projects/bsc83/Data/gene_annotations/gencode/v47/gencode.v47.primary_assembly.annotation.gtf
FASTA=/gpfs/projects/bsc83/Data/assemblies/GRCh38/modified/GRCh38.primary_assembly.sirvset4.genome.fa 
GTF=$(echo $(basename $PATHGTF) | sed 's/\.gtf//')
SAMPLENAME=$(echo $(basename $FASTQ) | sed 's/\._preprocessed_Q7.fastq.gz//')


echo "CREATING TRANSCRIPTOME"
/gpfs/projects/bsc83/utils/gffread/gffread \
    -w ref/$GTF.fa \
    -g $FASTA \
    $PATHGTF

# map
echo "MAPPING TO TRANSCRIPTOME"
mkdir -p data/transcriptome_mapped
(set -x;minimap2 -ax map-ont -t 112 ref/$GTF.fa $FASTQ |\
    samtools view -q 10 -Sb | \
    samtools sort -@112 - -o data/transcriptome_mapped/${SAMPLENAME}_transcriptome.bam)

samtools index -@112 -b data/transcriptome_mapped/${SAMPLENAME}_transcriptome.bam


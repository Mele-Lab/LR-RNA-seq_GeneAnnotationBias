#!/bin/bash

module load java-openjdk/22.0.1
module load gatk

REF_FASTA=/gpfs/projects/bsc83/Data/assemblies/GRCh38/modified/GRCh38.primary_assembly.sirvset4.genome.fa 
INPUT_BAM=/gpfs/projects/bsc83/Projects/pantranscriptome/pclavell/03_mapping/data/genomic/20240212_HS_2_PY2_GM10493.bam
OUTPUT_BAM=/gpfs/projects/bsc83/Projects/pantranscriptome/pclavell/05_variant_calling/data/20240212_HS_2_PY2_GM10493_splitN.bam
THREADS=112

### need to create a sequence dictionary for a reference FASTA
gatk CreateSequenceDictionary -R $REF_FASTA

### SplitNCigarReads
gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=$THREADS" SplitNCigarReads \
  -R $REF_FASTA \
  -I $INPUT_BAM \
  -O $OUTPUT_BAM
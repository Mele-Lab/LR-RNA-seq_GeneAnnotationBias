#!/bin/bash

module load minimap2 samtools \
    gatk R/4.3.2 vcftools bcftools \
    oneapi/2024.1 htslib/1.19.1 tabixpp/1.1.2

INPUT_ORIGINAL_BAM=/gpfs/projects/bsc83/Projects/pantranscriptome/pclavell/03_mapping/data/genomic/20240212_HS_2_PY2_GM10493.bam
INPUT_SPLIT_BAM=data/20240212_HS_2_PY2_GM10493_splitN.bam
OUTPUT_BAM=data/flagcorrected/20240212_HS_2_PY2_GM10493_splitN_corrected.bam
THREADS=112

Rscript /gpfs/projects/bsc83/utils/lrRNAseqVariantCalling/tools/flagCorrection.r \
  $INPUT_ORIGINAL_BAM \
  $INPUT_SPLIT_BAM \
  $OUTPUT_BAM \
  $THREADS
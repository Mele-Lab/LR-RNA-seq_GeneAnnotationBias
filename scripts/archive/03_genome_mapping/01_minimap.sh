#!/bin/bash

module load minimap2

QUERY=$1
#REFERENCE=$2

NAME=$(echo "$QUERY" | sed 's/_with_extracted.*//' | sed 's/.*\///')

echo "$NAME"

# index
#minimap2 -d 03_genome_mapping/gencodev44_transcriptome.mmi $REFERENCE

# map
minimap2 \
    -ax splice \
    --junc-bed /gpfs/projects/bsc83/Data/gene_annotation/gencode/release_44/gencode.v44.chr_patch_hapl_scaff.annotation.bed \
    -t 48 \
    -o 03_genome_mapping/01_alignment_results/"${NAME}"_juncbed_transcriptome.sam \
    03_genome_mapping/gencodev44_transcriptome.mmi \
    $QUERY

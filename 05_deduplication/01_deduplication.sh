#!/bin/bash

module load python/3.6.1

INPUT_BAM=$1
OUT_DIR=$2
EDIT_DISTANCE=2
NAME=$(echo "$INPUT_BAM" | sed 's/.bam//' | sed 's/.*\///')
echo $NAME

umi_tools dedup \
    --extract-umi-method read_id \
    --umi-separator _ \
    --method adjacency \
    --edit-distance-threshold $EDIT_DISTANCE \
    --per-contig \
    --per-gene \
    --gene-transcript-map 04_assess_dedup/gencodev44_transcript_map.tsv \
    --stdin $INPUT_BAM \
    --stdout $OUT_DIR/"$NAME"_"$EDIT_DISTANCE"_dedup.bam \
    --log $OUT_DIR/"$NAME"_"$EDIT_DISTANCE"_dedup.log
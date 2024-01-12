#!/bin/bash

module load python/3.6.1

QUERY=$1
EDIT_DISTANCE=$2
NAME=$(echo "$QUERY" | sed 's/.bam//' | sed 's/.*\///')
echo $NAME

umi_tools group \
    -i \
    --method adjacency \
    --edit-distance-threshold=$EDIT_DISTANCE \
    --per-contig \
    --per-gene \
    --gene-transcript-map 04_dedup/gencodev44_transcript_map.tsv \
    -I $QUERY \
    --group-out=04_dedup/01_dedup_results/"$NAME"_"$EDIT_DISTANCE"ED_percontig.tsv \
    --log=04_dedup/01_dedup_results/"$NAME"_"$EDIT_DISTANCE"ED_percontig.log
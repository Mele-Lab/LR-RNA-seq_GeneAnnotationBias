#!/bin/bash

module load python/3.6.1
module load samtools

QUERY=$1
NAME=$(echo "$QUERY" | sed 's/.sam//' | sed 's/.*\///')
echo $NAME

awk '$1 ~ /^@/ || $10 !~ /\*/' $QUERY > $TMPDIR/$NAME.filtered.sam

samtools view -b $TMPDIR/$NAME.filtered.sam > $TMPDIR/$NAME.bam
samtools sort $TMPDIR/$NAME.bam > $TMPDIR/$NAME.sorted.bam
samtools index -b $TMPDIR/$NAME.sorted.bam

umi_tools group \
    -i \
    --edit-distance-threshold=4 \
    --per-contig \
    --per-gene \
    --gene-transcript-map 04_dedup/gencodev44_transcript_map.tsv \
    -I $TMPDIR/$NAME.sorted.bam \
    --group-out=04_dedup/01_dedup_results/"$NAME"_4ED_percontig.tsv \
    --log=04_dedup/01_dedup_results/"$NAME"_4ED_percontig.log
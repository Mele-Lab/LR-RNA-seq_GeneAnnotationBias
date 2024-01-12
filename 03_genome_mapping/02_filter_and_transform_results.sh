#!/bin/bash

module load samtools

QUERY=$1

NAME=$(echo "$QUERY" | sed 's/.sam//' | sed 's/.*\///')
echo $NAME

# keep lines starting with @ or that do not have a * in 10th position
#awk '$1 ~ /^@/ || $10 !~ /\*/' $QUERY | awk '$1 ~ /^@/ || $2<2048 {print}' | awk '{split($1,a,"_"); if(length(a[2])==16|| $0 ~ /^@/) print $0}'> $TMPDIR/$NAME.filtered.sam
awk '$1 ~/^@/ {print $0}' $QUERY > 03_genome_mapping/01_alignment_results/$NAME.filtered.sam
awk '$2<2048 && $10 !~ /\*/' $QUERY | awk '{split($1,a,"_"); if(length(a[2])==16) print $0}'>> 03_genome_mapping/01_alignment_results/$NAME.filtered.sam


samtools view -b 03_genome_mapping/01_alignment_results/$NAME.filtered.sam > 03_genome_mapping/01_alignment_results/$NAME.bam
samtools sort 03_genome_mapping/01_alignment_results/$NAME.bam > 03_genome_mapping/01_alignment_results/$NAME.sorted.bam
samtools index -b 03_genome_mapping/01_alignment_results/$NAME.sorted.bam
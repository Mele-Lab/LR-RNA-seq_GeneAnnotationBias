#!/bin/bash

module load samtools
SAMPLE=$1

if [[ "$SAMPLE" == *"assembly"* ]]; then
  TYPE="assembly"
else
  TYPE="general"
fi
echo $TYPE

COUNT=$(samtools view -@112 $SAMPLE | cut -f1 | sort | uniq | wc -l)

echo $SAMPLE $COUNT >> data/${TYPE}_mapping/number_mapped_reads.txt
#!/bin/bash

SAMPLE=$1
SAMPLE_NAME=$(echo "$SAMPLE" | sed 's/.*\///' |sed 's/_.ED_percontig.tsv//')

for ed in {0,1,2,3,4}
    do cut -f8 04_dedup/01_dedup_results/"$SAMPLE_NAME"_"$ed"ED_percontig.tsv | sort | uniq -c >> 04_dedup/01_dedup_results/"$SAMPLE_NAME"_statistics.tsv
done

#!/bin/bash

for sample in $(ls -d data/2024*/*transcript_models.gtf)
    do echo -e "$sample\t$(awk '$3 == "transcript" {print $0}' $sample | wc -l)">> stats/count_discovered_transcripts.tsv
done
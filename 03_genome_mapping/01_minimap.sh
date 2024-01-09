#!/bin/bash

module load minimap2/2.24-r1122

REFERENCE=$1
QUERY=$2

# index
#minimap2 -d 03_genome_mapping/ref.mmi $REFERENCE

# map
minimap2 \
    -ax map-ont \
    -o 03_genome_mapping/sample.bam \
    03_genome_mapping/ref.mmi \
    $QUERY
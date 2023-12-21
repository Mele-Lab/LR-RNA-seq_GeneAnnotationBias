#!/bin/bash

module load minimap2

REF=$1
QUERY=$2

minimap2 $REF $QUERY \
    -map-ont \
    -a \
    -o output.sam \
    -t 48 \
    -A 2 -B 3 -O 3
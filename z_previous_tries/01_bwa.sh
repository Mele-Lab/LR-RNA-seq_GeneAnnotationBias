#!/bin/bash

module load bwa

REF=$1
QUERY=$2
OUTPUT=$3

bwa mem -t 48 -o $OUTPUT.sam \
    $REF $QUERY
#    -A 2 -B 4 -O [5,5] -E [1,1] -L [5,5] \

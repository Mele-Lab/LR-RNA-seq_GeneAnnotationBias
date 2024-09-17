#!/bin/bash

GTF=$1
OUTPUT=$(echo $GTF | sed 's/\.gtf//')
module load python

python /gpfs/projects/bsc83/utils/SUPPA-2.4/suppa.py generateEvents \
    -i $GTF \
    -o $OUTPUT \
    -f ioe \
    -e {SE,SS,MX,RI,FL} \
    --boundary V

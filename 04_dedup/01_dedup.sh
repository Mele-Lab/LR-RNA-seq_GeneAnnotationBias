#!/bin/bash

module load python/3.6.1

QUERY=$1

umi_tools group \
    -I $QUERY \
    --group-out=04_dedup/results.tsv \
    --log=04_dedup/group.log
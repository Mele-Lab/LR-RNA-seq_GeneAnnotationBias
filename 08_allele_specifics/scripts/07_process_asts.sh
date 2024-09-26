#!/bin/bash

INPUTS=$(ls data/06_calc_asts/*_asts.tsv)
OUTDIR=/data/07_processed_asts
mkdir -p $OUTDIR
process_asts \
    -i $INPUTS \
    -g genes.tsv \ ###### LIST OF GENES???
    -o $OUTDIR
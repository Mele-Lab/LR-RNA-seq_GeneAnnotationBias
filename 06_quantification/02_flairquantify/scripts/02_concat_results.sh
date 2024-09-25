#!/bin/bash

# add column with sample anem

INPUTDIR=data/240917_merge_geneEntry_correctedScaffolds_nochrEBV
for file in $(ls $INPUTDIR/*_stringent.counts.tsv)
    do OUTPUTFILE=$(basename $file | sed 's/\.counts\.tsv//')
    awk 'BEGIN{OFS="\t"} NR==1 {split($2, arr, "_"); print $0, arr[1]} NR>1 {print $0, arr[1]}' $file > $INPUTDIR/${OUTPUTFILE}_stringent_samplecol.tsv
    done

# concatenate
awk 'FNR==1 && NR!=1 { next } { print }' $INPUTDIR/*_stringent_samplecol.tsv > $INPUTDIR/concat_flair_stringent_counts.tsv

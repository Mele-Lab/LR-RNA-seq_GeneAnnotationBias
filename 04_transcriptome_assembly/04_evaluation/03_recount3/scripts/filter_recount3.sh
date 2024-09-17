#!/bin/bash

#wget https://snaptron.cs.jhu.edu/data/srav3h/junctions.bgz

zcat /gpfs/projects/bsc83/Projects/pantranscriptome/pclavell/04_transcriptome_assembly/04_evaluation/03_recount3/data/download_database/junctions.bgz |\
 cut -f2,3,4,6,7,8,9,13,14 |\
 grep -v "?" |\
 awk '$8>10{print $0"\t"$1":"$2"-"$3":"$4}' > recount3_srv3h_morethan10counts.tsv2

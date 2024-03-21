#!/bin/bash

module load blast

FASTA=$1
OUTPUT=$2

makeblastdb -in $FASTA -input_type fasta -dbtype nucl -out $OUTPUT
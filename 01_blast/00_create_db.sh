#!/bin/bash

module load blast

INPUT=$1
OUTPUT=$2

makeblastdb -in $INPUT -input_type fasta -dbtype nucl -out $OUTPUT
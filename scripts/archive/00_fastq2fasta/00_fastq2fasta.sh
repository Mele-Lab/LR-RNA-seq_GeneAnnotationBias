#!bin/bash

INPUT=$1
OUTPUT=$(echo "$INPUT" | sed 's/\.fastq/\.fasta/' | sed 's/.*\///')

echo "$OUTPUT"

sed -n '1~4s/^@/>/p;2~4p' $INPUT > 00_data_deduplication/$OUTPUT

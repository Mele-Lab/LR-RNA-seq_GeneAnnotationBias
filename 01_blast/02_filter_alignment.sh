#!/bin/bash

# Input file
INPUT=$1

# Output file
OUTPUT=$2

# Remove @ lines and rm BS so the number can be compared
#grep -v @ $INPUT | awk 'BEGIN{OFS="\t"}sub("BS:f:","",$16) '> $TMPDIR/input.sam

# Sort the file by query name (column 3) and alignment score (column 16)
sort -k1,1 -k11,11nr $INPUT | \
  awk 'BEGIN{OFS="\t"}{
    if (!seen[$1]++) {
      print;
    }
  }' | awk '$12>40 {print}'> $OUTPUT

echo "Filtered file saved to $OUTPUT"

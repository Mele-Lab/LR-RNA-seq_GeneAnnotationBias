#!/bin/bash

GTF=$1
MYPREFIX=$(echo $1 | sed 's-.*/--' | sed 's/\.isoforms\.gtf_classification\.txt//')
echo $MYPREFIX
#COUNT=$(cat $GTF| cut -f 6 | sort  | uniq -c | grep -v "category")
cat $GTF| cut -f 6 | sort  | uniq -c | grep -v "category" |\
    awk -v prefix="$MYPREFIX" '{print $0,prefix}' > data/${MYPREFIX}_count.txt


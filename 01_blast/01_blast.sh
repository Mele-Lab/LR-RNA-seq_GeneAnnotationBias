#!/bin/bash
module load blast/2.11.0
module load samtools

QUERY=$1
DATABASE=$2
OUTPUT_NAME=$3

blastn -query $QUERY \
    -db $DATABASE \
    -outfmt 5 \
    -out $OUTPUT_NAME \
    -num_threads 48 \
    -strand 'both' \
    -task 'blastn' \
    -gapopen 3 -gapextend 1 -penalty -2 -reward 1
    #default:    -gapopen 5 -gapextend 2 -penalty -3 -reward 1
    #-outfmt "17 SQ"\
        #-outfmt "6 qseqid qlen slen qstart qend sstart send qseq sseq evalue bitscore length pident nident mismatch positive gapopen gaps sstrand btop"\



#!/bin/bash

module load samtools

GENOMICBAM=$1

SAMPLE=$(echo $GENOMICBAM | sed 's-.*/--' | sed 's/\.bam//')
SIRVOMICBAM=data/assembly_mapping/sirv/${SAMPLE}_SIRV.bam


samtools view $GENOMICBAM | cut -f1 | awk '{print "@"$0}'> $TMPDIR/${SAMPLE}_genomic.readids
samtools view $SIRVOMICBAM | cut -f1 | awk '{print "@"$0}'> $TMPDIR/${SAMPLE}_sirvomic.readids
zcat ../02_ONT_preprocessing/data/q10/${SAMPLE}_preprocessed_Q10.fastq.gz |\
    awk 'NR % 4 == 1' > $TMPDIR/${SAMPLE}_all.readids



cat $TMPDIR/${SAMPLE}_genomic.readids $TMPDIR/${SAMPLE}_sirvomic.readids > $TMPDIR/${SAMPLE}_mapped.readids

grep -Fvxf $TMPDIR/${SAMPLE}_mapped.readids $TMPDIR/${SAMPLE}_all.readids > /gpfs/projects/bsc83/Projects/pantranscriptome/pclavell/03_mapping/data/assembly_mapping/unmapped_reads/${SAMPLE}_unmapped.txt

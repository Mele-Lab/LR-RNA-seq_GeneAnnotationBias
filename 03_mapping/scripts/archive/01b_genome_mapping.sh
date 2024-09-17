#!/bin/bash

module load minimap2
module load samtools

INPUTFASTQ=$1
PSEUDOGENEMASKED=_unmasked_unfiltered
JUNCBED=/gpfs/projects/bsc83/Data/gene_annotations/gencode/v47/modified/gencode.v47.primary_assembly.sirvset4.annotation.junc.bed
# INDEX=ref/GRCh38.primary_assembly.sirvset4.genome.masked_pseudogenes.mmi
INDEX=ref/GRCh38.primary_assembly.sirvset4.genome.mmi



# automatically parse previous info for the outputs
SAMPLE=$(echo $INPUTFASTQ | sed 's-.*/--' | sed 's/_preprocessed\.fastq\.gz//')
OUTSAM=data$PSEUDOGENEMASKED/$SAMPLE.sam
OUTBAM=data$PSEUDOGENEMASKED/$SAMPLE.bam
OUTSIRVBAM=data$PSEUDOGENEMASKED/sirv/${SAMPLE}_SIRV.bam


minimap2 \
    -ax splice \
    --junc-bed $JUNCBED \
    -u n \
    -t 112 \
    --MD \
    -o $OUTSAM \
    -a $INDEX \
    $INPUTFASTQ


# convert to bam, filter by mapping Q >=10 and included in genome contigs
samtools view -@112 -b $OUTSAM |\
    samtools sort -@112 -o $OUTBAM
samtools index -@112 -b $OUTBAM


rm $OUTSAM















# glinos -ax splice -uf -k14 –secondary=no parameters
# inamo 2024: default parameters and reference to the splice junctions in the GENCODE annotation
# lrgasp: minimap2 -ax splice

#full-length cDNAs, it would be desired to apply -u f to force minimap2 to consider the forward transcript strand only.
# # from minimap2 paper
# If RNA-seq reads are not sequenced from stranded libraries, 
# the read strand relative to the underlying transcript is unknown. 
# By default, minimap2 aligns each chain twice, 
# first assuming GT–AG as the splicing signal and 
# then assuming CT–AC, the reverse complement of GT–AG,
#  as the splicing signal. The alignment with a higher score 
#  is taken as the final alignment. This procedure also infers
#   the relative strand of reads that span canonical splicing sites.

# I decide to don't check only one strand although it is slower (its the default also)

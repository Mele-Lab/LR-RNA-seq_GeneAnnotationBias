#!/bin/bash

module load minimap2
module load samtools

INPUTFASTQ=$1
PSEUDOGENEMASKED=_pseudomasked_mapq10filtered
JUNCBED=/gpfs/projects/bsc83/Data/gene_annotations/gencode/v47/modified/gencode.v47.primary_assembly.sirvset4.annotation.junc.bed
INDEX=ref/GRCh38.primary_assembly.sirvset4.genome.masked_pseudogenes.mmi
#INDEX=ref/GRCh38.primary_assembly.sirvset4.genome.mmi



# automatically parse previous info for the outputs
SAMPLE=$(echo $INPUTFASTQ | sed 's-.*/--' | sed 's/_preprocessed\.fastq\.gz//')
OUTSAM=data$PSEUDOGENEMASKED/$SAMPLE.sam
OUTBAM=data$PSEUDOGENEMASKED/genomic/$SAMPLE.bam
OUTSIRVBAM=data$PSEUDOGENEMASKED/sirv/${SAMPLE}_SIRV.bam

# create dirs
mkdir data$PSEUDOGENEMASKED
mkdir data$PSEUDOGENEMASKED/genomic
mkdir data$PSEUDOGENEMASKED/sirv

minimap2 \
    -ax splice \
    --junc-bed $JUNCBED \
    -t 112 \
    --MD \
    -o $OUTSAM \
    -a $INDEX \
    $INPUTFASTQ


# # count number of reads mapped
# COUNT=$(samtools view -@112 -F 4 $OUTSAM | cut -f1 | sort | uniq | wc -l)
# echo $SAMPLE $COUNT >> data$PSEUDOGENEMASKED/number_of_mapped_reads.txt

# convert to bam, filter by mapping Q >=10 and included in genome contigs
samtools view -@112 -b --min-MQ 10 -L ref/GRCh38_primary_assembly_contigs.bed $OUTSAM |\
    samtools sort -@112 -o $OUTBAM
samtools index -@112 -b $OUTBAM
# convert to bam, filter by mapping Q >=10 and included in sirv contigs
samtools view -@112 -b --min-MQ 10 -L ref/sirv_erc_contigs.bed $OUTSAM |\
    samtools sort -@112 -o $OUTSIRVBAM
samtools index -@112 -b $OUTSIRVBAM

rm $OUTSAM





# # Count reads in each subsetted BAM file
# READS_IN_GENOMICBAM=$(samtools view -@112 $OUTBAM | cut -f1 | sort | uniq | wc -l)
# READS_IN_SIRVBAM=$(samtools view -@112 $OUTSIRVBAM | cut -f1 | sort | uniq | wc -l)

# # Count reads in FASTQ.GZ file
# # Every 4 lines represent one read in a FASTQ file # fastq count is not working because of the sample names
# FASTQ_FILE="/gpfs/projects/bsc83/Projects/pantranscriptome/pclavell/02_ONT_preprocessing/scripts/data/$SAMPLE/${SAMPLE}_preprocessed.fastq.gz"
# READS_IN_FASTQ=$(( $(zcat $FASTQ_FILE | wc -l) / 4 ))


# # Generate output
# echo -e "$SAMPLE\t$READS_IN_FASTQ\t$COUNT\t$READS_IN_GENOMICBAM\t$READS_IN_SIRVBAM" >> data$PSEUDOGENEMASKED/count_mapped_reads.tsv

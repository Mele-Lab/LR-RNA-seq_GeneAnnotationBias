#!/bin/bash

module load anaconda samtools
source activate base
conda init
conda activate /gpfs/projects/bsc83/utils/conda_envs/flair/

GENOMEFASTA=/gpfs/projects/bsc83/Data/assemblies/GRCh38/GRCh38.primary_assembly.genome.fa
ANNOTATIONGTF=/gpfs/projects/bsc83/Data/gene_annotations/gencode/v47/gencode.v47.primary_assembly.annotation.gtf

BAM=$1
ITERATION=$2
SAMPLENAME=$(echo $1 | sed 's-.*/--' | sed 's/\.bam//')
ORIGINALSAMPLE=$(echo $SAMPLENAME | sed 's/_downsampled.*//')
mkdir data/flair
# convert alignment bam to bed so flair can correct the splicing junctions
mkdir data/flair/00_beds
mkdir data/flair/00_beds/iter_$ITERATION

# sort and index bam
samtools sort -@112 -o /gpfs/projects/bsc83/Projects/pantranscriptome/pclavell/04_transcriptome_assembly/05_downsampling/data/bam/${SAMPLENAME}_sorted.bam $BAM
samtools index -@112 /gpfs/projects/bsc83/Projects/pantranscriptome/pclavell/04_transcriptome_assembly/05_downsampling/data/bam/${SAMPLENAME}_sorted.bam

SORTEDBAM=/gpfs/projects/bsc83/Projects/pantranscriptome/pclavell/04_transcriptome_assembly/05_downsampling/data/bam/${SAMPLENAME}_sorted.bam

/gpfs/projects/bsc83/utils/conda_envs/flair/bin/bam2Bed12 -i $SORTEDBAM > data/flair/00_beds/iter_$ITERATION/$SAMPLENAME.bed12




# FLAIR CORRECT
echo "FLAIR CORRECT"

mkdir data/flair/01_correct
mkdir data/flair/01_correct/iter_$ITERATION

/gpfs/projects/bsc83/utils/conda_envs/flair/bin/flair correct \
    -q data/flair/00_beds/iter_$ITERATION/$SAMPLENAME.bed12 \
    -f $ANNOTATIONGTF \
    -g $GENOMEFASTA \
    --ss_window 8 \
    --threads 112 \
    --output data/flair/01_correct/iter_$ITERATION/$SAMPLENAME

# FLAIR COLLAPSE
echo "FLAIR COLLAPSE"

mkdir data/flair/02_collapse
mkdir data/flair/02_collapse/iter_$ITERATION


#extract readID
samtools view -@112 $BAM | cut -f1 > $TMPDIR/read_ids$SAMPLENAME.txt

mkdir data/fastq

if [ ! -e data/fastq/${ORIGINALSAMPLE}.fastq ]; then
  zcat /gpfs/projects/bsc83/Projects/pantranscriptome/pclavell/02_ONT_preprocessing/data/q10/${ORIGINALSAMPLE}_preprocessed_Q10.fastq.gz > data/fastq/${ORIGINALSAMPLE}.fastq
fi

# subset fastq
/gpfs/projects/bsc83/utils/seqtk/seqtk subseq \
    data/fastq/${ORIGINALSAMPLE}.fastq \
    $TMPDIR/read_ids$SAMPLENAME.txt > $TMPDIR/${SAMPLENAME}_subset.fastq




/gpfs/projects/bsc83/utils/conda_envs/flair/bin/flair collapse \
    -g $GENOMEFASTA \
    --gtf $ANNOTATIONGTF \
    -q data/flair/01_correct/iter_$ITERATION/${SAMPLENAME}_all_corrected.bed \
    -r $TMPDIR/${SAMPLENAME}_subset.fastq \
    --stringent \
    --check_splice \
    --generate_map \
    --annotation_reliant generate \
    --threads 112 \
    --output data/flair/02_collapse/iter_$ITERATION/$SAMPLENAME

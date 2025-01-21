#!/bin/bash

module load anaconda bedops
source activate base
conda init
conda activate /gpfs/projects/bsc83/utils/conda_envs/flair/
#PATHGTF=/gpfs/projects/bsc83/Projects/pantranscriptome/novelannotations/merged/240917_merge_geneEntry_correctedScaffolds_nochrEBV.gtf
PATHGTF=/gpfs/projects/bsc83/Projects/pantranscriptome/novelannotations/merged/poder_v1.sorted.gtf
FASTA=/gpfs/projects/bsc83/Data/assemblies/GRCh38/modified/GRCh38.primary_assembly.sirvset4.genome.fa 
GTF=$(echo $(basename $PATHGTF) | sed 's/\.gtf//')

# echo "CREATING TRANSCRIPTOME"
# /gpfs/projects/bsc83/utils/gffread/gffread \
#     -w ref/$GTF.fa \
#     -g $FASTA \
#     $PATHGTF

# echo "CREATE ISOFORM.bed"
# perl scripts/gtf2bed12.pl $PATHGTF > ref/$GTF.bed


echo "QUANTIFYING"
OUTNAME=$(echo $(basename $1) | sed 's/\.txt//')

mkdir -p data/poder


/gpfs/projects/bsc83/utils/conda_envs/flair/bin/flair quantify \
    -r $1 \
    --isoforms ref/$GTF.fa \
    --output data/poder/${OUTNAME}_stringent \
    --threads 112 \
    --sample_id_only \
    --isoform_bed ref/$GTF.bed \
    --stringent






# Rcode to create array file to launch in mn5
# test <- fread("../02_ONT_preprocessing/arrayq10_fastqgz", header = F)
# test[, lab_sampleid:=tstrsplit(V1,"_")[[4]]]
# myfile <-metadata[, .(lab_sampleid, sample)][test, on="lab_sampleid"]
# fwrite(myfile[, `:=`(condition="condition1", batch="batch1")][, .(sample, condition, batch, V1)],
#        "../06_quantification/01_isoquantify/array_flair", row.names = F, col.names = F, quote = F, sep="\t")
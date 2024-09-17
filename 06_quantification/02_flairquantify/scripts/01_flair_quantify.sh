#!/bin/bash

module load anaconda
source activate base
conda init
conda activate /gpfs/projects/bsc83/utils/conda_envs/flair/
PATHGTF=/gpfs/projects/bsc83/Projects/pantranscriptome/novelannotations/merged/240917_merge_geneEntry_correctedScaffolds_nochrEBV.gtf
FASTA=/gpfs/projects/bsc83/Data/assemblies/GRCh38/modified/GRCh38.primary_assembly.sirvset4.genome.fa 
GTF=$(echo $(basename $PATHGTF) | sed 's/\.gtf//')

# echo "CREATING TRANSCRIPTOME"
# /gpfs/projects/bsc83/utils/gffread/gffread \
#     -w ref/$GTF.fa \
#     -g $FASTA \
#     $PATHGTF

echo "QUANTIFYING"
OUTNAME=$(echo $(basename $1) | sed 's/\.txt//')

mkdir -p data/240917_merge_geneEntry_correctedScaffolds_nochrEBV


/gpfs/projects/bsc83/utils/conda_envs/flair/bin/flair quantify \
    -r $1 \
    --isoforms ref/$GTF.fa \
    --output data/240917_merge_geneEntry_correctedScaffolds_nochrEBV/all_samples_counts \
    --threads 112
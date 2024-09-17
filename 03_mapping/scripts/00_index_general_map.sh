#!/bin/bash

module load minimap2

INDEX=ref/GRCh38.primary_assembly.sirvset4.genome.mmi
FASTA=/gpfs/projects/bsc83/Data/assemblies/GRCh38/modified/GRCh38.primary_assembly.sirvset4.genome.fa 


minimap2 -x map-ont -t 112 -d $INDEX $FASTA

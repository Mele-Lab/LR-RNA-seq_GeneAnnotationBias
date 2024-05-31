#!/bin/bash

module load miniconda
source activate sqanti3-snakemake
umi_tools group \
    --method adjacency \
    --edit-distance-threshold=2 \
    --umi-separator , \
    --per-gene \
    --per-contig \
    -I data/tenpercentbam/minimap/tenpercentbam_sorted.bam \
    --group-out umitest1_groupout \
    --log umitest2_log
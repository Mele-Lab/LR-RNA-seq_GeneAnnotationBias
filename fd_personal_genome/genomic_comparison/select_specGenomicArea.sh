#!/bin/bash
# Select genomic areas with a size superior or equal to a threshold (thr_size)
# and collapse the overlapping regions.

# Load the BEDTools module
module load BEDTools/2.30.0-GCC-11.3.0

# Set the threshold size
thr_size=1000

# Select the regions with a size superior or equal to the threshold
for file in $(ls */unique_regions.bed)
do
    # Get the directory path
    path="$(dirname $file)" 
    echo $file
    # Filter the regions with a size superior or equal to the threshold
    awk -v thr=${thr_size} '{if ($4 >= thr) print $0}' $file > $path/"unique_regions_sup${thr_size}.bed"
done

module load BEDTools/2.30.0-GCC-11.3.0
# Collapse the overlapping regions
for file in $(ls */unique_regions_sup${thr_size}.bed)
do
    # Get the directory path
    path="$(dirname $file)" 
    echo $file
    # Sort the regions by chrom, start and end
    sort -k1,1 -k2,2n $file > tmp.sorted.bed
    # Collapse the overlapping regions and sum their size
    bedtools merge -i tmp.sorted.bed  -c 4 -o sum > $path/"unique_regions_sup${thr_size}.collapsed.bed"
done

#!/bin/bash

POD5_DIR=$1
SAMPLE_NAME=$(dirname ${POD5_DIR} | awk '{split($1, parsedpath,"/"); print parsedpath[9]}')

module load dorado/0.5.3-linux-x64 cuda/12.2 nvidia-hpc-sdk/23.11

dorado duplex \
	-t 80 \
	scripts/dna_r10.4.1_e8.2_400bps_sup@v4.3.0 \
	$POD5_DIR > data/$SAMPLE_NAME.bam
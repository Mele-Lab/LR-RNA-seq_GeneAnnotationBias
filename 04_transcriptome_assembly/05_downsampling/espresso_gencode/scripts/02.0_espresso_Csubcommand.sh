#!/bin/bash

module load anaconda
source activate espresso


SAMPLENAME=$1
NUM=$2

perl /apps/GPP/ANACONDA/2023.07/envs/espresso/bin/ESPRESSO_C.pl \
    -I data/espresso_c/$SAMPLENAME/$NUM \
    -F data/espresso_c/$SAMPLENAME/fastas/$NUM.fa \
    -X 0 -T 5
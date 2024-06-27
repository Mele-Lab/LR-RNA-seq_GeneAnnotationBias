#!/bin/bash

module load java-openjdk/17.0.11+9 gatk

gatk MergeVcfs \
    -I $INPUT_VCFs \
    -O $OUTPUT_VCF \
    --REFERENCE_SEQUENCE $REFERENCE

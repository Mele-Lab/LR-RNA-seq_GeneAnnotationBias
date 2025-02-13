#!/bin/bash

SAMPLE=$1
zcat $SAMPLE | wc -l > ${SAMPLE}_count
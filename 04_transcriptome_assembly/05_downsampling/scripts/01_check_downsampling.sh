#!/bin/bash

SAMPLE=$1

module load samtools
samtools view $1 | wc -l
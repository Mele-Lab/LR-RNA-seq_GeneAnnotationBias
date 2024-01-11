#!/bin/bash

INPUT=$1 # sam file

module load python/3.10.2

python3 02_extract_UMI/01_extract_UMI.py $INPUT
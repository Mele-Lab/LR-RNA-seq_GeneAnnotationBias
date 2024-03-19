#!/bin/bash

INPUT=$1 # filtered sam file
NAME=$(echo "$XML" | sed 's/.*\///' | sed 's/\.filtered.sam//')

module load python/3.10.2

python3 02_extract_UMI/01_extract_UMI.py "$INPUT" "$NAME"
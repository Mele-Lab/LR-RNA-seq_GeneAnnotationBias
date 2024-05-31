#!/bin/bash

module load dorado

dorado trim \
	--threads 80 \
	--emit-fastq \
	-v \
       data/tenpercentbam/tenpercentbam.fastq > outputtrim.fastq	

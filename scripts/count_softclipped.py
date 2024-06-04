
import re
import pandas as pd
import sys
import numpy as np
import os

os.getcwd()

sam_path="/gpfs/projects/bsc83/Projects/pantranscriptome/pclavell/ONT_preprocessing/scripts/minimap_output_postporechop_split_params.sam"
import pysam
samfile = pysam.AlignmentFile(sam_path, "r")


softclipping=[]
for read in samfile:
    mycigar=read.cigartuples
    for i in mycigar:
        if i[0]==4:
            softclipping.append(i[1])

np.savetxt('/gpfs/projects/bsc83/Projects/pantranscriptome/pclavell/ONT_preprocessing/scripts/softclipping_post_split_params.txt', softclipping, fmt='%d')
#!/bin/python3

# This code assumes that the SAM file is filtered to include ONLY:
    # 1 alignment per read
    # 16 nt UMI containing alignments

import re
import pandas as pd
import numpy as np
import sys

#-DEFINITIONS--------------------------------------------------------------

# Parse cigar
def parse_cigar(cigar):
    pattern = re.compile(r'(\d+)([S=XHID])')
    parsed_cigar = pattern.findall(cigar)
    return parsed_cigar

# detect potential concatamers
def is_concatamer(cigar):
    parsed_cigar = parse_cigar(cigar)
    logical = ['S' in alignment_type for alignment_type in parsed_cigar]
    sub_cigar = [item for item, include in zip(parsed_cigar, logical) if include]
    count = int(0)
    for item in sub_cigar: 
        if int(item[0])>40:
            count += 1
    if count > 1:
        return True
    return False

cigar='444S33=16X20=718S'
extract_UMI(cigar)


# Functions definitions
def find_UMI_location(cigar):
    # parse cigar
    parsed_cigar = parse_cigar(cigar)

    start_position = 0
    end_position = 0

    # assess how many read bases are before the UMI
    for count, alignment_type in parsed_cigar:
        count = int(count)
        if alignment_type in ['S', '=', 'H', 'I']:
            start_position += count
        elif alignment_type == 'X':
            if count > 10:
                break
            elif count <= 10:
                start_position += count
    end_position = start_position+16
    return start_position, end_position

# extract UMI sequence
def extract_UMI(df):
    UMI_location = find_UMI_location(df['cigar'])
    UMI = df['seq'][UMI_location[0]:UMI_location[1]]
    result = [[df['read_name']],
                [df['seq']],
                [UMI],
                [UMI_location[0]],
                [UMI_location[1]],
                [df['cigar']],
                [is_concatamer(df['cigar'])]]
    return(result)

#---------------------------------------------------------------


# Read filtered sam file (unmapped reads have been filtered out)
sam_path = "01_blast/01_results/top100res.filtered.sam"
sam_path = sys.argv[0]

df = pd.read_csv(sam_path, sep='\t', skiprows=2)
df.columns=['read_name', 'flag', 'ref_name', 'pos', 'mapq', 'cigar', 'rnext', 'pnext', 'template_length', 'seq', 'qual', 'x1', 'x2', 'x3', 'x4']
print(df)

# for each read, obtain the sequence of the UMI
umi_series= df.apply(extract_UMI, axis=1)

# convert series to df
umi_table = umi_series.apply(lambda x: pd.Series({'read_name': x[0], 'read_seq': x[1], 'umi_seq': x[2], 'umi_start' : x[3], 'umi_end' : x[4], 'cigar': x[5], 'potential_concatamer': x[6]}))

# convert lists inside of dataframe to strings
umi_df = umi_table.map(lambda x: x[0] if len(x) > 0 else x)

# remove rows that for some reason miss the read name
umi_df = umi_df[umi_df['read_seq'] != '*']

# drop potential concatamers
umi_df = umi_df[~ umi_df['potential_concatamer']]

# add column with read_name + UMIseq
umi_df.insert(1, 'read_name_UMI', '>' + umi_df['read_name'] + '_' + umi_df['umi_seq'])

# keep only unique reads (reads with only one alignment containing a UMI)
umi_df = umi_df[~umi_df['read_name'].duplicated()]



# save table
umi_df.to_csv('02_extract_UMI/reads_and_umis.tsv', sep='\t', index=False)

# now create a fasta file
data_to_create_fasta= umi_df[['read_name_UMI', 'read_seq']]
data_to_create_fasta.to_csv('/home/pclavell/mounts/projects/Projects/gencode_diversity/deduplication/02_extract_UMI/umi_extracted.fasta', sep='\n', index=False, header=False)



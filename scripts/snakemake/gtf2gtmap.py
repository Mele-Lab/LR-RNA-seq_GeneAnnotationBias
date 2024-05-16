import pyranges as pr
import sys

inputgtf = sys.argv[1]
outputgt_map = sys.argv[2]

df = pr.read_gtf(inputgtf).df
df = df.loc[df.transcript_id.notnull()]
df = df[['gene_id', 'transcript_id']].drop_duplicates()
df.to_csv(outputgt_map, header=None, index=False, sep='\t')
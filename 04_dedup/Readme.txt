# Pau Clavell

|--gencodev44_transcript_map.tsv
awk '/>/{split($0,a,"|"); print a[2],a[1]}' /gpfs/projects/bsc83/Data/gene_annotation/gencode/release_44/gencode.v44.transcripts.fa | sed 's/ >/\t/' 
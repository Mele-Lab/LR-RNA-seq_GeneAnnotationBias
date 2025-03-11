#!/usr/bin/sh


conda activate sqanti3

sq3="[PATH]/tools/SQANTI3/sqanti3_qc.py"

genome="HG002"

for genome in HG002 HG00621 HG01928 NA18906 NA19240 HG01952
do
echo $genome

ref="hg38"

GTF="[PATH]/mapping_minimap/isoquant_analysis/mapIsoTrToRef/results/${genome}_${ref}/${genome}_${ref}.gtf"
reference="[PATH]/data/annotation/241018_v47_poder_merge.placeholder_gene_name.withCDS.trBiotype.gtf"
fasta="[PATH]/data/genome/hg38/hg38.fa"
output="test"

python $sq3 \
    $GTF \
    $reference \
    $fasta \
    -d ${genome}_${ref} \
    --report skip \
    --force_id_ignore \
    --aligner_choice minimap2 \
    --skipORF \
    -o $output

done
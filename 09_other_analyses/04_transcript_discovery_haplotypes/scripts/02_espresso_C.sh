#!/bin/bash
module load anaconda
source activate espresso

FILEPATH=$1
SAMPLENAME=$2

GENOMEFASTA=/gpfs/projects/bsc83/Data/assemblies/GRCh38/GRCh38.primary_assembly.genome.fa
ANNOTATIONGTF=/gpfs/projects/bsc83/Data/gene_annotations/gencode/v47/gencode.v47.primary_assembly.annotation.gtf

mkdir -p data/espresso_c/$SAMPLENAME

#
echo "0---FIRST STEP---0"
python /gpfs/projects/bsc83/utils/conda_envs/espresso/snakemake/scripts/split_espresso_s_output_for_c.py \
    --orig-work-dir data/espresso_s/$SAMPLENAME \
    --new-base-dir data/espresso_c/$SAMPLENAME \
    --target-reads-per-c 500000 \
    --num-threads-per-c 5 \
    --genome-fasta $GENOMEFASTA

for i in $(seq 0 $(ls data/espresso_c/$SAMPLENAME | grep -Eo '[1-9]*' | wc -l)); do
    echo -e "bash scripts/02.0_espresso_Csubcommand.sh ${SAMPLENAME} ${i}" >> data/espresso_c/$SAMPLENAME/greasyfile.tsv
done

THREADS=5
bash /gpfs/projects/bsc83/utils/slurm_wrapper_pau_mn5.sh \
    -j espressoC_${SAMPLENAME} \
    -N $(echo "($(cat data/espresso_c/$SAMPLENAME/greasyfile.tsv | wc -l)*$THREADS + 112 - 1) / 112" | bc) \
    -tn $(echo $((112/${THREADS}))) \
    -u $THREADS \
    -g True \
    -gf data/espresso_c/$SAMPLENAME/greasyfile.tsv \
    -c "" \
    -t 3:20:00


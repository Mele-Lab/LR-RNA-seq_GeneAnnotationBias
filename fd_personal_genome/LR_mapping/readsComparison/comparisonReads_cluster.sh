#!/bin/bash

##################
# slurm settings #
##################

# For array
#SBATCH --array=1-12

# where to put stdout / stderr
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# time limit in minutes
#SBATCH --time=02:30:00

# queue
#SBATCH --qos=shorter

# memory (MB)
#SBATCH --mem=40G

# job name
#SBATCH --job-name comparisonReads

# cpu slots
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

#################
# start message #
#################
start_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] starting on $(hostname)

# make bash behave more robustly
set -e
set -u
set -o pipefail

###################
# set environment #
###################


###############################################
# Submit array according to a list/file       #
# slurm arrays are 0 based, sed line counting #
# starts from 1                               #
###############################################

line=`sed "4q;d" mergePgRef.list`
line=`sed "$((SLURM_ARRAY_TASK_ID))q;d" mergePgRef.list`

genome=$(echo $line | cut -f1 -d" ")
ref=$(echo $line | cut -f2 -d" ")

echo $genome
echo $ref

file="[PATH]/mapping_minimap/readFollowing/fullReadsInfo/mergePgRef/${genome}_${ref}.alignment.primary.sorted.allInfo.tsv"
mkdir mergePgRef/diffStats/${genome}_${ref}
cd mergePgRef/diffStats/${genome}_${ref}

###############
# run command #
###############
echo "statsAlignement"

cat $file | awk '
BEGIN {
    print "readId \t nbBlocks_diff \t sizeBlock_diff \t sizeGap_diff \t FLAG_diff \t MAPQ_diff \t nbMismatchAndGaps_diff \t dpScoreMaxScoringAlgt_diff \t dpAlgtScore_diff \t typeOfAlgt_diff  \t chainingScore1_diff \t chainingScore2_diff \t gapCompressedPerBaseSeqDiv_diff \t transcriptPart"
}

{
    if ($10 != "NA" && $27 != "NA") {
        readId= $1
        nbBlocks_diff= $2-$19
        sizeBlock_diff= $3-$20
        sizeGap_diff= $4-$21
        FLAG_diff= $5-$22
        MAPQ_diff= $8-$25
        nbMismatchAndGaps_diff=$10-$27
        dpScoreMaxScoringAlgt_diff=$11-$28
        dpAlgtScore_diff=$12-$29
        typeOfAlgt_diff=$13-$30
        chainingScore1_diff=$14-$31
        chainingScore2_diff=$15-$32
        gapCompressedPerBaseSeqDiv_diff=$16-$33
        transcriptPart="none"
        if ($18 != "*" && $35 != "*"){transcriptPart="both"} else if ($18 != "*" && $35 == "*"){transcriptPart="PG"} else if ($18 == "*" && $35 != "*"){transcriptPart="ref"}

        print readId "\t" nbBlocks_diff "\t" sizeBlock_diff "\t" sizeGap_diff "\t" FLAG_diff "\t" MAPQ_diff "\t" nbMismatchAndGaps_diff "\t" dpScoreMaxScoringAlgt_diff "\t" dpAlgtScore_diff "\t" typeOfAlgt_diff  "\t" chainingScore1_diff "\t" chainingScore2_diff "\t" gapCompressedPerBaseSeqDiv_diff "\t" transcriptPart
    }
}' > statsAlignment.tsv

#Focus on CIGAR
echo "CIGAR"
cat $file | cut -f9 > tmp1.tmp
python3 [PATH]/mapping_minimap/readFollowing/fullReadsInfo/cigar_parser.py tmp1.tmp > tmp1_cigar.tmp
cat $file | cut -f26 > tmp2.tmp
python3 /[PATH]/mapping_minimap/readFollowing/fullReadsInfo/cigar_parser.py tmp2.tmp > tmp2_cigar.tmp

echo "readId" > readId.tmp
cat $file | cut -f1 >> readId.tmp

paste tmp1_cigar.tmp tmp2_cigar.tmp readId.tmp > tmp3.tmp

cat tmp3.tmp | awk '
BEGIN {
  print "readId\tM_diff\tI_diff\tD_diff\tN_diff\tS_diff\tH_diff\tP_diff\teq_diff\tX_diff\ttotal_diff\tunaligned_diff\tcoverage_diff"
}

NR > 1{ 
    if ($12 != "NA" && $24 != "NA") {
        M_diff=$1-$13
        I_diff=$2-$14
        D_diff=$3-$15
        N_diff=$4-$16
        S_diff=$5-$17
        H_diff=$6-$18
        P_diff=$7-$19
        eq_diff=$8-$20
        X_diff=$9-$21
        total_diff=$10-$22
        unaligned_diff=$11-$23
        coverage_diff=$12-$24
    }

    print $25 "\t" M_diff "\t" I_diff "\t" D_diff "\t" N_diff "\t" S_diff "\t" H_diff "\t" P_diff "\t" eq_diff "\t" X_diff "\t" total_diff "\t" unaligned_diff "\t" coverage_diff }' > statsCIGAR.tsv


rm *tmp*
echo "merge"
python3 [PATH]/mapping_minimap/readFollowing/fullReadsInfo/mergeDiffTable.py statsAlignment.tsv statsCIGAR.tsv statsAll.tsv

for type in none both PG ref; do
echo $type
cat statsAll.tsv | cut -f1,2,3,4,26 > statsAll_reduced.tsv
grep "$type" statsAll.tsv | cut -f1,2,3,4,26 > statsAll_${type}.tsv
done


cut -f1,7 statsAlignment.tsv > tmp1.tsv
cut -f1,13 statsCIGAR.tsv > tmp2.tsv

## Merge with python

grep -v "NA" statsAll.tsv  | cut -f1,7,26 > missmatchesGapCovDiff.tsv
cat  missmatchesGapCovDiff.tsv | cut -f2,3 | sort | uniq -c > missmatchesGapCovDiffTable.tsv

#0034a20e-13db-4b82-a70f-037288e34963:0
#463S95M466S  1024tot
###############
# end message #
###############
echo "##########################################################"
echo "##########################################################"
echo "##########################################################"

cgroup_dir=$(awk -F: '{print $NF}' /proc/self/cgroup)
peak_mem=`cat /sys/fs/cgroup$cgroup_dir/memory.peak`
echo [$(date +"%Y-%m-%d %H:%M:%S")] peak memory is $(echo $peak_mem | numfmt --to=iec) \($peak_mem\) bytes
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds

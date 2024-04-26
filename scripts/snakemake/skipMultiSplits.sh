#!/usr/bin/bash

#input=/users/project/gencode_006070/gkaur/masterTablev2/noHiSS/table/artifacts/duplex_tools/
#splitStep=/users/project/gencode_006070_no_backup/gkaur/splitSamples_ONT_JointCadapter

input=$1 # it is like a wd
splitStep=$2 # directory basename with the splittings
sample=$3

for num in 1step 2step; do
	python getMultiReads.py $splitStep\_$num/${sample}/split_multiple_times.pkl $splitStep\_$num/${sample}/detected_multiSPLIT_IDs
	perl -pi -e "s/', '/\n/g" $splitStep\_$num/${sample}/detected_multiSPLIT_IDs
	perl -pi -e "s/{'|'}//g" $splitStep\_$num/${sample}/detected_multiSPLIT_IDs
	echo "DONE $sample"
done

mkdir $input/strategy3_analysis/${sample}/
cat $splitStep\_1step/${sample}/detected_multiSPLIT_IDs $splitStep\_2step/${sample}/detected_multiSPLIT_IDs \| sort | uniq > $input/strategy3_analysis/${sample}/detected_multiSPLIT_IDs_step1_2

zgrep "None" $splitStep\_2step/${sample}/*_split_split.fastq.gz | grep "_._" | awk -F"_" '{print $1}' | sort --parallel=4 -T $input/TEMP/ | uniq | awk '{print substr($1,2); }' > $input/strategy3_analysis/${sample}/UNdetected_multiSPLIT_IDs_step1_2


cat $input/strategy3_analysis/${sample}/detected_multiSPLIT_IDs_step1_2 $input/strategy3_analysis/${sample}/UNdetected_multiSPLIT_IDs_step1_2 | sort --parallel=4 -T $input/TEMP/ | uniq > $input/strategy3_analysis/${sample}/FINAL_multiSPLIT_IDs_step1_2_unq

gunzip $splitStep\_2step/${sample}/*_split_split.fastq.gz
perl $input/script/removeMultiSplitReads.pl $input/strategy3_analysis/${sample}/FINAL_multiSPLIT_IDs_step1_2_unq $splitStep\_2step/${sample}/*_split_split.fastq $splitStep\_2step/${sample}/multiSplitReads.fastq $splitStep\_2step/${sample}/splitReads.fastq
echo "DONE $sample"


gzip $splitStep\_2step/${sample}/*_split_split.fastq $splitStep\_2step/${sample}/multiSplitReads.fastq $splitStep\_2step/${sample}/splitReads.fastq
echo "gzipped"

unzipped_lines=`wc -l $splitStep\_2step/${sample}/splitReads.fastq`
echo -e "${unzipped_lines}"


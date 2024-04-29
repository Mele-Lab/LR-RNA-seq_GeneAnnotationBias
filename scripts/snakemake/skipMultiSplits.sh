#!/usr/bin/bash

#scriptwd=/users/project/gencode_006070/gkaur/masterTablev2/noHiSS/table/artifacts/duplex_tools/
#splitStep=/users/project/gencode_006070_no_backup/gkaur/splitSamples_ONT_JointCadapter

scriptwd=$1 # it is like a wd, were are all the scripts
outdir=$2
splitStep=$3 # directory basename with the splittings
sample=$4 # sample name
threads=$5

for num in 1step 2step; do
	python $scriptwd/getMultiReads.py $splitStep\_$num/${sample}/split_multiple_times.pkl \
	$splitStep\_$num/${sample}/detected_multiSPLIT_IDs.txt
	perl -pi -e "s/', '/\n/g" $splitStep\_$num/${sample}/detected_multiSPLIT_IDs.txt
	perl -pi -e "s/{'|'}//g" $splitStep\_$num/${sample}/detected_multiSPLIT_IDs.txt
	echo "DONE $sample"
done

mkdir $outdir/${sample}/
cat $splitStep\_1step/${sample}/detected_multiSPLIT_IDs.txt $splitStep\_2step/${sample}/detected_multiSPLIT_IDs.txt |\
 sort | uniq > $outdir/${sample}/detected_multiSPLIT_IDs_step1_2.txt

zgrep "None" $splitStep\_2step/${sample}/*_split_split.fastq.gz | grep "_._" |\
 awk -F"_" '{print $1}' | sort --parallel=$threads -T $scriptwd/TEMP/ |\
 uniq | awk '{print substr($1,2); }' > $outdir/${sample}/UNdetected_multiSPLIT_IDs_step1_2.txt

mkdir $outdir/TEMP
cat $outdir/${sample}/detected_multiSPLIT_IDs_step1_2.txt $outdir/${sample}/UNdetected_multiSPLIT_IDs_step1_2.txt | sort --parallel=$threads -T $outdir/TEMP | uniq > $outdir/${sample}/FINAL_multiSPLIT_IDs_step1_2_unq.txt

gunzip $splitStep\_2step/${sample}/*_split_split.fastq.gz

perl $scriptwd/removeMultiSplitReads.pl $outdir/${sample}/FINAL_multiSPLIT_IDs_step1_2_unq.txt $splitStep\_2step/${sample}/*_split_split.fastq $splitStep\_2step/${sample}/multiSplitReads.fastq $splitStep\_2step/${sample}/splitReads.fastq
echo "DONE $sample"


gzip $splitStep\_2step/${sample}/*_split_split.fastq $splitStep\_2step/${sample}/multiSplitReads.fastq $splitStep\_2step/${sample}/splitReads.fastq
echo "gzipped"

unzipped_lines=`zcat $splitStep\_2step/${sample}/splitReads.fastq.gz | wc -l`
echo -e "${unzipped_lines}"


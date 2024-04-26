#!/usr/bin/bash
#$ -N getMultiSplits
#$ -t 1-44
#$ -tc 44
#$ -q rg-el7,long-centos79,short-centos79
#$ -cwd
#$ -l virtual_free=32G,h_rt=2:00:00
#$ -pe smp 1
#$ -V
#$ -e logs_lines/$TASK_ID.err
#$ -o logs_lines/$TASK_ID.out

#readme_joint.sh #split using strategy3; split1 (COMPLETE adapter) -> result fed to split2 (JointCapTrapAdapter)
input=/users/project/gencode_006070/gkaur/masterTablev2/noHiSS/table/artifacts/duplex_tools/
splitStep=/users/project/gencode_006070_no_backup/gkaur/splitSamples_ONT_JointCadapter
#splitStep_2=/users/project/gencode_006070_no_backup/gkaur/splitSamples_ONT_JointCadapter_1step/

#rm -r ${splitStep1}/tagged_multisplits/
#mkdir ${splitStep1}/tagged_multisplits/ ${splitStep2}/tagged_multisplits/

file=$(ls -1  /users/project/gencode_006070/gkaur/masterTablev2/noHiSS/table/artifacts/duplex_tools/INPUT_FILES/*.gz | sed -n ${SGE_TASK_ID}p)

sample=`basename $file`

for num in 1step 2step; do
#	rm -r $splitStep\_$num/${sample}/detected_multiSPLIT_IDs
	python ${input}/script/getMultiReads.py $splitStep\_$num/${sample}/split_multiple_times.pkl $splitStep\_$num/${sample}/detected_multiSPLIT_IDs
	perl -pi -e "s/', '/\n/g" $splitStep\_$num/${sample}/detected_multiSPLIT_IDs
	perl -pi -e "s/{'|'}//g" $splitStep\_$num/${sample}/detected_multiSPLIT_IDs
	echo "DONE $sample"
done

mkdir $input/strategy3_analysis/${sample}/
cat $splitStep\_1step/${sample}/detected_multiSPLIT_IDs $splitStep\_2step/${sample}/detected_multiSPLIT_IDs | sort | uniq > $input/strategy3_analysis/${sample}/detected_multiSPLIT_IDs_step1_2

zgrep "None" $splitStep\_2step/${sample}/*_split_split.fastq.gz | grep "_._" | awk -F"_" '{print $1}' | sort --parallel=4 -T $input/TEMP/ | uniq | awk '{print substr($1,2); }' > $input/strategy3_analysis/${sample}/UNdetected_multiSPLIT_IDs_step1_2

#cat $input/strategy3_analysis/${sample}/detected_multiSPLIT_IDs_step1_2 $input/strategy3_analysis/${sample}/UNdetected_multiSPLIT_IDs_step1_2 > $input/strategy3_analysis/${sample}/FINAL_multiSPLIT_IDs_step1_2

cat $input/strategy3_analysis/${sample}/detected_multiSPLIT_IDs_step1_2 $input/strategy3_analysis/${sample}/UNdetected_multiSPLIT_IDs_step1_2 | sort --parallel=4 -T $input/TEMP/ | uniq > $input/strategy3_analysis/${sample}/FINAL_multiSPLIT_IDs_step1_2_unq

gunzip $splitStep\_2step/${sample}/*_split_split.fastq.gz
perl $input/script/removeMultiSplitReads.pl $input/strategy3_analysis/${sample}/FINAL_multiSPLIT_IDs_step1_2_unq $splitStep\_2step/${sample}/*_split_split.fastq $splitStep\_2step/${sample}/multiSplitReads.fastq $splitStep\_2step/${sample}/splitReads.fastq
echo "DONE $sample"

##TESTS
#toBeRemoved=`cat $input/strategy3_analysis/${sample}/FINAL_multiSPLIT_IDs_step1_2_unq | awk -F"_" '{print $1}' | sort --parallel=4 -T $input/TEMP/ | uniq | wc -l`
#getRunID=`head -1 $splitStep\_2step/${sample}/multiSplitReads.fastq | awk -F"/" '{print $2}' | awk -F"_" '{print $1}'`
#actuallyRemoved=`grep $getRunID $splitStep\_2step/${sample}/multiSplitReads.fastq | awk -F"_" '{print $1}' | sort --parallel=4 -T $input/TEMP/ | uniq | wc -l`
#echo -e "Tested remove:$toBeRemoved $getRunID removed:$actuallyRemoved"

gzip $splitStep\_2step/${sample}/*_split_split.fastq $splitStep\_2step/${sample}/multiSplitReads.fastq $splitStep\_2step/${sample}/splitReads.fastq
echo "gzipped"

unzipped_lines=`wc -l $splitStep\_2step/${sample}/splitReads.fastq`
echo -e "${unzipped_lines}"

#perl -pi -e "s/', '/\n/g" | perl -pi -e "s/{'|'}//g" ${splitStep1}/${sample}/detected_multiSPLITs 

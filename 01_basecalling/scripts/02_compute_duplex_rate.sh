module load samtools

BAM=$1
FILENAME=$(echo $BAM | sed 's-.*data/--' | sed 's/\.bam//')

mkdir data/stats

echo $FILENAME	"duplex"	$(samtools view -c -@112 --tag dx:1 $BAM) > data/stats/${FILENAME}_duplexrate.tsv
echo $FILENAME	"all"	$(samtools view -c -@112 $BAM) >> data/stats/${FILENAME}_duplexrate.tsv
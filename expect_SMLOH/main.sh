#!/usr/bin/bash

#$ -S /usr/bin/bash
#$ -l s_vmem=12G
#$ -cwd
#$ -j y
#$ -t 1-500:1 
#$ -o ~/tcga_data/somatic_conversion/log2/

#SGE_TASK_ID=4


line=$(($SGE_TASK_ID+1))
pid=`head -n $line mut_num_top500_with_bam.tsv|tail -n 1|awk '{print $1}'`
tid=`head -n $line mut_num_top500_with_bam.tsv|tail -n 1|awk '{print $2}'`
tbam=`head -n $line mut_num_top500_with_bam.tsv|tail -n 1|awk '{print $3}'`
nid=`head -n $line mut_num_top500_with_bam.tsv|tail -n 1|awk '{print $4}'`
nbam=`head -n $line mut_num_top500_with_bam.tsv|tail -n 1|awk '{print $5}'`

echo "start $pid"
mkdir $pid
cd $pid
if [ -e count_triplet_num.txt ]; then
	if [ -e coveraged.txt ]; then
		gzip -f coveraged.txt
	fi
	exit 1
fi

~/softs/gdc-client download -t ~/tcga_data/gdc-user-token.2022-03-18T06_54_36.824Z.txt $tid $nid

TUMOR="${tid}/${tbam}"
NORM="${nid}/${nbam}"
echo "start samtools depth"
/usr/local/package/samtools/0.1.20/bin/samtools depth $NORM $TUMOR|awk '$2>8&&$3>6' >coveraged.txt
#perl ../depth2bed.pl coveraged.txt >coverage.bed
#sh ../exon_list/mk_intersected.sh coverage.bed

echo "start intersect"
perl ../count_coverage.pl $pid

gzip coveraged.txt
#rm -rf $tid $nid

exit 0

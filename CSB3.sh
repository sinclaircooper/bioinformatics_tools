#! /bin/sh

FILES=$1/*".fastq.gz"

for f in $FILES
do 
	
	CSB3=$(zcat $f | egrep 'GGGGTTGGTGT|ACACCAACCCC|GGGGTTAGTGT|ACACTAACCCC' -c)
	total_reads=$(zcat $f | wc -l | awk '{print $1/4}')
	echo "$f\t$CSB3\t$total_reads"
		
done



#!/bin/bash
#set -e
#echo "Usage: edup.sh <alignment.bam> <threads> <outdir> <final_bam>"

bam=$1
threads=$2
out_dir=$3
final_bam=$4

bam_file=`basename $bam`
bam_prefix=`echo $bam_file | sed 's/.bam//g'`
sort_bam=${out_dir}/sort.tmp.bam
dd_bam=${out_dir}/dd.tmp.bam
resort_bam=${out_dir}/${bam_prefix}.dedup.bam
echo "# Sorting by coordinates"
echo "\
samtools sort -@$threads -T $sort_bam.tmp -m2G -O bam -o $sort_bam $bam"
samtools sort -@$threads -T $sort_bam.tmp -m2G -O bam -o $sort_bam $bam

#echo "samtools index $sort_bam"
#samtools index $sort_bam
#echo ""

echo "\
picard MarkDuplicates REMOVE_DUPLICATES=true I=$sort_bam O=$dd_bam \
M=${out_dir}/dd.metrics.tmp ASSUME_SORT_ORDER=coordinate MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1024"
#picard MarkDuplicates REMOVE_DUPLICATES=true I=$sort_bam O=$dd_bam \
#	M=${out_dir}/dd.metrics.tmp ASSUME_SORT_ORDER=coordinate \
#	MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1024
picard MarkDuplicates REMOVE_SEQUENCING_DUPLICATES=true I=$sort_bam O=$dd_bam \
	M="${out_dir}/dd.metrics.tmp" ASSUME_SORT_ORDER=coordinate \
	MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1024
echo ""
#END

echo "# Resort to name order"
echo "\
samtools sort -@$threads -n -T $resort_bam.tmp -m2G -O bam -o $resort_bam $dd_bam"
samtools sort -@$threads -n -T $resort_bam.tmp -m2G -O bam -o $resort_bam $dd_bam
ln -s $resort_bam $final_bam
echo ""

echo "# Cleanup dedup intermediate files"
echo "\
rm $sort.bam; rm $dd_bam; rm ${out_dir}/dd.metrics.tmp"
rm $sort_bam
rm $dd_bam
rm ${out_dir}/dd.metrics.tmp


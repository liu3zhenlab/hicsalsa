#!/bin/bash
#echo "Usage: dedup.sh <alignment.bam> <threads> <outdir> <final_bam> <script_log>"

bam=$1
threads=$2
out_dir=$3
final_bam=$4
script_log=$5

bam_file=`basename $bam`
bam_prefix=`echo $bam_file | sed 's/.bam//g'`
sort_bam=${out_dir}/sort.tmp.bam
dd_bam=${out_dir}/dd.tmp.bam
resort_bam=${out_dir}/${bam_prefix}.dedup.bam

echo "# Sorting by coordinates"
echo "\
samtools sort -@$threads -T $sort_bam.tmp -m2G -O bam -o $sort_bam $bam" >> $script_log
samtools sort -@$threads -T $sort_bam.tmp -m2G -O bam -o $sort_bam $bam
echo >> $script_log

echo "\
picard MarkDuplicates REMOVE_DUPLICATES=true I=$sort_bam O=$dd_bam \
M=${out_dir}/dd.metrics.tmp ASSUME_SORT_ORDER=coordinate MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1024" >> $script_log
picard MarkDuplicates REMOVE_DUPLICATES=true I=$sort_bam O=$dd_bam \
	M="${out_dir}/dd.metrics.tmp" ASSUME_SORT_ORDER=coordinate \
	MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1024
echo "" >> $script_log
#END

echo "# Resort to name order"
echo "\
samtools sort -@$threads -n -T $resort_bam.tmp -m2G -O bam -o $resort_bam $dd_bam" >> $script_log
samtools sort -@$threads -n -T $resort_bam.tmp -m2G -O bam -o $resort_bam $dd_bam
echo >> $script_log

echo "ln -s $resort_bam $final_bam" >> $script_log
ln -s $resort_bam $final_bam
echo "" >> $script_log

echo "# Cleanup dedup intermediate files"
echo "rm $sort.bam; rm $dd_bam; rm ${out_dir}/dd.metrics.tmp" >> $script_log
rm $sort_bam
rm $dd_bam
rm ${out_dir}/dd.metrics.tmp


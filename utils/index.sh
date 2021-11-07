#!/bin/sh

asm_path=$1
asm_file=$2
index_dir=$3
script_log=$4

cd ${index_dir}

ln -s -f ${asm_path}/${asm_file} .

echo "### step 0. genome indexing" >> $script_log

echo "bwa index $asm_file" >>$script_log
bwa index $asm_file
echo "" >>$script_log

echo "samtools faidx $asm_file" >>$script_log
samtools faidx $asm_file
echo "" >>$script_log

cd $wd


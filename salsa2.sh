#!/bin/sh
#Usage: salsa.sh <asm_path> <asm_file> <bam> <iteration> <minilen> <enzbase> <genomesize> <outdir>

asm_path=$1 # path to the assembly, only the path
asm_file=$2 # assembly fasta file, no path
bam=$3 # bam, including path
iter=$4 # iteration
minilen=$5 # minimal length of contigs for scaffolding
enzbase=$6  # enzbaseyme base
genomesize=$7 # estimated genome size
gfa=$8
filter_ctg=$9 # file to contain congigs to be removed
index_dir="${10}" # path of indexed assembly
outdir="${11}"    # output directory
salsa_path="${12}" # path to salsa scripts
script_log="${13}" # script log file

bam_prefix=`basename $bam | sed 's/.bam//g'`
bed=${outdir}/${bam_prefix}.bed

echo "### convert BAM to BED" >>$script_log
if ! [ -s $bed ]; then
	echo "bedtools bamtobed -i $bam > $bed" >>$script_log
	bedtools bamtobed -i $bam > $bed
fi

option_para=""
if ! [ "$genomesize" = false ]; then
	option_para="${option_para} -s $genomesize"
fi

if ! [ "$enzbase" == "false" ]; then
	option_para="${option_para} -e $enzbase"
fi

if ! [ "$gfa" == "false" ]; then
	option_para="${option_para} -g $gfa"
fi

if ! [ "$filter_ctg" == "false" ]; then
	option_para="${option_para} -f $filter_ctg"
fi

echo "### run salsa run_pipeline.py" >>$script_log
echo "\
python $salsa_path/run_pipeline.py -a ${index_dir}/${asm_file} \
	-l ${index_dir}/${asm_file}.fai -i $iter \
	-b $bed -c $minilen -o $outdir $option_para -m yes -p yes" >>$script_log

python $salsa_path/run_pipeline.py -a ${index_dir}/${asm_file} \
	-l ${index_dir}/${asm_file}.fai -i $iter \
	-b $bed -c $minilen -o $outdir $option_para -m yes -p yes

echo "" >>$script_log

echo "### softlink the final scaffolds" >>$script_log
echo "ln -s -f ${outdir}/scaffolds_FINAL.fasta ${outdir}/../${bam_prefix}.hicsalsa.fasta" >>$script_log
ln -s -f ${outdir}/scaffolds_FINAL.fasta ${outdir}/../${bam_prefix}.hicsalsa.fasta
echo "ln -s -f ${outdir}/scaffolds_FINAL.agp ${outdir}/../${bam_prefix}.hicsalsa.agp" >>$script_log
ln -s -f ${outdir}/scaffolds_FINAL.agp ${outdir}/../${bam_prefix}.hicsalsa.agp

echo "### Salsa is ready." >>$script_log


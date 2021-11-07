#!/bin/bash

##########################################
# modified ARIMA GENOMICS MAPPING PIPELINE
##########################################
# the script maps paired-end HiC reads to the assembly

##########################################
# Commands #
##########################################

ASMPATH=$1  # path to the assembly fasta
ASMFILE=$2  # file name fo the assembly fasta
PREFIX=$3	# prefix for output files
BWA_OPTS=$4
FASTQ1=$5 # read 1
FASTQ2=$6 # read 2
IDXWD=$7    # index path
REF=${IDXWD}/${ASMFILE}	# ref fasta file with indexes
OUTDIR=$8   # aln path
FINAL_BAM=$9 # final alignment BAM
utils_path="${10}" # path to util scripts
THREADS="${11}"  # threads
SCRIPT_LOG="${12}"  # script log file

BWA='bwa'
SAMTOOLS='samtools'
FILTER="$utils_path/filter_five_end.pl"
COMBINER="$utils_path/two_read_bam_combiner.pl"
DEDUP="$utils_path/dedup.sh"

RAW_DIR="$OUTDIR/raw"
FILT_DIR="$OUTDIR/filtered"
COMB_DIR="$OUTDIR/combined"
DEDUP_DIR="$OUTDIR/dedup"

if [ ! -d $RAW_DIR ]; then
	mkdir -p $RAW_DIR
fi

if [ ! -d $FILT_DIR ]; then
	mkdir -p $FILT_DIR
fi

if [ ! -d $COMB_DIR ]; then
	mkdir -p $COMB_DIR
fi

if [ ! -d $DEDUP_DIR ]; then
	mkdir -p $DEDUP_DIR
fi

echo "### Step 1.A: FASTQ to BAM (1st)" >>$SCRIPT_LOG
echo "\
$BWA mem -t$THREADS $BWA_OPTS $REF $FASTQ1 | $SAMTOOLS view -Sb - > $RAW_DIR/${PREFIX}_1.bam" >>$SCRIPT_LOG
$BWA mem -t$THREADS $BWA_OPTS $REF $FASTQ1 | $SAMTOOLS view -Sb - > $RAW_DIR/${PREFIX}_1.bam
echo ""

echo "### Step 2.A: Filter 5' end (1st)"
echo "\
$SAMTOOLS view -h $RAW_DIR/${PREFIX}_1.bam | perl $FILTER | $SAMTOOLS view -@$THREADS -Sb - > $FILT_DIR/${PREFIX}_1.bam" >>$SCRIPT_LOG
$SAMTOOLS view -h $RAW_DIR/${PREFIX}_1.bam | perl $FILTER | $SAMTOOLS view -@$THREADS -Sb - > $FILT_DIR/${PREFIX}_1.bam &&
echo "" || exit -1

echo "### Remove intermediate file" >>$SCRIPT_LOG
echo "rm $RAW_DIR/${PREFIX}_1.bam" >>$SCRIPT_LOG
rm $RAW_DIR/${PREFIX}_1.bam
echo ""  >>$SCRIPT_LOG

echo "### Step 1.B: FASTQ to BAM (2nd)" >>$SCRIPT_LOG
echo "$BWA mem -t$THREADS $BWA_OPTS $REF $FASTQ2 | $SAMTOOLS view -Sb - > $RAW_DIR/${PREFIX}_2.bam" >>$SCRIPT_LOG
$BWA mem -t$THREADS $BWA_OPTS $REF $FASTQ2 | $SAMTOOLS view -Sb - > $RAW_DIR/${PREFIX}_2.bam
echo "" >>$SCRIPT_LOG

echo "### Step 2.B: Filter 5' end (2nd)" >>$SCRIPT_LOG
echo "\
$SAMTOOLS view -h $RAW_DIR/${PREFIX}_2.bam | perl $FILTER | $SAMTOOLS view -@$THREADS -Sb - > $FILT_DIR/${PREFIX}_2.bam" >>$SCRIPT_LOG
$SAMTOOLS view -h $RAW_DIR/${PREFIX}_2.bam | perl $FILTER | $SAMTOOLS view -@$THREADS -Sb - > $FILT_DIR/${PREFIX}_2.bam &&
echo "" || exit -1

echo "### Remove intermediate file" >>$SCRIPT_LOG
echo "rm $RAW_DIR/${PREFIX}_2.bam" >>$SCRIPT_LOG
rm $RAW_DIR/${PREFIX}_2.bam
echo >>$SCRIPT_LOG

echo "### Step 3.A: Filter Combiner" >>$SCRIPT_LOG
echo "\
perl $COMBINER $FILT_DIR/${PREFIX}_1.bam $FILT_DIR/${PREFIX}_2.bam | $SAMTOOLS view -@$THREADS -Sb > $COMB_DIR/$PREFIX.bam" >>$SCRIPT_LOG
perl $COMBINER $FILT_DIR/${PREFIX}_1.bam $FILT_DIR/${PREFIX}_2.bam | $SAMTOOLS view -@$THREADS -Sb > $COMB_DIR/$PREFIX.bam
echo "" >>$SCRIPT_LOG

echo "### Start to dedup" >>$SCRIPT_LOG
echo "$DEDUP $COMB_DIR/$PREFIX.bam $THREADS $DEDUP_DIR $FINAL_BAM" >>$SCRIPT_LOG
$DEDUP $COMB_DIR/$PREFIX.bam $THREADS $DEDUP_DIR $FINAL_BAM

echo "#### Finished Mapping!" >>$SCRIPT_LOG
echo "" >>$SCRIPT_LOG


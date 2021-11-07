#!/bin/bash
seqtk subseq /bulk/liu3zhen/research/projects/Chang7-2/main_sl/1_hifiasm/1b_asm2_ok/Ch72asm2.p_ctg.fasta \
	mt_list > mt.fasta
# reads
gunzip -c /bulk/liu3zhen/LiuRawData/collaboration/CAAS/Hi-C/Cleandata/1884/1884_R1.fq.gz \
	| head -n 100000 > r1.fastq
gunzip -c /bulk/liu3zhen/LiuRawData/collaboration/CAAS/Hi-C/Cleandata/1884/1884_R2.fq.gz \
	| head -n 100000 > r2.fastq


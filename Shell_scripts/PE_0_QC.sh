#!/usr/bin/env bash

#This script runs fastQC on all your raw fastq files and outputs them in the fastQC directory
#The script is written for paired end reads. If your reads are single ended you need to use SE_0_QC.sh

#read in variables
source common_variables.sh

#run fastQC
for filename in $Totals_filenames
do
fastqc $fastq_dir/${filename}_R1.fastq.gz --outdir=$fastqc_dir &
fastqc $fastq_dir/${filename}_R2.fastq.gz --outdir=$fastqc_dir &
done
wait

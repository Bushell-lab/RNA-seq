#!/usr/bin/env bash

###filenames
#These are the filenames for all Total RNA-seq samples (without the <.fastq> or any alternative extension)

Totals_filenames='Totals_1 Totals_2 Totals_3'

###set the sumber of threads available to use
#It is recomended to use one or two less than what is available and also consider whether any else is being run at the same time
#Some of the packages used do not support multi-threading and so the for loops run in parallel, so that all files are run at the same time. For this reason do not run on more files than the number of cores available to use
threadN=16

###adaptors
#Totals adaptors
Totals_adaptor='AGATCGGAAGAG' #this is the adaptor used in the LEXOGEN CORALL Total RNA-Seq Library Prep Kit

###paths
parent_dir='/Path/to/data' #This is the path to the parent directory that contains all the data and where all the processed data will be saved

#The following directories are where all the processed data will be saved. These all need to be created prior to starting the analysis

#set the directory where the raw bcl data is. the directory that contains the raw sequencing data in bcl format. This is what you get from a sequencing run and needs to be demulitplexed to write the <.fastq> files.
#If you have more than one bcl directory (you will get one for each sequencing run), then hash one out and write a new one below, each time you re-run the demultiplex.sh script script, so that this acts as a log for all the bcl directories associated with this project
bcl_dir='Path/to/bcl/data'

fastq_dir=${parent_dir}/fastq_files
fastqc_dir=${parent_dir}/fastQC_files
SAM_dir=${parent_dir}/SAM_files
BAM_dir=${parent_dir}/BAM_files
log_dir=${parent_dir}/logs

STAR_dir=${parent_dir}/STAR
rsem_dir=${parent_dir}/rsem

#The following directories are where all the csv files that are used as input into R will be saved
analysis_dir=${parent_dir}/Analysis

most_abundant_transcripts_dir=${analysis_dir}/most_abundant_transcripts
DESeq2_dir=${analysis_dir}/DESeq2_output
reads_summary_dir=${analysis_dir}/reads_summary
fgsea_dir=${analysis_dir}/fgsea

#The following directories are where all the plots generated in R will be saved
plots_dir=${parent_dir}/plots

DE_analysis_dir=${plots_dir}/DE_analysis
PCA_dir=${plots_dir}/PCAs
Interactive_scatters_dir=${plots_dir}/Interactive_scatters
fgsea_plots_dir=${plots_dir}/fgsea
fgsea_scatters_dir=${plots_dir}/fgsea/scatters
fgsea_interactive_scatters_dir=${plots_dir}/fgsea/Interactive_scatters

#Fastas
fasta_dir='/Path/to/FASTAs'

pc_fasta=${fasta_dir}/GENCODE/v38/filtered/gencode.v38.pc_transcripts_filtered.fa
rsem_index=${fasta_dir}/GENCODE/v38/filtered/rsem_bowtie2_index/gencode.v38.pc_transcripts_filtered
STAR_index=${fasta_dir}/GENCODE/v38/original/STAR_index
STAR_GTF=${fasta_dir}/GENCODE/v38/original/gencode.v38.annotation.gtf

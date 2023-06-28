#!/usr/bin/env bash

#-s specifies strand-specific read counting. 0 for unstranded reads, 1 for stranded reads and 2 for reversely stranded reads. This depends on the library used in the sequencing protocol.
#-p species that fragments (or templates) will be counted instead of reads. This is only applicable for paired-end reads.

#read in variables
source common_variables.sh

#count gene alignments
for filename in $Totals_filenames
do
featureCounts -a $STAR_GTF -s 1 -o $fCounts_dir/${filename}_fCounts_output.txt $BAM_dir/${filename}_genome_sorted.bam -T $threadN -p 2> $log_dir/${filename}_fCounts_log.txt &
done
wait

for filename in $Totals_filenames
do
cut -f1,7-8 $fCounts_dir/${filename}_fCounts_output.txt > $fCounts_dir/${filename}_fCounts_output_matrix.txt
done

#!/bin/bash

#SAMTOOLS latest directory
export PATH=$PATH:/usr/local/bin/bin
#Feature counts version 2.0.3


#featureCounts for qualtrim deduplicated (150bp)
#input of BAM files for duplicated aligned reads are in GRAMPUS
#With multimapping -M
#-O assign reads to all overlapping features... (TE can occur in multiple places)
#--fraction multimapped reads are fractioned


cd /media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/03_STAR

featureCounts \
 -a /media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/00_META_DATA/Potato_DM_v6.1_reference/DM_1-3_516_R44_potato_genome_assembly.v6.1.hm.out.gtf \
 -t similarity \
 -T 10 \
 -M \
 -O \
 -g Target \
 --fraction \
 -p \
 -o ../05_FEATURECOUNTS/RepetitiveElements_multimapped_counts.txt \
 --countReadPairs \
 *Run2_Aligned.sortedByCoord.out.bam
 

exit

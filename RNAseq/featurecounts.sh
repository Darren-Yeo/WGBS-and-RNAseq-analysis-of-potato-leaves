#!/bin/bash

#SAMTOOLS latest directory
export PATH=$PATH:/usr/local/bin/bin
#Feature counts version 2.0.3


#featureCounts for qualtrim deduplicated (150bp)
#input of BAM files for duplicated aligned reads are in GRAMPUS
#Without multimapping 

cd /media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/03_STAR

 featureCounts\
 -t gene\
 -T 6\
 -p\
 --countReadPairs\
 -g ID\
 -a /media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/00_META_DATA/Potato_DM_v6.1_reference/download.20230626.090008/Phytozome/PhytozomeV13/Stuberosum/v6.1/annotation/Stuberosum_686_v6.1.gene_exons.gff3\
 -o ../05_FEATURECOUNTS/UniqueCountsDuplicatedHighStress.txt\
  *Run2_Aligned.sortedByCoord.out.bam

exit

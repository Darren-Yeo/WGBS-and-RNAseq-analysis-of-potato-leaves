#!/bin/bash
###run TEtranscripts in conda env

path_to_bam=/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/03_STAR/

##Annabelle
TEtranscripts -t ${path_to_bam}RNA4H_Run2_Aligned.sortedByCoord.out.bam ${path_to_bam}RNA5H_Run2_Aligned.sortedByCoord.out.bam ${path_to_bam}RNA15H_Run2_Aligned.sortedByCoord.out.bam\
 -c ${path_to_bam}RNA4K_Run2_Aligned.sortedByCoord.out.bam ${path_to_bam}RNA5K_Run2_Aligned.sortedByCoord.out.bam ${path_to_bam}RNA15K_Run2_Aligned.sortedByCoord.out.bam\
 --GTF /media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/00_META_DATA/Potato_DM_v6.1_reference/download.20230626.090008/Phytozome/PhytozomeV13/Stuberosum/v6.1/annotation/Stuberosum_686_v6.1.gene_exons.gtf\
 --TE /media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/10_TE_script/DM_1-3_516_R44_potato_genome_assembly.v6.1.hm.out.gtf\
 --sortByPos\
 --project AA_TEtranscripts_out\
 --outdir ./Results_AA/

#camel
TEtranscripts -t ${path_to_bam}RNA2H_Run2_Aligned.sortedByCoord.out.bam ${path_to_bam}RNA3H_Run2_Aligned.sortedByCoord.out.bam ${path_to_bam}RNA11H_Run2_Aligned.sortedByCoord.out.bam\
 -c ${path_to_bam}RNA2K_Run2_Aligned.sortedByCoord.out.bam ${path_to_bam}RNA3K_Run2_Aligned.sortedByCoord.out.bam ${path_to_bam}RNA11K_Run2_Aligned.sortedByCoord.out.bam\
 --GTF /media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/00_META_DATA/Potato_DM_v6.1_reference/download.20230626.090008/Phytozome/PhytozomeV13/Stuberosum/v6.1/annotation/Stuberosum_686_v6.1.gene_exons.gtf\
 --TE /media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/10_TE_script/DM_1-3_516_R44_potato_genome_assembly.v6.1.hm.out.gtf\
 --sortByPos\
 --project CA_TEtranscripts_out\
 --outdir ./Results_CA/



exit
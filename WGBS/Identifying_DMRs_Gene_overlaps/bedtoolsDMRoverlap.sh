#!/bin/bash
#bedtools 2.29

##CG
CG=/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/WGBS_HIGH_HEAT/ANALYSIS/CpG_DMRs/FISHER

##CHG
CHG=/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/CHG_DMRs/FISHER

##CHH
CHH=/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/CHH_DMRs/FISHER

path_to_database=/media/rna/Epipotato16TB/00_META_DATA

path_before_After_HS=/BEDTOOLSHeatStressBeforeAfter/

path_Controls=/BEDTOOLS_Controls/

##CG Control
#bedtools intersect -a ${CG}/AA28_CvsCA28_C.bed -b ${path_to_database}/Downstream2kb_genes_only.bed -wa -wb -f 0.8 > ${CG}/${path_Controls}/CG_DownstreamBEDTOOLS.bed
#bedtools intersect -a ${CG}/AA28_CvsCA28_C.bed -b ${path_to_database}/Upstream2kb_genes_only.bed -wa -wb -f 0.8 > ${CG}/${path_Controls}/CG_UpstreamBEDTOOLS.bed
#bedtools intersect -a ${CG}/AA28_CvsCA28_C.bed -b ${path_to_database}/Stuberosum_686_v6.1.gene_only.gff3 -wa -wb -f 0.8 > ${CG}/${path_Controls}/CG_GenebodyBEDTOOLS.bed
#
#
##CHG Control
#bedtools intersect -a ${CHG}/CHG_AA28_CvsCA28_C.bed -b ${path_to_database}/Downstream2kb_genes_only.bed -wa -wb -f 0.8 > ${CHG}/${path_Controls}/CHG_DownstreamBEDTOOLS.bed
#bedtools intersect -a ${CHG}/CHG_AA28_CvsCA28_C.bed -b ${path_to_database}/Upstream2kb_genes_only.bed -wa -wb -f 0.8 > ${CHG}/${path_Controls}/CHG_UpstreamBEDTOOLS.bed
#bedtools intersect -a ${CHG}/CHG_AA28_CvsCA28_C.bed -b ${path_to_database}/Stuberosum_686_v6.1.gene_only.gff3 -wa -wb -f 0.8 > ${CHG}/${path_Controls}/CHG_GenebodyBEDTOOLS.bed
#
#
##CHH Control
#bedtools intersect -a ${CHH}/CHH_AA28_CvsCA28_C.bed -b ${path_to_database}/Downstream2kb_genes_only.bed -wa -wb -f 0.8 > ${CHH}/${path_Controls}/CHH_DownstreamBEDTOOLS.bed
#bedtools intersect -a ${CHH}/CHH_AA28_CvsCA28_C.bed -b ${path_to_database}/Upstream2kb_genes_only.bed -wa -wb -f 0.8 > ${CHH}/${path_Controls}/CHH_UpstreamBEDTOOLS.bed
#bedtools intersect -a ${CHH}/CHH_AA28_CvsCA28_C.bed -b ${path_to_database}/Stuberosum_686_v6.1.gene_only.gff3 -wa -wb -f 0.8 > ${CHH}/${path_Controls}/CHH_GenebodyBEDTOOLS.bed


##### AA before vs AA  after ##
##CG
#bedtools intersect -a ${CG}/AA28vsAA42_H.bed -b ${path_to_database}/Downstream2kb_genes_only.bed -wa -wb -f 0.8 > ${CG}/${path_before_After_HS}/CG_AA_DownstreamBEDTOOLS.bed
#bedtools intersect -a ${CG}/AA28vsAA42_H.bed -b ${path_to_database}/Upstream2kb_genes_only.bed -wa -wb -f 0.8 > ${CG}/${path_before_After_HS}/CG_AA_UpstreamBEDTOOLS.bed
#bedtools intersect -a ${CG}/AA28vsAA42_H.bed -b ${path_to_database}/Stuberosum_686_v6.1.gene_only.gff3 -wa -wb -f 0.8 > ${CG}/${path_before_After_HS}/CG_AA_GenebodyBEDTOOLS.bed
#
##CHG
#bedtools intersect -a ${CHG}/CHG_AA28vsAA42_H.bed -b ${path_to_database}/Downstream2kb_genes_only.bed -wa -wb -f 0.8 > ${CHG}/${path_before_After_HS}/CHG_AA_DownstreamBEDTOOLS.bed
#bedtools intersect -a ${CHG}/CHG_AA28vsAA42_H.bed -b ${path_to_database}/Upstream2kb_genes_only.bed -wa -wb -f 0.8 > ${CHG}/${path_before_After_HS}/CHG_AA_UpstreamBEDTOOLS.bed
#bedtools intersect -a ${CHG}/CHG_AA28vsAA42_H.bed -b ${path_to_database}/Stuberosum_686_v6.1.gene_only.gff3 -wa -wb -f 0.8 > ${CHG}/${path_before_After_HS}/CHG_AA_GenebodyBEDTOOLS.bed
#
#
##CHH
#bedtools intersect -a ${CHH}/CHH_AA28vsAA42_H.bed -b ${path_to_database}/Downstream2kb_genes_only.bed -wa -wb -f 0.8 > ${CHH}/${path_before_After_HS}/CHH_AA_DownstreamBEDTOOLS.bed
#bedtools intersect -a ${CHH}/CHH_AA28vsAA42_H.bed -b ${path_to_database}/Upstream2kb_genes_only.bed -wa -wb -f 0.8 > ${CHH}/${path_before_After_HS}/CHH_AA_UpstreamBEDTOOLS.bed
#bedtools intersect -a ${CHH}/CHH_AA28vsAA42_H.bed -b ${path_to_database}/Stuberosum_686_v6.1.gene_only.gff3 -wa -wb -f 0.8 > ${CHH}/${path_before_After_HS}/CHH_AA_GenebodyBEDTOOLS.bed
#
#
##### CA before vs CA  after ##

#CG
bedtools intersect -a ${CG}/CA28vsCA42_H.bed -b ${path_to_database}/Downstream2kb_genes_only.bed -wa -wb -f 0.8 > ${CG}/${path_before_After_HS}/CG_CA_DownstreamBEDTOOLS.bed
bedtools intersect -a ${CG}/CA28vsCA42_H.bed -b ${path_to_database}/Upstream2kb_genes_only.bed -wa -wb -f 0.8 > ${CG}/${path_before_After_HS}/CG_CA_UpstreamBEDTOOLS.bed
bedtools intersect -a ${CG}/CA28vsCA42_H.bed -b ${path_to_database}/Stuberosum_686_v6.1.gene_only.gff3 -wa -wb -f 0.8 > ${CG}/${path_before_After_HS}/CG_CA_GenebodyBEDTOOLS.bed

#CHG
bedtools intersect -a ${CHG}/CHG_CA28vsCA42_H.bed -b ${path_to_database}/Downstream2kb_genes_only.bed -wa -wb -f 0.8 > ${CHG}/${path_before_After_HS}/CHG_CA_DownstreamBEDTOOLS.bed
bedtools intersect -a ${CHG}/CHG_CA28vsCA42_H.bed -b ${path_to_database}/Upstream2kb_genes_only.bed -wa -wb -f 0.8 > ${CHG}/${path_before_After_HS}/CHG_CA_UpstreamBEDTOOLS.bed
bedtools intersect -a ${CHG}/CHG_CA28vsCA42_H.bed -b ${path_to_database}/Stuberosum_686_v6.1.gene_only.gff3 -wa -wb -f 0.8 > ${CHG}/${path_before_After_HS}/CHG_CA_GenebodyBEDTOOLS.bed

#CHH
bedtools intersect -a ${CHH}/CHH_CA28vsCA42_H.bed -b ${path_to_database}/Downstream2kb_genes_only.bed -wa -wb -f 0.8 > ${CHH}/${path_before_After_HS}/CHH_CA_DownstreamBEDTOOLS.bed
bedtools intersect -a ${CHH}/CHH_CA28vsCA42_H.bed -b ${path_to_database}/Upstream2kb_genes_only.bed -wa -wb -f 0.8 > ${CHH}/${path_before_After_HS}/CHH_CA_UpstreamBEDTOOLS.bed
bedtools intersect -a ${CHH}/CHH_CA28vsCA42_H.bed -b ${path_to_database}/Stuberosum_686_v6.1.gene_only.gff3 -wa -wb -f 0.8 > ${CHH}/${path_before_After_HS}/CHH_CA_GenebodyBEDTOOLS.bed





exit
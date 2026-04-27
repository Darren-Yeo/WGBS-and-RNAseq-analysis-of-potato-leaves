#!/bin/bash
### look for correlation between selected DMRs and SNPs,.. only bi allelic snps

#all bi alleic snps, QUAL >=30
path_to_filtered_SNPS=/media/rna/Epipotato16TB/WGS/05_VARIANT_CALL/filtered_output_parallel_NEWQUAL30.vcf.gz

##all DMR DEGs 
path_to_selected_DMRs=/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/SNP_DMRs/BED_files_DMR_DEGs

Save_path=/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/SNP_DMRs/BED_files_DMR_DEG_overlaps

cd ${path_to_selected_DMRs}

for bedfiles in *.bed; do

echo ${bedfiles}

filename="${bedfiles%.bed}"


bedtools intersect -a ${path_to_filtered_SNPS} -b ${bedfiles} -wa -wb > ${Save_path}/${filename}_DMRDEG_overlaps.tsv


done





exit
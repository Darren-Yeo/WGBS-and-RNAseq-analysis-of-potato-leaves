#!/bin/bash
### look for correlation between selected DMRs and SNPs,.. only bi allelic snps

path_to_filtered_SNPS=/media/rna/Epipotato16TB/WGS/05_VARIANT_CALL/filtered_output_parallel_NEWQUAL30.vcf.gz

path_to_selected_DMRs=/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/SNP_DMRs/BED_files_DMR_GENEs

Save_path=/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/SNP_DMRs/BED_files_DMR_GENEs_Overlaps

cd ${path_to_selected_DMRs}

for bedfiles in *.bed; do

filename="${bedfiles%.bed}"
echo ${filename:8}

bedtools intersect -a ${path_to_filtered_SNPS} -b ${bedfiles} -wa -wb > ${Save_path}/${filename:8}_SNP_overlaps.tsv

done

exit
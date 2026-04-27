#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 13 16:05:59 2025

@author: rna
"""

import bsxplorer as bsx
import matplotlib.pyplot as plt

levels_AA = bsx.ChrLevels.from_bismark("/media/rna/Epipotato16TB/Epipotato_methylome/04_BISMARK/METHYLATION_CX_REPORT/Unfiltered/CX_PARENTS_MERGED/methylationDataList_Parents_pooled_AA_C.CX.txt", 
                                    chr_min_length=10**6, 
                                    window_length=500000)

levels_CA = bsx.ChrLevels.from_bismark("/media/rna/Epipotato16TB/Epipotato_methylome/04_BISMARK/METHYLATION_CX_REPORT/Unfiltered/CX_PARENTS_MERGED/methylationDataList_Parents_pooled_CA_C.CX.txt", 
                                    chr_min_length=10**6, 
                                    window_length=500000)



###load single dataset
levels_AA_1 = bsx.ChrLevels.from_bismark("/media/rna/Epipotato16TB/Epipotato_methylome/04_BISMARK/METHYLATION_CX_REPORT/Unfiltered/CHR_only/CHR_NAME_SORTED_Qualtrim_AA_A_1_bismark_bt2_pe.deduplicated.CX_report.txt", 
                                    chr_min_length=10**6, 
                                    window_length=500000)


levels_AA_2 = bsx.ChrLevels.from_bismark("/media/rna/Epipotato16TB/Epipotato_methylome/04_BISMARK/METHYLATION_CX_REPORT/Unfiltered/CHR_only/CHR_NAME_SORTED_Qualtrim_AA_B_1_bismark_bt2_pe.deduplicated.CX_report.txt", 
                                    chr_min_length=10**6, 
                                    window_length=500000)

levels_CA_1 = bsx.ChrLevels.from_bismark("/media/rna/Epipotato16TB/Epipotato_methylome/04_BISMARK/METHYLATION_CX_REPORT/Unfiltered/CHR_only/CHR_NAME_SORTED_Qualtrim_CA_A_1_bismark_bt2_pe.deduplicated.CX_report.txt", 
                                    chr_min_length=10**6, 
                                    window_length=500000)

levels_CA_2 = bsx.ChrLevels.from_bismark("/media/rna/Epipotato16TB/Epipotato_methylome/04_BISMARK/METHYLATION_CX_REPORT/Unfiltered/CHR_only/CHR_NAME_SORTED_Qualtrim_CA_B_1_bismark_bt2_pe.deduplicated.CX_report.txt", 
                                    chr_min_length=10**6, 
                                    window_length=500000)

###save the merged files
levels_AA.filter(context="CG").save_plot_rds("/media/rna/Epipotato16TB/Epipotato_methylome/04_BISMARK/METHYLATION_CX_REPORT/Unfiltered/CX_PARENTS_MERGED/CG_AA_chr_level.rds")
levels_AA.filter(context="CHG").save_plot_rds("/media/rna/Epipotato16TB/Epipotato_methylome/04_BISMARK/METHYLATION_CX_REPORT/Unfiltered/CX_PARENTS_MERGED/CHG_AA_chr_level.rds")
levels_AA.filter(context="CHH").save_plot_rds("/media/rna/Epipotato16TB/Epipotato_methylome/04_BISMARK/METHYLATION_CX_REPORT/Unfiltered/CX_PARENTS_MERGED/CHH_AA_chr_level.rds")

levels_CA.filter(context="CG").save_plot_rds("/media/rna/Epipotato16TB/Epipotato_methylome/04_BISMARK/METHYLATION_CX_REPORT/Unfiltered/CX_PARENTS_MERGED/CG_CA_chr_level.rds")
levels_CA.filter(context="CHG").save_plot_rds("/media/rna/Epipotato16TB/Epipotato_methylome/04_BISMARK/METHYLATION_CX_REPORT/Unfiltered/CX_PARENTS_MERGED/CHG_CA_chr_level.rds")
levels_CA.filter(context="CHH").save_plot_rds("/media/rna/Epipotato16TB/Epipotato_methylome/04_BISMARK/METHYLATION_CX_REPORT/Unfiltered/CX_PARENTS_MERGED/CHH_CA_chr_level.rds")



levels_AA_1.filter(context="CG").save_plot_rds("/media/rna/Epipotato16TB/Epipotato_methylome/04_BISMARK/METHYLATION_CX_REPORT/Unfiltered/CX_PARENTS_SINGLE/CG_AA_1_chr_level.rds")
levels_AA_1.filter(context="CHG").save_plot_rds("/media/rna/Epipotato16TB/Epipotato_methylome/04_BISMARK/METHYLATION_CX_REPORT/Unfiltered/CX_PARENTS_SINGLE/CHG_AA_1_chr_level.rds")
levels_AA_1.filter(context="CHH").save_plot_rds("//media/rna/Epipotato16TB/Epipotato_methylome/04_BISMARK/METHYLATION_CX_REPORT/Unfiltered/CX_PARENTS_SINGLE/CHH_AA_1_chr_level.rds")

levels_CA_1.filter(context="CG").save_plot_rds("/media/rna/Epipotato16TB/Epipotato_methylome/04_BISMARK/METHYLATION_CX_REPORT/Unfiltered/CX_PARENTS_SINGLE/CG_CA_1_chr_level.rds")
levels_CA_1.filter(context="CHG").save_plot_rds("/media/rna/Epipotato16TB/Epipotato_methylome/04_BISMARK/METHYLATION_CX_REPORT/Unfiltered/CX_PARENTS_SINGLE/CHG_CA_1_chr_level.rds")
levels_CA_1.filter(context="CHH").save_plot_rds("/media/rna/Epipotato16TB/Epipotato_methylome/04_BISMARK/METHYLATION_CX_REPORT/Unfiltered/CX_PARENTS_SINGLE/CHH_CA_1_chr_level.rds")






levels.save_plot_rds()

levels_AA.filter(context="CG").draw_mpl(smooth=10)

levels_AA.filter(context="CHG").draw_mpl(smooth=10)

levels_AA.filter(context="CHH").draw_mpl(smooth=10)

levels_CA.filter(context="CG").draw_mpl(smooth=10)

levels_CA.filter(context="CHG").draw_mpl(smooth=10)

levels_CA.filter(context="CHH").draw_mpl(smooth=10)


levels = bsx.ChrLevels.from_bismark("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/CX_Files/CHR_only/CHR_Qualtrim_CA11_1_bismark_bt2_pe.deduplicated.CX_report.txt", 
                                    chr_min_length=10**6, 
                                    window_length=10**6)


levels.draw_mpl(smooth=5)
levels.filter(context="CG").draw_mpl(smooth=5)

levels.draw_mpl(smooth=5)
levels.filter(context="CHG").draw_mpl(smooth=5)

levels.draw_mpl(smooth=5)
levels.filter(context="CHH").draw_mpl(smooth=5)


levels.draw_mpl(smooth=5)
levels.filter(context="CG").draw_mpl(smooth=5)


levels_metagene = bsx.ChrLevels([
      levels.filter(context="CG"),
      levels.filter(context="CHG"),
      levels.filter(context="CHH")])

#### FROM high heat-stress experiment ###


###load single dataset

path="/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/WGBS_HIGH_HEAT/03_BISMARK/CX_FILES/"
path_to_save="/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/WGBS_HIGH_HEAT/03_BISMARK/CX_FILES/Single_Chr_level_RDS/"

levels_AA_1 = bsx.ChrLevels.from_bismark(f"{path}Qualtrim_AA4K_1_bismark_bt2_pe.deduplicated.CX_report.txt", chr_min_length=10**6, window_length=500000)
levels_AA_2 = bsx.ChrLevels.from_bismark(f"{path}Qualtrim_AA5K_1_bismark_bt2_pe.deduplicated.CX_report.txt", chr_min_length=10**6, window_length=500000)
levels_AA_3 = bsx.ChrLevels.from_bismark(f"{path}Qualtrim_AA15K_1_bismark_bt2_pe.deduplicated.CX_report.txt", chr_min_length=10**6, window_length=500000)


levels_CA_1 = bsx.ChrLevels.from_bismark(f"{path}Qualtrim_CA2K_1_bismark_bt2_pe.deduplicated.CX_report.txt", chr_min_length=10**6, window_length=500000)
levels_CA_2 = bsx.ChrLevels.from_bismark(f"{path}Qualtrim_CA3_1_bismark_bt2_pe.deduplicated.CX_report.txt", chr_min_length=10**6, window_length=500000)
levels_CA_3 = bsx.ChrLevels.from_bismark(f"{path}Qualtrim_CA11_1_bismark_bt2_pe.deduplicated.CX_report.txt", chr_min_length=10**6, window_length=500000)


levels_CA_1.filter(context="CG").draw_mpl(smooth=10)


levels_AA_1.filter(context="CG").save_plot_rds(f"{path_to_save}CG_levels_AA_1.rds")
levels_AA_1.filter(context="CHG").save_plot_rds(f"{path_to_save}CHG_levels_AA_1.rds")
levels_AA_1.filter(context="CHH").save_plot_rds(f"{path_to_save}CHH_levels_AA_1.rds")

levels_CA_1.filter(context="CG").save_plot_rds(f"{path_to_save}CG_levels_CA_1.rds")
levels_CA_1.filter(context="CHG").save_plot_rds(f"{path_to_save}CHG_levels_CA_1.rds")
levels_CA_1.filter(context="CHH").save_plot_rds(f"{path_to_save}CHH_levels_CA_1.rds")
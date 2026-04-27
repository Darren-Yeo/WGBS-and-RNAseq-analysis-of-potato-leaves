#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 30 11:17:59 2025

@author: rna
"""

import bsxplorer as bsx
import matplotlib.pyplot as plt


genome = bsx.Genome.from_custom(
    file="/media/rna/Epipotato16TB/00_META_DATA/v6.1_DM_Phytozome/DM_1-3_516_R44_potato_genome_assembly.v6.1.hm.out.gff",
    chr_col=0,
    start_col=3,
    end_col=4,
    strand_col=6,
    comment_char="#",
    has_header=False
)

genome = bsx.Genome.from_custom(
    file="/media/rna/Epipotato16TB/00_META_DATA/v6.1_DM_Phytozome/DM_1-3_516_R44_potato_genome_assembly.v6.1.hm.out.gtf",
    chr_col=0,
    start_col=3,
    end_col=4,
    strand_col=6,
    comment_char="#",
    has_header=False
)


genome.all(min_length=0, flank_length=2000)



#define base path
path="/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/CX_Files/CHR_only/"

##path_to_save_results
path_to_save="/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/MethylationPopStruct_1st2nd_Exp/BISMARKPLOT/merged/"
###FULL dataset
file_names={
"AA153W37C" : "CHR_Qualtrim_AA153W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt",
"AA15K" : "CHR_Qualtrim_AA15K_1_bismark_bt2_pe.deduplicated.CX_report.txt",
"AA43W37C" : "CHR_Qualtrim_AA43W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt",
"AA4K" : "CHR_Qualtrim_AA4K_1_bismark_bt2_pe.deduplicated.CX_report.txt",
"AA53W37C" : "CHR_Qualtrim_AA53W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt",
"AA5K" : "CHR_Qualtrim_AA5K_1_bismark_bt2_pe.deduplicated.CX_report.txt",
"CA11" : "CHR_Qualtrim_CA11_1_bismark_bt2_pe.deduplicated.CX_report.txt",
"CA113W37C" : "CHR_Qualtrim_CA113W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt",
"CA23W37C" : "CHR_Qualtrim_CA23W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt",
"CA2K" : "CHR_Qualtrim_CA2K_1_bismark_bt2_pe.deduplicated.CX_report.txt",
"CA3" : "CHR_Qualtrim_CA3_1_bismark_bt2_pe.deduplicated.CX_report.txt",
"CA33W37C" : "CHR_Qualtrim_CA33W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt"
    }

Metagene_list={}

for key,file in file_names.items():
    metagene = bsx.Metagene.from_bismark(
        f"{path}{file}",
        genome=genome.all(min_length=0, flank_length=2000),
        up_windows=100, body_windows=200, down_windows=100
        )
    
    Metagene_list[key]=metagene
    
    print(f"{file} has been loaded and stored with key {key}")


###merge AA samples togehter for control and heat-stress
metagene_AA_merge_list=bsx.MetageneFiles(samples=[Metagene_list["AA15K"], Metagene_list["AA4K"], Metagene_list["AA5K"]])
metagene_AA_merge = metagene_AA_merge_list.merge()

metagene_AA_H37C_merge_list=bsx.MetageneFiles(samples=[Metagene_list["AA153W37C"], Metagene_list["AA43W37C"], Metagene_list["AA53W37C"]])
metagene_AA_H37C_merge=metagene_AA_H37C_merge_list.merge()

#merge CA samples together for control and heat-stress
metagene_CA_merge_list=bsx.MetageneFiles(samples=[Metagene_list["CA11"], Metagene_list["CA2K"], Metagene_list["CA3"]])
metagene_CA_merge=metagene_CA_merge_list.merge()

metagene_CA_H37C_merge_list=bsx.MetageneFiles(samples=[Metagene_list["CA113W37C"], Metagene_list["CA23W37C"], Metagene_list["CA33W37C"]])
metagene_CA_H37C_merge=metagene_CA_H37C_merge_list.merge()


metagene_AA_merge.save_tsv("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/MethylationPopStruct_1st2nd_Exp/BISMARKPLOT/merged/TE/Metagene_TE_AA.tsv")
metagene_AA_H37C_merge.save_tsv("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/MethylationPopStruct_1st2nd_Exp/BISMARKPLOT/merged/TE/Metagene_TE_AA37C.tsv")
metagene_CA_merge.save_tsv("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/MethylationPopStruct_1st2nd_Exp/BISMARKPLOT/merged/TE/Metagene_TE_CA.tsv")
metagene_CA_H37C_merge.save_tsv("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/MethylationPopStruct_1st2nd_Exp/BISMARKPLOT/merged/TE/Metagene_TE_CA37C.tsv")


cat_metagene = bsx.MetageneFiles([
    metagene_AA_merge.filter(context="CHH"),
    metagene_AA_H37C_merge.filter(context="CHH"),
    metagene_CA_merge.filter(context="CHH"),
    metagene_CA_H37C_merge.filter(context="CHH")
], labels=["AA_28dap_C", "AA_H_42dap", "CA_28dap_C", "CA_H_42dap"])

cat_metagene.line_plot(stat="wmean").draw_mpl(minor_labels=Minor_labels)

#plot metagene for all context for control samples, only one replicate
Minor_labels = ["-2000kb","Repeat body", "+2000kb"]
contexts = ["CG", "CHG", "CHH"]

for ctx in contexts:
    cat_metagene = bsx.MetageneFiles([
        metagene_AA_merge.filter(context=ctx),
        metagene_AA_H37C_merge.filter(context=ctx),
        metagene_CA_merge.filter(context=ctx),
        metagene_CA_H37C_merge.filter(context=ctx)
    ], labels=["AA_28dap_C", "AA_H_42dap", "CA_28dap_C", "CA_H_42dap"])

    # Create the plot for this context
    ax = cat_metagene.line_plot(stat="wmean").draw_mpl(minor_labels=Minor_labels)

    # Get the figure from the Axes
    fig = ax.get_figure()

    # Save the figure with context name in filename
    fig.savefig(f"{path_to_save}Metagene_TE_{ctx}.svg", dpi=300, format='svg')

    # Close the figure to free memory
    plt.close(fig)
    
    
    metagene_AA_merge.filter(context=ctx).save_tsv(f"{path_to_save}Metagene_TE_AA_{ctx}.tsv")
    metagene_AA_H37C_merge.filter(context=ctx).save_tsv(f"{path_to_save}Metagene_TE_AA37C_{ctx}.tsv")
    metagene_CA_merge.filter(context=ctx).save_tsv(f"{path_to_save}Metagene_TE_CA_{ctx}.tsv")
    metagene_CA_H37C_merge.filter(context=ctx).save_tsv(f"{path_to_save}Metagene_TE_CA37C0_{ctx}.tsv")
    
    print(f"Saved Metagene plot for context {ctx}")

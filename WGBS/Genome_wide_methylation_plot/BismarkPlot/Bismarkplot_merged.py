#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 27 10:09:26 2025

@author: rna
"""

import bsxplorer as bsx
import matplotlib.pyplot as plt


#initialize the genome
genome = bsx.Genome.from_gff("/media/rna/Epipotato16TB/00_META_DATA/v6.1_DM_Phytozome/download.20230626.090008/Phytozome/PhytozomeV13/Stuberosum/v6.1/annotation/Stuberosum_686_v6.1.gene_exons.gff3")
genome.gene_body(min_length=0, flank_length=2000)

metagene = bsx.Metagene.from_bismark(
    "/media/rna/Epipotato16TB/Epipotato_methylome/04_BISMARK/METHYLATION_CX_REPORT/Unfiltered/NAME_SORTED_Qualtrim_AA_A_1_bismark_bt2_pe.deduplicated.CX_report.txt",
    genome=genome.gene_body(min_length=0, flank_length=2000),
    up_windows=100, body_windows=200, down_windows=100
)


window_kwargs = dict(up_windows=200, body_windows=400, down_windows=200)
metagene_AA=bsx.Metagene.from_bismark(
    "/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/CX_Files/Qualtrim_AA4K_1_bismark_bt2_pe.deduplicated.CX_report.txt",
    genome=genome.gene_body(min_length=0, flank_length=2000),**window_kwargs)


metagene_CA=bsx.Metagene.from_bismark(
    "/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/CX_Files/Qualtrim_CA2K_1_bismark_bt2_pe.deduplicated.CX_report.txt",
    genome=genome.gene_body(min_length=0, flank_length=2000),**window_kwargs)


cat_metagene = bsx.MetageneFiles([
    metagene_AA.filter(context="CHH"),
    metagene_CA.filter(context="CHH")
], labels=["AA", "CA" ])

tick_labels = ["-2000kb", "TSS", "", "TES", "+2000kb"]
cat_metagene.line_plot().draw_mpl()


######try to merge the replicates together...
#initialize the genome
genome = bsx.Genome.from_gff("/media/rna/Epipotato16TB/00_META_DATA/v6.1_DM_Phytozome/download.20230626.090008/Phytozome/PhytozomeV13/Stuberosum/v6.1/annotation/Stuberosum_686_v6.1.gene_exons.gff3")
genome.gene_body(min_length=0, flank_length=2000)

path = "/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/CX_Files/CHR_only/"
metagene_AA_C1 = bsx.Metagene.from_bismark(
    f"{path}CHR_Qualtrim_AA4K_1_bismark_bt2_pe.deduplicated.CX_report.txt",
    genome=genome.gene_body(min_length=0, flank_length=2000),
    up_windows=100, body_windows=200, down_windows=100
)

metagene_AA_C2 = bsx.Metagene.from_bismark(
    f"{path}CHR_Qualtrim_AA5K_1_bismark_bt2_pe.deduplicated.CX_report.txt",
    genome=genome.gene_body(min_length=0, flank_length=2000),
    up_windows=100, body_windows=200, down_windows=100
)

metagene_AA_C3 = bsx.Metagene.from_bismark(
    f"{path}CHR_Qualtrim_AA15K_1_bismark_bt2_pe.deduplicated.CX_report.txt",
    genome=genome.gene_body(min_length=0, flank_length=2000),
    up_windows=100, body_windows=200, down_windows=100
)

# Step 2: Merge them
metagene_files = bsx.MetageneFiles(samples=[metagene_AA_C1, metagene_AA_C2, metagene_AA_C3])
metagene_AA_merge = metagene_files.merge()

metagene_AA_merge.save_tsv("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/MethylationPopStruct_1st2nd_Exp/BISMARKPLOT/merged/metagene_AA_merge.tsv")
dir(metagene_AA_merge)


cat_metagene = bsx.MetageneFiles([
    metagene_AA_merge.filter(context="CG"),
    metagene_AA_merge.filter(context="CHG"),
    metagene_AA_merge.filter(context="CHH")
] )

cat_metagene.line_plot(stat="wmean").draw_mpl()


####NOW do or all samples

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
        genome=genome.gene_body(min_length=0, flank_length=2000),
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

 metagene_AA_merge.save_tsv(f"{path_to_save}GB/Metagene_AA.tsv")
 metagene_AA_H37C_merge.save_tsv(f"{path_to_save}GB/Metagene_AA37C.tsv")
 metagene_CA_merge.save_tsv(f"{path_to_save}GB/Metagene_CA.tsv")
 metagene_CA_H37C_merge.save_tsv(f"{path_to_save}GB/Metagene_CA37C.tsv")


###now plot
Minor_labels = ["-2000kb","Gene body", "+2000kb"]
contexts = ["CG", "CHG", "CHH"]

 cat_metagene = bsx.MetageneFiles([
     metagene_AA_merge.filter(context="CHG"),
     metagene_AA_H37C_merge.filter(context="CHG"),
     metagene_CA_merge.filter(context="CHG"),
     metagene_CA_H37C_merge.filter(context="CHG")
 ], labels=["AA_28dap_C", "AA_H_42dap", "CA_28dap_C", "CA_H_42dap"])

cat_metagene.line_plot(stat="wmean").draw_mpl(minor_labels=Minor_labels)

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
    fig.savefig(f"{path_to_save}Metagene_{ctx}.svg", dpi=300, format='svg')

    # Close the figure to free memory
    plt.close(fig)
    
    metagene_AA_merge.filter(context=ctx).save_tsv(f"{path_to_save}Metagene_AA_{ctx}.tsv")
    metagene_AA_H37C_merge.filter(context=ctx).save_tsv(f"{path_to_save}Metagene_AA37C_{ctx}.tsv")
    metagene_CA_merge.filter(context=ctx).save_tsv(f"{path_to_save}Metagene_CA_{ctx}.tsv")
    metagene_CA_H37C_merge.filter(context=ctx).save_tsv(f"{path_to_save}Metagene_CA37C0_{ctx}.tsv")

    print(f"Saved Metagene plot for context {ctx}")

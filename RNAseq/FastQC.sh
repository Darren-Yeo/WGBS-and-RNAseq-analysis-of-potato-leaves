#!/bin/bash


path=/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/02_BBduk
Output_path=/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/01_FASTQC_TRIMMED


fastqc -o "$Output_path" "$path"/*.fq.gz -t 11

multiqc -o "$Output_path" "$Output_path" .


exit
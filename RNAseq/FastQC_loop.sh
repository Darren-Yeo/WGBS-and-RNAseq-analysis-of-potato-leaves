#!/bin/bash


path=/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/X208SC24118853-Z03-F001/01.RawData
Output_path=/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/01_FASTQC


for dir in "$path"/*;do

    if [ -d "$dir" ] ; then

        echo "$dir exist, running fastqc"

        fastqc -o "$Output_path" "$dir"/*.fq.gz -t 11

    fi

done

multiqc -o "$Output_path" "$Output_path" .


exit
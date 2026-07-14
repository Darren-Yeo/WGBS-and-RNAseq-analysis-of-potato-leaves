#!/bin/bash


path=/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/02_BBduk

ls $path/*_1.fq.gz | xargs -n 1 basename | sort > $path/samples_1.txt

ls $path/*_2.fq.gz | xargs -n 1 basename | sort > $path/samples_2.txt

ls -d /media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/X208SC24118853-Z03-F001/01.RawData/RNA* | xargs -n 1 basename | sort > $path/Name.txt

paste -d'\t' $path/samples_1.txt $path/samples_2.txt $path/Name.txt > $path/output.txt


exit
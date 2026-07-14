#!/bin/bash
### All promoters are already in the 5' to 3'
##hence --norc
path=/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/06_ANALYSIS/EpiPotato_mRNA_HighStress/AnnabelleVsCamel

cd ${path}

mkdir -p ./FIMO/up
mkdir -p ./FIMO/down

## run for upregulated
echo "Run FIMO Upregulated"
fimo \
 --oc ./FIMO/up \
 --verbosity 1 \
 --bfile ./FIMO/Database_file.txt \
 --thresh 0.05 \
 --qv-thresh \
 --norc \
 ./FIMO/JASPAR2024_CORE_plants_non-redundant_pfms_meme.txt \
 ./AAvsCA_DEGs_2kb_Upregulated.fa &

# run for downregulated
echo "Run FIMO Downregulated"
 fimo \
 --oc ./FIMO/down \
 --verbosity 1 \
 --bfile ./FIMO/Database_file.txt \
 --thresh 0.05 \
 --qv-thresh \
 --norc \
 ./FIMO/JASPAR2024_CORE_plants_non-redundant_pfms_meme.txt \
 ./AAvsCA_DEGs_2kb_Downregulated.fa &

# wait for both background processes
wait
echo "Both FIMO runs finished!"
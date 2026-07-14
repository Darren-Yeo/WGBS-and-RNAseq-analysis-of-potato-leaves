#!/bin/bash
#bbmap 39.01

path=/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/X208SC24118853-Z03-F001/01.RawData
output_path=../../02_BBduk

cd $path

for dir in *;do

    if [ -d "$dir" ] ; then

         fwd=${dir}_1.fq.gz
         rvs=${dir}_2.fq.gz

         echo "Trimming $fwd $rvs"

#Trimming of left end ktrim=r
         echo "Right Trimming ${dir}" 
         /home/rna/Softwares/BBMap_39.01/bbmap/bbduk.sh \
         in=./$dir/${fwd} in2=./$dir/${rvs} \
         out=$output_path/TrimR_${fwd} out2=$output_path/TrimR_${rvs} \
         ref=/home/rna/Softwares/BBMap_39.01/bbmap/resources/adapters.fa \
         ktrim=r k=23 mink=11 hdist=1 threads=11 \
         tpe tbo


#Trimming of left end ktrim=l
         echo "Left Trimming ${dir}"
         /home/rna/Softwares/BBMap_39.01/bbmap/bbduk.sh \
         in=${output_path}/TrimR_${fwd} in2=${output_path}/TrimR_${rvs} \
         out=${output_path}/TrimL_${fwd} out2=${output_path}/TrimL_${rvs} \
         ref=/home/rna/Softwares/BBMap_39.01/bbmap/resources/adapters.fa \
         ktrim=l k=23 mink=11 hdist=1 threads=11 \
         tpe tbo
 
#Base Trimming, Quality filtiering, length filtering
          echo "Quality Trimming ${dir}"
          /home/rna/Softwares/BBMap_39.01/bbmap/bbduk.sh \
          in=${output_path}/TrimL_${fwd} in2=${output_path}/TrimL_${rvs} \
          out=${output_path}/Qual_${fwd} out2=${output_path}/Qual_${rvs} \
          qtrim=rl trimq=30 maq=25 minlen=35 threads=11;


          rm ${output_path}/TrimR_* ${output_path}/TrimL_*
    fi

done
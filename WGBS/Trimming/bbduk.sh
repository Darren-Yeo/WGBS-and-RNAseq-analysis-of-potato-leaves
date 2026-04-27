#!/bin/bash 
#bbmap 39.01
#fastqc v0.11.9

kmerlength=23
minklength=11
hammingdist=1
Threads=20
TrimQualilty=30
minQuality=25
minlength=35

RawData=/path/to/raw/data
BBduk_path=/path/to/trimmed/data



#GO in to methylome raw data for Annabelle
cd $RawData

ls *_1.fq.gz | sort > samples_1.txt

ls *_2.fq.gz | sort > samples_2.txt

paste -d'\t' samples_1.txt samples_2.txt > output.txt

date > $BBduk_path/BbdukTime.txt

cat output.txt | while read -s fwd rvs ;do

echo ${fwd} ${rvs}

#Trimming right end
/media/rna/Epipotato16TB/Softwares/BBMap_39.01/bbmap/bbduk.sh\
 in=./${fwd}\
 in2=./${rvs}\
 out=$BBduk_path/TrimR_${fwd}\
 out2=$BBduk_path/TrimR_${rvs}\
 ref=/media/rna/Epipotato16TB/Softwares/BBMap_39.01/bbmap/resources/adapters.fa\
 ktrim=r\
 k=$kmerlength\
 mink=$minklength\
 hdist=$hammingdist\
 threads=$Threads\
 tpe\
 tbo

#Trimming of left end
/media/rna/Epipotato16TB/Softwares/BBMap_39.01/bbmap/bbduk.sh\
 in=$BBduk_path/TrimR_${fwd}\
 in2=$BBduk_path/TrimR_${rvs}\
 out=$BBduk_path/TrimL_${fwd}\
 out2=$BBduk_path/TrimL_${rvs}\
 ref=/media/rna/Epipotato16TB/Softwares/BBMap_39.01/bbmap/resources/adapters.fa\
 ktrim=l\
 k=$kmerlength\
 mink=$minklength\
 hdist=$hammingdist\
 threads=$Threads\
 tpe\
 tbo
 
#Trimming of the first 10 base pairs
/media/rna/Epipotato16TB/Softwares/BBMap_39.01/bbmap/bbduk.sh\
 in=$BBduk_path/TrimL_${fwd}\
 in2=$BBduk_path/TrimL_${rvs}\
 out=$BBduk_path/TrimFTL_${fwd}\
 out2=$BBduk_path/TrimFTL_${rvs}\
 ftl=10

#Base Trimming, Quality filtiering, length filtering
/media/rna/Epipotato16TB/Softwares/BBMap_39.01/bbmap/bbduk.sh\
 in=$BBduk_path/TrimFTL_${fwd}\
 in2=$BBduk_path/TrimFTL_${rvs}\
 out=$BBduk_path/Qualtrim_${fwd}\
 out2=$BBduk_path/Qualtrim_${rvs}\
 qtrim=rl\
 trimq=$TrimQualilty\
 maq=$minQuality\
 minlen=$minlength\
 threads=$Threads

rm $BBduk_path/TrimR_*\
 $BBduk_path/TrimL_*\
 $BBduk_path/TrimFTL_*

done

exit
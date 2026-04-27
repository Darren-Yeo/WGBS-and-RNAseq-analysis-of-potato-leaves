#!/bin/bash
##Bismark version v0.24.1
#script ran till bismark

#FOR SAMTOOLS version 1.17, installed from samtools github... location in /usr/local/bin/bin/samtools, when uninstalling, use sudo rm (maybe also the extra bin folder)
#export PATH=$PATH:/usr/local/bin/bin




#
##generate Genome folder from the fasta files
#/path/to//bismark_genome_preparation\
# /path/to/bismark/genome/folder\
# --parallel 11
#

##Run Bismark

path=path/to/trimmed/samples

cd ${path}

ls *_1.fq.gz | sort > samples_1.txt

ls *_2.fq.gz | sort > samples_2.txt

paste -d'\t' samples_1.txt samples_2.txt > output.txt

date > /path/to/output/TimeBismarkRun.txt

cat output.txt | while read -s fwd rvs ;do

echo ${fwd} ${rvs}


/path/to//bismark\
 --genome_folder /path/to/bismark/genome/folder\
 -1 ${fwd}\
 -2 ${rvs}\
 --bowtie2\
 --parallel 11\
 -o /path/to/output/
 #--local\

done

date >> /path/to/output/TimeBismarkRun.txt



cd /path/to/output/

date > ./deduplicated/Bismarkdeduplicate.txt

ls *bismark_bt2_pe.bam > samples.txt


#loop all bam files in MORE_BISMARK_DARA
cat samples.txt | while read -s out;do

/path/to//deduplicate_bismark\
 --output_dir ./deduplicated\
 --bam\
 ${out}

done

date >> ../deduplicated/Bismarkdeduplicate.txt



exit
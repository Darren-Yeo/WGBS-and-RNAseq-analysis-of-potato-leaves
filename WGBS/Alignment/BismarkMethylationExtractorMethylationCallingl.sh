#!/bin/bash
#Bismark caller  v0.24.1
###--CX/--CX_context## not used hence only bedgraph and coverage of CpG context. HAVE TO USE IT IN NEXT RUN!!!

#SAMTOOLS latest directory
export PATH=$PATH:/usr/local/bin/bin

path=path/to/samples

cd ${path}

#cat SamplesNamesSorted_second_run.txt | while read -s bam ;do #SamplesNamesSorted.txt ran till, G_365, everything screwed up.. yeah... Continue from there, screwed up again  327_G sample run 327_G again and others

#echo ${bam}

/path/to/bismark_methylation_extractor\
 --paired-end\
 --bedGraph\
 --cytosine_report\
 --CX\
 --ucsc\
 --gzip\
 --comprehensive\
 --mbias_off\
 --parallel 10\
 -o /path/to/output/directory\
 --genome_folder /path/to/bismark/genome/folder\
 --ignore 10\
 --ignore_r2 10\
 *samples


#done

###STILL NEEDS TO RUN SUMMARY###


exit
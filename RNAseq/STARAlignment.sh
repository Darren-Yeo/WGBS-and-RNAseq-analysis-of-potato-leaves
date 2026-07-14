#!/bin/bash
#2.7.9a_2021-06-25
#sjdbOverhang 120, max read length was 121, thus max length -1, always dependent on max length of reads, but doesnt affect result, only computational resources
STAR\
 --runThreadN 11\
 --runMode genomeGenerate\
 --genomeDir /media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/03_STAR/genome_directory\
 --genomeFastaFiles /media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/00_META_DATA/Potato_DM_v6.1_reference/download.20230626.090008/Phytozome/PhytozomeV13/Stuberosum/v6.1/assembly/Stuberosum_686_v6.1.fa\
 --sjdbGTFfile /media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/00_META_DATA/Potato_DM_v6.1_reference/download.20230626.090008/Phytozome/PhytozomeV13/Stuberosum/v6.1/annotation/Stuberosum_686_v6.1.gene_exons.gff3\
 --sjdbOverhang 120\
 --sjdbGTFfeatureExon exon\
 --sjdbGTFtagExonParentTranscript ID\
 --sjdbGTFtagExonParentGene Parent\
 --genomeSAindexNbases 13 ###!!!FOR STAR 2.7.9a  --genomeSAindexNbases 14 is too large for the genome size=639586700, which may cause seg-fault at the mapping step. Re-run genome generation with recommended --genomeSAindexNbases 13

#FOR SAMTOOLS version 1.17, installed from samtools github... location in /usr/local/bin/bin/samtools, when uninstalling, use sudo rm (maybe also the extra bin folder)
export PATH=$PATH:/usr/local/bin/bin

path=/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/02_BBduk

cd $path

cat output.txt | while read -s fwd rvs out ;do

echo ${fwd} ${rvs} ${out}
#--outFilterMismatchNmax 2\
STAR\
	 --runThreadN 11\
	 --outFilterMultimapScoreRange 1\
	 --genomeDir /media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/03_STAR/genome_directory\
	 --readFilesIn ./${fwd} ./${rvs}\
	 --readFilesCommand zcat\
	 --outFilterMultimapNmax 20\
	 --outFilterMismatchNoverReadLmax 0.02\
	 --alignSJoverhangMin 8\
	 --alignSJDBoverhangMin 1\
	 --alignIntronMin 20\
	 --alignIntronMax 30000\
	 --alignMatesGapMax 1000000\
	 --sjdbOverhang 120\
	 --outFileNamePrefix ../03_STAR/${out}_Run1_\
	 --outSAMtype BAM SortedByCoordinate

STAR\
	 --runThreadN 11\
	 --outFilterMultimapScoreRange 1\
	 --genomeDir /media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/03_STAR/genome_directory\
	 --readFilesIn ./${fwd} ./${rvs}\
	 --readFilesCommand zcat\
	 --outFilterMultimapNmax 20\
	 --outFilterMismatchNoverReadLmax 0.02\
	 --alignSJoverhangMin 8\
	 --alignSJDBoverhangMin 1\
	 --alignIntronMin 20\
	 --alignIntronMax 30000\
	 --alignMatesGapMax 1000000\
	 --sjdbOverhang 120\
	 --outFileNamePrefix ../03_STAR/${out}_Run2_\
	 --outSAMtype BAM SortedByCoordinate\
	 --sjdbFileChrStartEnd ../03_STAR/${out}_Run1_SJ.out.tab\



samtools index ../03_STAR/${out}_Run2_Aligned.sortedByCoord.out.bam

rm ../03_STAR/${out}_Run1_Aligned.sortedByCoord.out.bam ../03_STAR/${out}_Run1_Log*

rm -rf ../03_STAR/${out}_Run1__STARtmp*

done

multiqc ../03_STAR/*Run2_Log.final.out* -f -o ../04_STAR

exit


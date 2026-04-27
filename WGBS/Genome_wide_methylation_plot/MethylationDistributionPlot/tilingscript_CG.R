library(DMRcaller)
library(GenomicRanges)
library(tidyverse)
library(genomation)
library(zoo)
path <- "/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/CX_Files/CHR_only/Context_specific/"

methylationDataList_pooled <- GRangesList(
  
  
  "CA_28dap_pooled" = readBismarkPool(c(paste0(path,"CG_CHR_Qualtrim_CA2K_1_bismark_bt2_pe.deduplicated.CX_report.txt"), 
                                        paste0(path,"CG_CHR_Qualtrim_CA3_1_bismark_bt2_pe.deduplicated.CX_report.txt"),
                                        paste0(path,"CG_CHR_Qualtrim_CA11_1_bismark_bt2_pe.deduplicated.CX_report.txt"))),
  
  "CA_42dap_H_pooled" = readBismarkPool(c(paste0(path,"CG_CHR_Qualtrim_CA23W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt"), 
                                          paste0(path,"CG_CHR_Qualtrim_CA33W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt"),
                                          paste0(path,"CG_CHR_Qualtrim_CA113W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt"))),
  
  ###PArents
  "AA_28dap_pooled" = readBismarkPool(c(paste0(path,"CG_CHR_Qualtrim_AA4K_1_bismark_bt2_pe.deduplicated.CX_report.txt"), 
                                        paste0(path,"CG_CHR_Qualtrim_AA5K_1_bismark_bt2_pe.deduplicated.CX_report.txt"),
                                        paste0(path,"CG_CHR_Qualtrim_AA15K_1_bismark_bt2_pe.deduplicated.CX_report.txt"))),
  
  "AA_42dap_H_pooled" = readBismarkPool(c(paste0(path,"CG_CHR_Qualtrim_AA43W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt"), 
                                          paste0(path,"CG_CHR_Qualtrim_AA53W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt"),
                                          paste0(path,"CG_CHR_Qualtrim_AA153W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt")))
  
)

filter_scaffolds <- function(gr) {
  return(subset(gr, !grepl("scaffold",seqnames(gr))))
}

methylationDataList_pooled_chr <-lapply(methylationDataList_pooled, 
                                        filter_scaffolds) 

methylationDataList_pooled_AllSamples_chr_GR<-joinReplicates(methylationDataList_pooled_chr[[1]], methylationDataList_pooled_chr[[2]]) 
for (h in 2:3) { #10 samples hence 10-1
  print(h)
  methylationDataList_pooled_AllSamples_chr_GR<-joinReplicates(methylationDataList_pooled_AllSamples_chr_GR, methylationDataList_pooled_chr[[h+1]])
}

####Define the column names for methylationDataList_Controls_all_single_chr_genotypes
column_names <- rep(names(methylationDataList_pooled_chr), each=2)
extra_names <- rep(c("meth", "cov"),times=4)
column_extra_names <- paste(column_names,extra_names, sep ="_")

colnames(mcols(methylationDataList_pooled_AllSamples_chr_GR))[2:9] <- column_extra_names



####24696564 ranges
meth_gr<- methylationDataList_pooled_AllSamples_chr_GR
min_cov <- 8
min_samples <- 6

HardFilterCoverage<- function(meth_gr,min_cov=8, min_samples){
  
  # Get coverage column names (adjust pattern if needed)
  cov_cols <- grep("_cov$", colnames(mcols(meth_gr)), value = TRUE)
  
  # Make a matrix of coverage values
  cov_mat <- as.matrix(mcols(meth_gr)[, cov_cols])
  
  # Logical matrix: TRUE where coverage ≥ 8
  pass_cov <- cov_mat >= min_cov
  
  # Count how many samples passed per site
  samples_passing <- rowSums(pass_cov)
  
  # Keep sites where at least 6 samples passed
  meth_gr_filtered <- meth_gr[samples_passing >= min_samples]
  
  return(meth_gr_filtered)
}

####19524431
methylationDataList_filtered_pooled_AllSamples_chr_GR<- HardFilterCoverage(methylationDataList_pooled_AllSamples_chr_GR,
                                                                           min_cov=8,
                                                                           min_samples = 3)


# Define the potato genome 6.1
soltub_genome_frame <- data.frame(chr = c("chr01", "chr02", "chr03", "chr04", "chr05", "chr06", 
                                          "chr07", "chr08", "chr09", "chr10", "chr11", "chr12"),
                                  start = rep(0, times = 12),
                                  end = c(88591686-1, 46102915-1, 60707570-1, 69236331-1, 55599697-1, 59091578-1,
                                          57639317-1, 59226000-1, 67600300-1, 61044151-1, 46777387-1, 59670755-1))

chr_lengths <- soltub_genome_frame$end
names(chr_lengths) <- soltub_genome_frame$chr

##define the seqlengths ghr lengths
seqlengths(methylationDataList_filtered_pooled_AllSamples_chr_GR) <- chr_lengths


TileGenomeToSnippets<- function(meth_gr,tilesize=200,last_tile=TRUE){
  
  
  # Tile genome into 200 bp bins
  tile<- tileGenome(seqlengths(meth_gr),
                    tilewidth = tilesize,
                    cut.last.tile.in.chrom = last_tile)
  
  # Step 2: find overlaps
  hits <- findOverlaps(tile, meth_gr)
  
  # Extract indexes
  tile_idx <- queryHits(hits)
  meth_idx <- subjectHits(hits)
  
  return(list(
    hits_GR = hits,
    tile_index = tile_idx,
    meth_index = meth_idx
  ))
  
}

TotalHits_Tile<- TileGenomeToSnippets(methylationDataList_filtered_pooled_AllSamples_chr_GR,
                                      tilesize = 200,
                                      last_tile = TRUE)

# Create a data.frame with overlap info and methylation data for aggregation
library(dplyr)
tile200bp_CX_DF <- data.frame(tile = TotalHits_Tile$tile_index,
                              AA_28dap_pooled_meth = mcols(methylationDataList_filtered_pooled_AllSamples_chr_GR)$AA_28dap_pooled_meth[TotalHits_Tile$meth_index],
                              AA_28dap_pooled_cov = mcols(methylationDataList_filtered_pooled_AllSamples_chr_GR)$AA_28dap_pooled_cov[TotalHits_Tile$meth_index],
                              
                              AA_42dap_H_pooled_meth = mcols(methylationDataList_filtered_pooled_AllSamples_chr_GR)$AA_42dap_H_pooled_meth[TotalHits_Tile$meth_index],
                              AA_42dap_H_pooled_cov = mcols(methylationDataList_filtered_pooled_AllSamples_chr_GR)$AA_42dap_H_pooled_cov[TotalHits_Tile$meth_index],
                              
                              CA_28dap_pooled_meth = mcols(methylationDataList_filtered_pooled_AllSamples_chr_GR)$CA_28dap_pooled_meth[TotalHits_Tile$meth_index],
                              CA_28dap_pooled_cov = mcols(methylationDataList_filtered_pooled_AllSamples_chr_GR)$CA_28dap_pooled_cov[TotalHits_Tile$meth_index],
                              
                              CA_42dap_H_pooled_meth = mcols(methylationDataList_filtered_pooled_AllSamples_chr_GR)$CA_42dap_H_pooled_meth[TotalHits_Tile$meth_index],
                              CA_42dap_H_pooled_cov = mcols(methylationDataList_filtered_pooled_AllSamples_chr_GR)$CA_42dap_H_pooled_cov[TotalHits_Tile$meth_index])

#####JUST TO verify if its binned correctly
region_of_interest <- GRanges(seqnames = "chr01", ranges = IRanges(start = 122801, end = 123000))

methyl_sites_in_region <- subsetByOverlaps(methylationDataList_filtered_pooled_AllSamples_chr_GR, 
                                           region_of_interest)
##use custom bash code to aggregate for wrighted counts by 200bp windows
fwrite(tile200bp_CX_DF, 
       "/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/MethylationPopStruct_1st2nd_Exp/MethylationDistributionPlot/tile200bp_CX_DF.tsv", 
       sep = "\t")

# Aggregate sums by tile
#tile200bp_CX_DF_agg <- tile200bp_CX_DF %>%
#  group_by(tile) %>%
#  summarise(across(ends_with("_meth"), sum, na.rm = TRUE),
#            across(ends_with("_cov"), sum, na.rm = TRUE))


###import the aggregated tiles
library(tidyverse)
tile200bp_CX_DF_agg <- read.delim("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/MethylationPopStruct_1st2nd_Exp/MethylationDistributionPlot/tile200bp_CX_DF_agg.tsv")

tile200bp_CX_DF_agg %>% summarise(
  zero_AA_28dap_cov = sum(AA_28dap_pooled_cov == 0),
  zero_AA_42dap_H_cov = sum(AA_42dap_H_pooled_cov == 0),
  zero_CA_28dap_cov = sum(CA_28dap_pooled_cov == 0),
  zero_CA_42dap_H_cov = sum(CA_42dap_H_pooled_cov == 0)
)


tile200bp_CpG_weightMeth <- tile200bp_CX_DF_agg %>% arrange(tile) %>% transmute(
  tile=tile,
  AA_28dap_C= (AA_28dap_pooled_meth / AA_28dap_pooled_cov) * 100,
  AA_42dap_H = (AA_42dap_H_pooled_meth / AA_42dap_H_pooled_cov) * 100,
  CA_28dap_C = (CA_28dap_pooled_meth / CA_28dap_pooled_cov) * 100,
  CA_42dap_H = (CA_42dap_H_pooled_meth / CA_42dap_H_pooled_cov) * 100
)

tile200bp_CpG_weightMeth %>%
  filter(if_any(everything(), is.na)) %>% nrow()

##Omit all NA values, likely come from at least one of 4 samples have 0 coverage
tile200bp_CpG_weightMeth_na_omit <- na.omit(tile200bp_CpG_weightMeth)

## test for normality
shapiro.test(sample(tile200bp_CpG_weightMeth_na_omit$AA_28dap_C, 5000))
shapiro.test(sample(tile200bp_CpG_weightMeth_na_omit$AA_42dap_H, 5000))
shapiro.test(sample(tile200bp_CpG_weightMeth_na_omit$CA_28dap_C, 5000))
shapiro.test(sample(tile200bp_CpG_weightMeth_na_omit$CA_42dap_H, 5000))

wilcox.test(tile200bp_CpG_weightMeth$AA_28dap_C, tile200bp_CpG_weightMeth$AA_42dap_H, paired = TRUE)
wilcox.test(tile200bp_CpG_weightMeth$CA_28dap_C, tile200bp_CpG_weightMeth$CA_42dap_H, paired = TRUE)




write.table(tile200bp_CpG_weightMeth_na_omit,
            "/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/MethylationPopStruct_1st2nd_Exp/MethylationDistributionPlot/tile200bp_CpG_weightMeth_na_omit.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE)
            
tile200bp_CpG_weightMeth_na_omit_long <- tile200bp_CpG_weightMeth_na_omit %>% 
  pivot_longer(cols = -tile,
               names_to = "Samples",
               values_to = "Weighted_Methylation")
            
            
            
ggplot(tile200bp_CpG_weightMeth_na_omit_long, aes(x = Samples, y = Weighted_Methylation, fill = Samples)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +   # Hide outlier dots
  theme_minimal() +
  labs(
    title = "Weighted Methylation Distribution per Sample",
    x = "Sample",
    y = "Weighted Methylation (%)"
  ) +
  theme(legend.position = "none")


#### Using compute methylation profile

# Create 200 bp genome-wide windows
genome_windows <- tileGenome(
  seqlengths = chr_lengths,
  tilewidth = 200,
  cut.last.tile.in.chrom = TRUE
)

###compute for AA
methylationDataList_pooled_chr_AA <- list()
for (chr in names(chr_lengths)) {
  
  methylationDataList_pooled_chr_AA[[chr]]<- computeMethylationProfile(methylationDataList_pooled_chr$AA_28dap_pooled,
                                                                region = GRanges(seqnames = Rle(chr), ranges = IRanges(1,chr_lengths[chr])),
                                                                windowSize=2000, 
                                                                context = "CG")
}

methylationDataList_pooled_chr_CA <- list()
for (chr in names(chr_lengths)) {
  
  methylationDataList_pooled_chr_CA[[chr]]<- computeMethylationProfile(methylationDataList_pooled_chr$CA_28dap_pooled,
                                                                       region = GRanges(seqnames = Rle(chr), ranges = IRanges(1,chr_lengths[chr])),
                                                                       windowSize=2000, 
                                                                       context = "CG")
}



#combine all together
##AA
methylationDataList_pooled_chr_AA_df<- rbind(as.data.frame(methylationDataList_pooled_chr_AA[["chr01"]]) %>% filter(sumReadsN >=8),
      as.data.frame(methylationDataList_pooled_chr_AA[["chr02"]]) %>% filter(sumReadsN >=8),
      as.data.frame(methylationDataList_pooled_chr_AA[["chr03"]]) %>% filter(sumReadsN >=8),
      as.data.frame(methylationDataList_pooled_chr_AA[["chr04"]]) %>% filter(sumReadsN >=8),
      as.data.frame(methylationDataList_pooled_chr_AA[["chr05"]]) %>% filter(sumReadsN >=8),
      as.data.frame(methylationDataList_pooled_chr_AA[["chr06"]]) %>% filter(sumReadsN >=8),
      as.data.frame(methylationDataList_pooled_chr_AA[["chr07"]]) %>% filter(sumReadsN >=8),
      as.data.frame(methylationDataList_pooled_chr_AA[["chr08"]]) %>% filter(sumReadsN >=8),
      as.data.frame(methylationDataList_pooled_chr_AA[["chr09"]]) %>% filter(sumReadsN >=8),
      as.data.frame(methylationDataList_pooled_chr_AA[["chr10"]]) %>% filter(sumReadsN >=8),
      as.data.frame(methylationDataList_pooled_chr_AA[["chr11"]]) %>% filter(sumReadsN >=8),
      as.data.frame(methylationDataList_pooled_chr_AA[["chr12"]]) %>% filter(sumReadsN >=8))
 
##CA
methylationDataList_pooled_chr_CA_df<- rbind(as.data.frame(methylationDataList_pooled_chr_CA[["chr01"]]) %>% filter(sumReadsN >=8),
                                             as.data.frame(methylationDataList_pooled_chr_CA[["chr02"]]) %>% filter(sumReadsN >=8),
                                             as.data.frame(methylationDataList_pooled_chr_CA[["chr03"]]) %>% filter(sumReadsN >=8),
                                             as.data.frame(methylationDataList_pooled_chr_CA[["chr04"]]) %>% filter(sumReadsN >=8),
                                             as.data.frame(methylationDataList_pooled_chr_CA[["chr05"]]) %>% filter(sumReadsN >=8),
                                             as.data.frame(methylationDataList_pooled_chr_CA[["chr06"]]) %>% filter(sumReadsN >=8),
                                             as.data.frame(methylationDataList_pooled_chr_CA[["chr07"]]) %>% filter(sumReadsN >=8),
                                             as.data.frame(methylationDataList_pooled_chr_CA[["chr08"]]) %>% filter(sumReadsN >=8),
                                             as.data.frame(methylationDataList_pooled_chr_CA[["chr09"]]) %>% filter(sumReadsN >=8),
                                             as.data.frame(methylationDataList_pooled_chr_CA[["chr10"]]) %>% filter(sumReadsN >=8),
                                             as.data.frame(methylationDataList_pooled_chr_CA[["chr11"]]) %>% filter(sumReadsN >=8),
                                             as.data.frame(methylationDataList_pooled_chr_CA[["chr12"]]) %>% filter(sumReadsN >=8))


AA2000bp_GR <- GRanges(
  seqnames = methylationDataList_pooled_chr_AA_df$seqnames,
  ranges = IRanges(start = methylationDataList_pooled_chr_AA_df$start, end = methylationDataList_pooled_chr_AA_df$end),
  strand = methylationDataList_pooled_chr_AA_df$strand,
  score = methylationDataList_pooled_chr_AA_df$Proportion
)

CA2000bp_GR <- GRanges(
  seqnames = methylationDataList_pooled_chr_CA_df$seqnames,
  ranges = IRanges(start = methylationDataList_pooled_chr_CA_df$start, end = methylationDataList_pooled_chr_CA_df$end),
  strand = methylationDataList_pooled_chr_CA_df$strand,
  score = methylationDataList_pooled_chr_CA_df$Proportion
)


library(GenomicRanges)
library(rtracklayer)
export(AA2000bp_GR, "/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/MethylationPopStruct_1st2nd_Exp/MethylationDistributionPlot/AA2000bp_GR.bed", 
       format = "BED")

export(CA2000bp_GR, "/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/MethylationPopStruct_1st2nd_Exp/MethylationDistributionPlot/CA2000bp_GR.bed", 
       format = "BED")

###import the overlaps

SNP_2000bp_meth_overlap_list<- list(
  AA = read.delim("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/MethylationPopStruct_1st2nd_Exp/MethylationDistributionPlot/2000bp_window_SNP/AA_2000bp_SNP_overlaps.tsv", header=FALSE),
  CA = read.delim("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/MethylationPopStruct_1st2nd_Exp/MethylationDistributionPlot/2000bp_window_SNP/CA_2000bp_SNP_overlaps.tsv", header=FALSE))

head(SNP_2000bp_meth_overlap_list$AA)

extract_dosage <- function(gt_string) {
  gt <- strsplit(gt_string, ":")[[1]][1]       # Get only genotype part
  alleles <- unlist(strsplit(gt, "/"))         # Split to alleles
  sum(as.numeric(alleles))                     # Count ALT alleles
}

CleanTables<- function(Table){
  
  ### select only key tables and rename them
  Table_edit<- Table[,c(1:7,10:16)]
  colnames(Table_edit) <- c("chr","POS","V3","REF","ALT","QUAL","V7","AA","CA","Chr_dmr","start_DMR","end_DMR","V15","Meth_Direction")
  
  ##calculate the dosages
  Table_edit$AA_DOS <- sapply(Table_edit$AA, extract_dosage)
  Table_edit$CA_DOS <- sapply(Table_edit$CA, extract_dosage)
  
  #clean table one last time
  Table_edit_clean<- Table_edit %>% select(chr,POS,REF,ALT,AA_DOS,CA_DOS,start_DMR,end_DMR,Meth_Direction)
  
  return(Table_edit_clean)
}


SNP_2000bp_meth_overlap_clean_list<- lapply(SNP_2000bp_meth_overlap_list, CleanTables)

SumMeanDosagesEachDMR<- function(DMR_SNP_TABLE){
  
  ##calculate SUM and Mean for all SNPs within each DMR
  DMR_SNP_TABLE_SUM_MEAN<- DMR_SNP_TABLE %>% 
    group_by(chr, start_DMR, end_DMR, Meth_Direction) %>%
    summarise(
      n_snps = n(),
      sum_AA_DOS = sum(AA_DOS),
      sum_CA_DOS = sum(CA_DOS)) %>%
    ungroup()

  return(DMR_SNP_TABLE_SUM_MEAN)
}


SNP_2000bp_meth_overlap_clean_SumDOS_list<- lapply(SNP_2000bp_meth_overlap_clean_list, SumMeanDosagesEachDMR)


AA_meth_SNP<- SNP_2000bp_meth_overlap_clean_SumDOS_list$AA %>% select(chr,start_DMR,end_DMR,n_snps,sum_AA_DOS, Meth_Direction) %>% rename(start=start_DMR,
                                                                                                            end=end_DMR,
                                                                                                            Methylation=Meth_Direction)
CA_meth_SNP<- SNP_2000bp_meth_overlap_clean_SumDOS_list$CA %>% select(chr,start_DMR,end_DMR,n_snps,sum_CA_DOS, Meth_Direction) %>% rename(start=start_DMR,
                                                                                                           end=end_DMR,
                                                                                                           Methylation=Meth_Direction)

cor.test(AA_meth_SNP$sum_AA_DOS,AA_meth_SNP$Methylation, method="pearson")
cor.test(CA_meth_SNP$sum_CA_DOS,CA_meth_SNP$Methylation, method="pearson")

# AA_meth_SNP: dataframe with sum_AA_DOS and Methylation
ggplot(AA_meth_SNP, aes(x = Methylation, y = sum_AA_DOS)) +
  geom_point(alpha = 0.6, color = "blue") +  # points
  geom_smooth(method = "lm", color = "red") + # linear regression line
  labs(
    title = "Correlation between SNP allele dosage and methylation (AA)",
    x = "Total SNP allele dosage per 2000 bp window",
    y = "Average methylation proportion"
  ) +
  theme_bw()


ggplot(CA_meth_SNP, aes(x = sum_CA_DOS, y = Methylation)) +
  geom_point(alpha = 0.6, color = "blue") +  # points
  geom_smooth(method = "lm", color = "red") + # linear regression line
  labs(
    title = "Correlation between SNP allele dosage and methylation (AA)",
    x = "Total SNP allele dosage per 2000 bp window",
    y = "Average methylation proportion"
  ) +
  theme_bw()



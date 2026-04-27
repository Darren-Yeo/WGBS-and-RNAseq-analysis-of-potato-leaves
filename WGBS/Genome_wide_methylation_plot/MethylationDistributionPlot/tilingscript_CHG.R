library(DMRcaller)
library(GenomicRanges)
library(tidyverse)
library(genomation)
library(zoo)
path <- "/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/CX_Files/CHR_only/Context_specific/"

methylationDataList_CHG_pooled <- GRangesList(
  
  
  "CA_28dap_pooled" = readBismarkPool(c(paste0(path,"CHG_CHR_Qualtrim_CA2K_1_bismark_bt2_pe.deduplicated.CX_report.txt"), 
                                        paste0(path,"CHG_CHR_Qualtrim_CA3_1_bismark_bt2_pe.deduplicated.CX_report.txt"),
                                        paste0(path,"CHG_CHR_Qualtrim_CA11_1_bismark_bt2_pe.deduplicated.CX_report.txt"))),
  
  "CA_42dap_H_pooled" = readBismarkPool(c(paste0(path,"CHG_CHR_Qualtrim_CA23W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt"), 
                                          paste0(path,"CHG_CHR_Qualtrim_CA33W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt"),
                                          paste0(path,"CHG_CHR_Qualtrim_CA113W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt"))),
  
  ###PArents
  "AA_28dap_pooled" = readBismarkPool(c(paste0(path,"CHG_CHR_Qualtrim_AA4K_1_bismark_bt2_pe.deduplicated.CX_report.txt"), 
                                        paste0(path,"CHG_CHR_Qualtrim_AA5K_1_bismark_bt2_pe.deduplicated.CX_report.txt"),
                                        paste0(path,"CHG_CHR_Qualtrim_AA15K_1_bismark_bt2_pe.deduplicated.CX_report.txt"))),
  
  "AA_42dap_H_pooled" = readBismarkPool(c(paste0(path,"CHG_CHR_Qualtrim_AA43W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt"), 
                                          paste0(path,"CHG_CHR_Qualtrim_AA53W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt"),
                                          paste0(path,"CHG_CHR_Qualtrim_AA153W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt")))
  
)

filter_scaffolds <- function(gr) {
  return(subset(gr, !grepl("scaffold",seqnames(gr))))
}

methylationDataList_CHG_pooled_chr <-lapply(methylationDataList_CHG_pooled, 
                                            filter_scaffolds)

methylationDataList_CHG_AllSamples_chr_GR<-joinReplicates(methylationDataList_CHG_pooled[[1]], methylationDataList_CHG_pooled[[2]]) 
for (h in 2:3) { #10 samples hence 10-1
  print(h)
  methylationDataList_CHG_AllSamples_chr_GR<-joinReplicates(methylationDataList_CHG_AllSamples_chr_GR, methylationDataList_CHG_pooled[[h+1]])
}

####Define the column names for methylationDataList_Controls_all_single_chr_genotypes
column_names <- rep(names(methylationDataList_CHG_pooled_chr), each=2)
extra_names <- rep(c("meth", "cov"),times=4)
column_extra_names <- paste(column_names,extra_names, sep ="_")

###28244801
colnames(mcols(methylationDataList_CHG_AllSamples_chr_GR))[2:9] <- column_extra_names

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

###23508496
methylationDataList_filtered_CHG_AllSamples_chr_GR<- HardFilterCoverage(methylationDataList_CHG_AllSamples_chr_GR,
                                                                        min_cov=8,
                                                                        min_samples = 3)

 <- data.frame(chr = c("chr01", "chr02", "chr03", "chr04", "chr05", "chr06", 
                                          "chr07", "chr08", "chr09", "chr10", "chr11", "chr12"),
                                  start = rep(0, times = 12),
                                  end = c(88591686-1, 46102915-1, 60707570-1, 69236331-1, 55599697-1, 59091578-1,
                                          57639317-1, 59226000-1, 67600300-1, 61044151-1, 46777387-1, 59670755-1))

chr_lengths <- soltub_genome_frame$end
names(chr_lengths) <- soltub_genome_frame$chr

##define the seqlengths ghr lengths
seqlengths(methylationDataList_filtered_CHG_AllSamples_chr_GR) <- chr_lengths

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


TotalHits_Tile<- TileGenomeToSnippets(methylationDataList_filtered_CHG_AllSamples_chr_GR,
                                      tilesize = 200,
                                      last_tile = TRUE)


# Create a data.frame with overlap info and methylation data for aggregation
library(dplyr)
tile200bp_CHG_CX_DF <- data.frame(tile = TotalHits_Tile$tile_index,
                                  AA_28dap_pooled_meth = mcols(methylationDataList_filtered_CHG_AllSamples_chr_GR)$AA_28dap_pooled_meth[TotalHits_Tile$meth_index],
                                  AA_28dap_pooled_cov = mcols(methylationDataList_filtered_CHG_AllSamples_chr_GR)$AA_28dap_pooled_cov[TotalHits_Tile$meth_index],
                                  
                                  AA_42dap_H_pooled_meth = mcols(methylationDataList_filtered_CHG_AllSamples_chr_GR)$AA_42dap_H_pooled_meth[TotalHits_Tile$meth_index],
                                  AA_42dap_H_pooled_cov = mcols(methylationDataList_filtered_CHG_AllSamples_chr_GR)$AA_42dap_H_pooled_cov[TotalHits_Tile$meth_index],
                                  
                                  CA_28dap_pooled_meth = mcols(methylationDataList_filtered_CHG_AllSamples_chr_GR)$CA_28dap_pooled_meth[TotalHits_Tile$meth_index],
                                  CA_28dap_pooled_cov = mcols(methylationDataList_filtered_CHG_AllSamples_chr_GR)$CA_28dap_pooled_cov[TotalHits_Tile$meth_index],
                                  
                                  CA_42dap_H_pooled_meth = mcols(methylationDataList_filtered_CHG_AllSamples_chr_GR)$CA_42dap_H_pooled_meth[TotalHits_Tile$meth_index],
                                  CA_42dap_H_pooled_cov = mcols(methylationDataList_filtered_CHG_AllSamples_chr_GR)$CA_42dap_H_pooled_cov[TotalHits_Tile$meth_index])
library(data.table)

##export to aggregate.. weighted meth reads and cov in bash
fwrite(tile200bp_CHG_CX_DF, 
       "/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/MethylationPopStruct_1st2nd_Exp/MethylationDistributionPlot/tile200bp_CHG_CX_DF.tsv", 
       sep = "\t")


###load the aggregate dataset
tile200bp_CHG_CX_DF_AGG <- read.delim("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/MethylationPopStruct_1st2nd_Exp/MethylationDistributionPlot/tile200bp_CHG_CX_DF_AGG.tsv")


##count total number of 0s
tile200bp_CHG_CX_DF_AGG %>% summarise(
  zero_AA_28dap_cov = sum(AA_28dap_pooled_cov == 0),
  zero_AA_42dap_H_cov = sum(AA_42dap_H_pooled_cov == 0),
  zero_CA_28dap_cov = sum(CA_28dap_pooled_cov == 0),
  zero_CA_42dap_H_cov = sum(CA_42dap_H_pooled_cov == 0)
)

tile200bp_CHG_weightMeth <- tile200bp_CHG_CX_DF_AGG %>% 
  arrange(tile) %>% transmute(
    tile=tile,
    AA_28dap_C= (AA_28dap_pooled_meth / AA_28dap_pooled_cov) * 100,
    AA_42dap_H = (AA_42dap_H_pooled_meth / AA_42dap_H_pooled_cov) * 100,
    CA_28dap_C = (CA_28dap_pooled_meth / CA_28dap_pooled_cov) * 100,
    CA_42dap_H = (CA_42dap_H_pooled_meth / CA_42dap_H_pooled_cov) * 100
  )

##determin number of rows with at least one NA
tile200bp_CHG_weightMeth %>%
  filter(if_any(everything(), is.na)) %>% nrow()

##Omit all NA values, likely come from at least one of 4 samples have 0 coverage
tile200bp_CHG_weightMeth_na_omit <- na.omit(tile200bp_CHG_weightMeth)

write.table(tile200bp_CHG_weightMeth_na_omit,
            "/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/MethylationPopStruct_1st2nd_Exp/MethylationDistributionPlot/tile200bp_CHG_weightMeth_na_omit.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE)


tile200bp_CHG_weightMeth_na_omit_long <- tile200bp_CHG_weightMeth_na_omit %>% 
  pivot_longer(cols = -tile,
               names_to = "Samples",
               values_to = "Weighted_Methylation")

ggplot(tile200bp_CHG_weightMeth_na_omit_long, aes(x = Samples, y = Weighted_Methylation, fill = Samples)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +   # Hide outlier dots
  theme_minimal() +
  labs(
    title = "Weighted Methylation Distribution per Sample",
    x = "Sample",
    y = "Weighted Methylation (%)"
  ) +
  theme(legend.position = "none")




###### testing DMRcaller compute methylation profile
regions <-  GRanges(seqnames = Rle("chr01"), ranges = IRanges(0,88591685))
profileCHG <- computeMethylationProfile(methylationDataList_CHG_pooled_chr$CA_28dap_pooled,
                                         regions,
                                         windowSize = 200,
                                         context = "CHG")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           

CalculateMethylationProfile<- function(GRANGES_Object, window=200){
  
  chr_methylationProfile <- list()
  for (CHR in c("chr01","chr02","chr03","chr04","chr05","chr06","chr07","chr08","chr09","chr10","chr11","chr12")) {
    
    print(CHR)
    chr_in_loop<- soltub_genome_frame %>% dplyr::filter(chr==CHR)
    
    regions <- GRanges(seqnames = Rle(CHR), ranges = IRanges(chr_in_loop$start,chr_in_loop$end))
    
    DF<- computeMethylationProfile(GRANGES_Object,
                              regions,
                              windowSize = window,
                              context = "CHG")         
    
    
    chr_methylationProfile[[CHR]] <- DF
  }
  return(chr_methylationProfile)
  
}

CHG_MethylationProfile<- lapply(methylationDataList_CHG_pooled_chr, CalculateMethylationProfile, window=200)

CHG_MethylationProfile_CA<- rbind(mcols(CHG_MethylationProfile$CA_28dap_pooled$chr01),
                                  mcols(CHG_MethylationProfile$CA_28dap_pooled$chr02),
                                  mcols(CHG_MethylationProfile$CA_28dap_pooled$chr03),
                                  mcols(CHG_MethylationProfile$CA_28dap_pooled$chr04),
                                  mcols(CHG_MethylationProfile$CA_28dap_pooled$chr05),
                                  mcols(CHG_MethylationProfile$CA_28dap_pooled$chr06),
                                  mcols(CHG_MethylationProfile$CA_28dap_pooled$chr07),
                                  mcols(CHG_MethylationProfile$CA_28dap_pooled$chr08),
                                  mcols(CHG_MethylationProfile$CA_28dap_pooled$chr09),
                                  mcols(CHG_MethylationProfile$CA_28dap_pooled$chr10),
                                  mcols(CHG_MethylationProfile$CA_28dap_pooled$chr11),
                                  mcols(CHG_MethylationProfile$CA_28dap_pooled$chr12) )

CHG_MethylationProfile_CA_H<- rbind(mcols(CHG_MethylationProfile$CA_42dap_H_pooled$chr01),
                                  mcols(CHG_MethylationProfile$CA_42dap_H_pooled$chr02),
                                  mcols(CHG_MethylationProfile$CA_42dap_H_pooled$chr03),
                                  mcols(CHG_MethylationProfile$CA_42dap_H_pooled$chr04),
                                  mcols(CHG_MethylationProfile$CA_42dap_H_pooled$chr05),
                                  mcols(CHG_MethylationProfile$CA_42dap_H_pooled$chr06),
                                  mcols(CHG_MethylationProfile$CA_42dap_H_pooled$chr07),
                                  mcols(CHG_MethylationProfile$CA_42dap_H_pooled$chr08),
                                  mcols(CHG_MethylationProfile$CA_42dap_H_pooled$chr09),
                                  mcols(CHG_MethylationProfile$CA_42dap_H_pooled$chr10),
                                  mcols(CHG_MethylationProfile$CA_42dap_H_pooled$chr11),
                                  mcols(CHG_MethylationProfile$CA_42dap_H_pooled$chr12) )

CHG_MethylationProfile_AA<- rbind(mcols(CHG_MethylationProfile$AA_28dap_pooled$chr01),
                                    mcols(CHG_MethylationProfile$AA_28dap_pooled$chr02),
                                    mcols(CHG_MethylationProfile$AA_28dap_pooled$chr03),
                                    mcols(CHG_MethylationProfile$AA_28dap_pooled$chr04),
                                    mcols(CHG_MethylationProfile$AA_28dap_pooled$chr05),
                                    mcols(CHG_MethylationProfile$AA_28dap_pooled$chr06),
                                    mcols(CHG_MethylationProfile$AA_28dap_pooled$chr07),
                                    mcols(CHG_MethylationProfile$AA_28dap_pooled$chr08),
                                    mcols(CHG_MethylationProfile$AA_28dap_pooled$chr09),
                                    mcols(CHG_MethylationProfile$AA_28dap_pooled$chr10),
                                    mcols(CHG_MethylationProfile$AA_28dap_pooled$chr11),
                                    mcols(CHG_MethylationProfile$AA_28dap_pooled$chr12) )

CHG_MethylationProfile_AA_H<- rbind(mcols(CHG_MethylationProfile$AA_42dap_H_pooled$chr01),
                                    mcols(CHG_MethylationProfile$AA_42dap_H_pooled$chr02),
                                    mcols(CHG_MethylationProfile$AA_42dap_H_pooled$chr03),
                                    mcols(CHG_MethylationProfile$AA_42dap_H_pooled$chr04),
                                    mcols(CHG_MethylationProfile$AA_42dap_H_pooled$chr05),
                                    mcols(CHG_MethylationProfile$AA_42dap_H_pooled$chr06),
                                    mcols(CHG_MethylationProfile$AA_42dap_H_pooled$chr07),
                                    mcols(CHG_MethylationProfile$AA_42dap_H_pooled$chr08),
                                    mcols(CHG_MethylationProfile$AA_42dap_H_pooled$chr09),
                                    mcols(CHG_MethylationProfile$AA_42dap_H_pooled$chr10),
                                    mcols(CHG_MethylationProfile$AA_42dap_H_pooled$chr11),
                                    mcols(CHG_MethylationProfile$AA_42dap_H_pooled$chr12) )

CHG_MethylationProfile_AA_noNA <- na.omit(CHG_MethylationProfile_AA)
CHG_MethylationProfile_AA_H_noNA <- na.omit(CHG_MethylationProfile_AA_H)
CHG_MethylationProfile_CA_noNA<- na.omit(CHG_MethylationProfile_CA)
CHG_MethylationProfile_CA_H_noNA<- na.omit(CHG_MethylationProfile_CA_H)

summary(CHG_MethylationProfile_AA_noNA$Proportion)
summary(CHG_MethylationProfile_AA_H_noNA$Proportion)
summary(CHG_MethylationProfile_CA_noNA$Proportion)
summary(CHG_MethylationProfile_CA_H_noNA$Proportion)

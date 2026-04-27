library(DMRcaller)
library(GenomicRanges)
library(tidyverse)
library(genomation)

path <- "/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/CX_Files/CHR_only/Context_specific/"



methylationDataList_pooled_CHH <- GRangesList(
  ###PArents
  "AA_28dap_pooled" = readBismarkPool(c(paste0(path,"CHH_CHR_Qualtrim_AA4K_1_bismark_bt2_pe.deduplicated.CX_report.txt"), 
                                    paste0(path,"CHH_CHR_Qualtrim_AA5K_1_bismark_bt2_pe.deduplicated.CX_report.txt"),
                                    paste0(path,"CHH_CHR_Qualtrim_AA15K_1_bismark_bt2_pe.deduplicated.CX_report.txt"))),
  
  "AA_42dap_H_pooled" = readBismarkPool(c(paste0(path,"CHH_CHR_Qualtrim_AA43W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt"), 
                                          paste0(path,"CHH_CHR_Qualtrim_AA53W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt"),
                                          paste0(path,"CHH_CHR_Qualtrim_AA153W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt"))),
  
  "CA_28dap_pooled" = readBismarkPool(c(paste0(path,"CHH_CHR_Qualtrim_CA2K_1_bismark_bt2_pe.deduplicated.CX_report.txt"), 
                                        paste0(path,"CHH_CHR_Qualtrim_CA3_1_bismark_bt2_pe.deduplicated.CX_report.txt"),
                                        paste0(path,"CHH_CHR_Qualtrim_CA11_1_bismark_bt2_pe.deduplicated.CX_report.txt"))),
  
  "CA_42dap_H_pooled" = readBismarkPool(c(paste0(path,"CHH_CHR_Qualtrim_CA23W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt"), 
                                          paste0(path,"CHH_CHR_Qualtrim_CA33W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt"),
                                          paste0(path,"CHH_CHR_Qualtrim_CA113W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt")))
  
)




filter_scaffolds <- function(gr) {
  return(subset(gr, !grepl("scaffold",seqnames(gr))))
}
##CHH
methylationDataList_pooled_CHH_chr <-lapply(methylationDataList_pooled_CHH, 
                                                     filter_scaffolds) 


computeDMRGaussian <- function(GrangesReference,GrangesComparison,Context,minProportionDiff,DMRsize,pvalue){
  
  DMRGaussian  <- computeDMRs(GrangesReference,
                              GrangesComparison,
                              regions = NULL,
                              context = Context,
                              method = "noise_filter",
                              windowSize = 100,
                              kernelFunction = "gaussian",
                              test = "fisher", #change back to fisher
                              pValueThreshold = pvalue,
                              minCytosinesCount = 4, # lower parameter, 
                              minProportionDifference = minProportionDiff, #0.4
                              minGap = 200,
                              minSize = DMRsize, #50 
                              minReadsPerCytosine = 8,
                              cores = 16)
  
  return(DMRGaussian)
}

###RUN CHH DMR
##AAvsAA_H
CHH_AA28vsAA42_H<- computeDMRGaussian(methylationDataList_pooled_CHH_chr[["AA_28dap_pooled"]],
                                      methylationDataList_pooled_CHH_chr[["AA_42dap_H_pooled"]],
                                  Context="CHH",
                                  minProportionDiff = 0.1,
                                  DMRsize = 25, ## accordingly size of CHH methylation is not huge, standard was 50
                                  pvalue=0.05)
write.table(as.data.frame(CHH_AA28vsAA42_H),"/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/CHH_DMRs/FISHER/CHH_AA28vsAA42_H.txt",
            col.names = TRUE)

##CAvsCA_H
CHH_CA28vsCA42_H <- computeDMRGaussian(methylationDataList_pooled_CHH_chr[["CA_28dap_pooled"]],
                                       methylationDataList_pooled_CHH_chr[["CA_42dap_H_pooled"]],
                                  Context="CHH",
                                  minProportionDiff = 0.1,
                                  DMRsize = 25, ## accordingly size of CHH methylation is not huge, standard was 50
                                  pvalue=0.05)
write.table(as.data.frame(CHH_CA28vsCA42_H),"/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/CHH_DMRs/FISHER/CHH_CA28vsCA42_H.txt",
            col.names = TRUE)


####DMRs CHH... ##convert to bed, cvh
CHH_AA28vsAA42_H <- read.csv("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/CHH_DMRs/CHH_AA28vsAA42_H.txt", sep="")
CHH_CA28vsCA42_H <- read.csv("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/CHH_DMRs/CHH_CA28vsCA42_H.txt", sep="")

CHH_AA28vsAA42_H_GR <- GRanges(
  seqnames = CHH_AA28vsAA42_H$seqnames,
  ranges = IRanges(start = CHH_AA28vsAA42_H$start, end = CHH_AA28vsAA42_H$end),
  strand = CHH_AA28vsAA42_H$strand,
  score = CHH_AA28vsAA42_H$direction
)

export(CHH_AA28vsAA42_H_GR, "/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/CHH_DMRs/FISHER/CHH_AA28vsAA42_H.bed", format = "BED")

CHH_CA28vsCA42_H_GR <- GRanges(
  seqnames = CHH_CA28vsCA42_H$seqnames,
  ranges = IRanges(start = CHH_CA28vsCA42_H$start, end = CHH_CA28vsCA42_H$end),
  strand = CHH_CA28vsCA42_H$strand,
  score = CHH_CA28vsCA42_H$direction
)

export(CHH_CA28vsCA42_H_GR, "/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/CHH_DMRs/FISHER/CHH_CA28vsCA42_H.bed", format = "BED")

library(rtracklayer)

##AA vs CA C
CHH_AA28_CvsCA28_C<- computeDMRGaussian(methylationDataList_pooled_CHH_chr[["AA_28dap_pooled"]],
                                        methylationDataList_pooled_CHH_chr[["CA_28dap_pooled"]],
                                  Context="CHH",
                                  minProportionDiff = 0.1,
                                  DMRsize = 25, ## accordingly size of CHH methylation is not huge, standard was 50
                                  pvalue=0.05)

write.table(as.data.frame(CHH_AA28_CvsCA28_C),"/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/CHH_DMRs/FISHER/CHH_AA28_CvsCA28_C.txt",
            col.names = TRUE)

##AA vs CA H
CHH_AA42_HvsCA42_H<- computeDMRGaussian(methylationDataList_pooled_CHH_chr[["AA_42dap_H_pooled"]],
                                        methylationDataList_pooled_CHH_chr[["CA_42dap_H_pooled"]],
                                  Context="CHH",
                                  minProportionDiff = 0.1,
                                  DMRsize = 25, ## accordingly size of CHH methylation is not huge, standard was 50
                                  pvalue=0.05)

write.table(as.data.frame(CHH_AA42_HvsCA42_H),"/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/CHH_DMRs/FISHER/CHH_AA42_HvsCA42_H.txt",
            col.names = TRUE)


S.tuberosum.v6.1_latest <- read.csv("/media/rna/Epipotato16TB/00_META_DATA/S.tuberosum.v6.1_latest.txt", 
                                    sep="")
GENES_df <- read.delim("/media/rna/Epipotato16TB/00_META_DATA/GENES_df.txt")
colnames(GENES_df)[1] <- "seqnames"
colnames(GENES_df)[12] <- "locusName"
GENES_df$locusName <- gsub("Name=","",GENES_df$locusName)

DMR_CHH_list <- list(
               AA_CvH = read.csv("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/CHH_DMRs/FISHER/CHH_AA28vsAA42_H.txt", sep=""),
               CA_CvH = read.csv("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/CHH_DMRs/FISHER/CHH_CA28vsCA42_H.txt", sep=""),
               AAvCA_C = read.csv("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/CHH_DMRs/FISHER/CHH_AA28_CvsCA28_C.txt", sep=""),
               AAvCA_H = read.csv("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/CHH_DMRs/FISHER/CHH_AA42_HvsCA42_H.txt", sep=""))

for (DMR_comparison in names(DMR_CHH_list)) {
  
  print(nrow(as.data.frame(DMR_CHH_list[[DMR_comparison]])))
}

FindDMRoverlaps <- function(GRANGESlist,Gene_annotation) {
  
  ##first create a list with both methods of DMRs results
  DMRslistsAllMethods<- GRANGESlist
  
  ###
  S.tuberosum_annotation <- S.tuberosum.v6.1_latest %>% 
    select(locusName, v6.1.Description,aktualisierte.annotation..bitte.ergänzen.)
  
  
  ###create emptylist
  DMRtypes_Genes <- list()
  for (MethodType in names(DMRslistsAllMethods)) {
    
    #Select the method
    DMRTableMethod <- DMRslistsAllMethods[[MethodType]]
    
    #full join with gene annotation for next step
    FulljoinDMRsGenes_DMRTableMethod<- full_join(as.data.frame(DMRTableMethod),
                                                 Gene_annotation, by="seqnames")
    
    ####Filter out the DMRs for relevant dmrs in the region
    Genes_DMRTableMethod_genebody<- FulljoinDMRsGenes_DMRTableMethod %>% 
      filter(start.x >= start.y & end.x <= end.y) %>% 
      inner_join(S.tuberosum_annotation, by="locusName")
    
    Genes_DMRTableMethod_US<- FulljoinDMRsGenes_DMRTableMethod %>% 
      filter(start.x >= Upstream2000 & end.x <= start.y) %>% 
      inner_join(S.tuberosum_annotation, by="locusName")
    
    Genes_DMRTableMethod_DS<- FulljoinDMRsGenes_DMRTableMethod %>% 
      filter(start.x >= end.y & end.x <= Downstream2000) %>% 
      inner_join(S.tuberosum_annotation, by="locusName")
    
    ##Save the results to empty list
    DMRtypes_Genes[[MethodType]] <- list(
      GeneBody = Genes_DMRTableMethod_genebody,
      Upstream = Genes_DMRTableMethod_US,
      Downstream = Genes_DMRTableMethod_DS
    )
    
  }
  
  return(DMRtypes_Genes)
  
}

DMR_list_CHH_annot<- FindDMRoverlaps(DMR_CHH_list,
                                 GENES_df)

CorrectForUpstreamDownstream<- function(DMR_Genebody,
                                        DMR_putative_Upstream,
                                        DMR_putative_downstream){
  
  DMR_corrected_Upstream <- rbind(DMR_putative_Upstream %>% dplyr::filter(str_detect(V7,"\\+")), ### upstream +
                                  DMR_putative_downstream %>% dplyr::filter(str_detect(V7,"\\-"))) ### downstream - (based on strand reading frame)
  
  DMR_corrected_Downstream <- rbind(DMR_putative_downstream %>% dplyr::filter(str_detect(V7,"\\+")),
                                    DMR_putative_Upstream %>% dplyr::filter(str_detect(V7,"\\-")))
    

  Final_Correct_DMR_annot <- list(
   Genebody = DMR_Genebody,
   UpstreamCorrected = DMR_corrected_Upstream,
   DownstreamCorrected = DMR_corrected_Downstream
  )
  
}

DMR_list_CHH_annot_corrected <- list()
for (DMR_comparison in names(DMR_list_CHH_annot)) {
  
  DMR_list_CHH_annot_corrected[[DMR_comparison]]  <- CorrectForUpstreamDownstream(
                                     DMR_Genebody = DMR_list_CHH_annot[[DMR_comparison]]$GeneBody,
                                     DMR_putative_Upstream = DMR_list_CHH_annot[[DMR_comparison]]$Upstream,
                                     DMR_putative_downstream = DMR_list_CHH_annot[[DMR_comparison]]$Downstream
  )
  
}

SaveDMRAnnotationTable <- function(DMR_ComparisonTableList,
                                   path){
  
  for (DMR_Comparison in names(DMR_ComparisonTableList)) {
    
    write.table(DMR_ComparisonTableList[[DMR_Comparison]],
                paste0(path,DMR_Comparison,".txt"),
                col.names = TRUE)
  }
}

##save each annotation table
SaveDMRAnnotationTable(DMR_list_CHH_annot_corrected$AA_CvH, "/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/CHH_DMRs/FISHER/Annabelle_CvH/")
SaveDMRAnnotationTable(DMR_list_CHH_annot_corrected$CA_CvH, "/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/CHH_DMRs/FISHER/Camel_CvH/")
SaveDMRAnnotationTable(DMR_list_CHH_annot_corrected$AAvCA_C, "/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/CHH_DMRs/FISHER/Control/")
SaveDMRAnnotationTable(DMR_list_CHH_annot_corrected$AAvCA_H, "/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/CHH_DMRs/FISHER/Heat/")



path <- "/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/CX_Files/CHR_only/Context_specific/"
path2 <- "/media/rna/Epipotato16TB/Epipotato_methylome/04_BISMARK/METHYLATION_CX_REPORT/"


###verify DMRs with first exp
methylationDataList_AnnabelleCHH <- GRangesList(
  ###PArents
  "AA_28dap_C_pooled" = readBismarkPool(c(paste0(path,"CHH_CHR_Qualtrim_AA4K_1_bismark_bt2_pe.deduplicated.CX_report.txt"), 
                                        paste0(path,"CHH_CHR_Qualtrim_AA5K_1_bismark_bt2_pe.deduplicated.CX_report.txt"),
                                        paste0(path,"CHH_CHR_Qualtrim_AA15K_1_bismark_bt2_pe.deduplicated.CX_report.txt"))),
  
  "AA_48dap_H_pooled" = readBismarkPool(c(paste0(path,"CHH_CHR_Qualtrim_AA43W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt"), 
                                          paste0(path,"CHH_CHR_Qualtrim_AA53W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt"),
                                          paste0(path,"CHH_CHR_Qualtrim_AA153W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt"))),
  
  "CA_28dap_C_pooled" = readBismarkPool(c(paste0(path,"CHH_CHR_Qualtrim_CA2K_1_bismark_bt2_pe.deduplicated.CX_report.txt"), 
                                        paste0(path,"CHH_CHR_Qualtrim_CA3_1_bismark_bt2_pe.deduplicated.CX_report.txt"),
                                        paste0(path,"CHH_CHR_Qualtrim_CA11_1_bismark_bt2_pe.deduplicated.CX_report.txt"))),
  
  "CA_48dap_H_pooled" = readBismarkPool(c(paste0(path,"CHH_CHR_Qualtrim_CA23W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt"), 
                                          paste0(path,"CHH_CHR_Qualtrim_CA33W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt"),
                                          paste0(path,"CHH_CHR_Qualtrim_CA113W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt"))),
  
  "AA_42dap_C_pooled" = readBismarkPool(c(paste0(path2,"_AA_A_1_bismark_bt2_pe.ded_CHH.txt"), 
                                    paste0(path2,"_AA_B_1_bismark_bt2_pe.ded_CHH.txt"))),
  
  "AA_42dap_H_pooled" = readBismarkPool(c(paste0(path2,"_AA_F_1_bismark_bt2_pe.ded_CHH.txt"), 
                                    paste0(path2,"_AA_G_1_bismark_bt2_pe.ded_CHH.txt"))),
  
  "CA_42dap_C_pooled" = readBismarkPool(c(paste0(path2,"_CA_A_1_bismark_bt2_pe.ded_CHH.txt"), 
                                    paste0(path2,"_CA_B_1_bismark_bt2_pe.ded_CHH.txt"))),
  
  "CA_42dap_H_pooled" = readBismarkPool(c(paste0(path2,"_CA_F_1_bismark_bt2_pe.ded_CHH.txt"), 
                                    paste0(path2,"_CA_G_1_bismark_bt2_pe.ded_CHH.txt")))
  
  
)

methylationDataList_control_CHH_Annabelle_single <- GRangesList(
  ###PArents
  "AA_C_28_1" = readBismark(paste0(path,"CHH_CHR_Qualtrim_AA4K_1_bismark_bt2_pe.deduplicated.CX_report.txt")), 
  "AA_C_28_2" = readBismark(paste0(path,"CHH_CHR_Qualtrim_AA5K_1_bismark_bt2_pe.deduplicated.CX_report.txt")),
  "AA_C_28_3" = readBismark(paste0(path,"CHH_CHR_Qualtrim_AA15K_1_bismark_bt2_pe.deduplicated.CX_report.txt")), 
  
  "AA_H_48_1" = readBismark(paste0(path,"CHH_CHR_Qualtrim_AA43W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt")), 
  "AA_H_48_2" = readBismark(paste0(path,"CHH_CHR_Qualtrim_AA53W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt")),
  "AA_H_48_3" = readBismark(paste0(path,"CHH_CHR_Qualtrim_AA153W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt")), 
  
  "AA_C_42_1" = readBismark(paste0(path2,"_AA_A_1_bismark_bt2_pe.ded_CHH.txt")), 
  "AA_C_42_2" = readBismark(paste0(path2,"_AA_B_1_bismark_bt2_pe.ded_CHH.txt")),
  
  "AA_H_42_1" = readBismark(paste0(path2,"_AA_F_1_bismark_bt2_pe.ded_CHH.txt")), 
  "AA_H_42_2" = readBismark(paste0(path2,"_AA_G_1_bismark_bt2_pe.ded_CHH.txt"))
  
  )

methylationDataList_CHH_Camel_single <- GRangesList(
  ###PArents
  "CA_C_28_1" = readBismark(paste0(path,"CHH_CHR_Qualtrim_CA2K_1_bismark_bt2_pe.deduplicated.CX_report.txt")), 
  "CA_C_28_2" = readBismark(paste0(path,"CHH_CHR_Qualtrim_CA3_1_bismark_bt2_pe.deduplicated.CX_report.txt")),
  "CA_C_28_3" = readBismark(paste0(path,"CHH_CHR_Qualtrim_CA11_1_bismark_bt2_pe.deduplicated.CX_report.txt")), 
  
  "CA_H_48_1" = readBismark(paste0(path,"CHH_CHR_Qualtrim_CA23W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt")), 
  "CA_H_48_2" = readBismark(paste0(path,"CHH_CHR_Qualtrim_CA33W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt")),
  "CA_H_48_3" = readBismark(paste0(path,"CHH_CHR_Qualtrim_CA113W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt")), 
  
  "CA_C_42_1" = readBismark(paste0(path2,"_CA_A_1_bismark_bt2_pe.ded_CHH.txt")), 
  "CA_C_42_2" = readBismark(paste0(path2,"_CA_B_1_bismark_bt2_pe.ded_CHH.txt")),
  
  "CA_H_42_1" = readBismark(paste0(path2,"_CA_F_1_bismark_bt2_pe.ded_CHH.txt")), 
  "CA_H_42_2" = readBismark(paste0(path2,"_CA_G_1_bismark_bt2_pe.ded_CHH.txt"))
  
)



filter_scaffolds <- function(gr) {
  return(subset(gr, !grepl("scaffold",seqnames(gr))))
}
##CHH
methylationDataList_control_CHH_Annabelle_single_chr <-lapply(methylationDataList_control_CHH_Annabelle_single, 
                                            filter_scaffolds) 



methylationDataList_control_CHH_Annabelle_single_chr_join<-joinReplicates(methylationDataList_control_CHH_Annabelle_single_chr[[1]], 
                                                                      methylationDataList_control_CHH_Annabelle_single_chr[[2]]) 
for (h in 2:9) { #24 samples including replicates hence 10-1
  print(h)
  methylationDataList_control_CHH_Annabelle_single_chr_join<-joinReplicates(methylationDataList_control_CHH_Annabelle_single_chr_join, 
                                                                            methylationDataList_control_CHH_Annabelle_single_chr[[h+1]])
}

test <- methylationDataList_control_CHH_Annabelle_single_chr_join

####Define the column names for methylationDataList_Controls_all_single_chr_genotypes
column_names <- rep(names(methylationDataList_control_CHH_Annabelle_single_chr), each=2)
extra_names <- rep(c("meth", "cov"),times=10)
column_extra_names <- paste(column_names,extra_names, sep ="_")

colnames(mcols(methylationDataList_control_CHH_Annabelle_single_chr_join))[2:21] <- column_extra_names




###Function to get positions of GrangesObject based on DMRs
GetGrangesPositionFromDMRs <-function(DMR_table,GrangesObject){ 
  
  DMR_Genebody_list <- list()
  for (i in 1:nrow(DMR_table)) {
    
    print(paste0("Getting Granges for DMR ",i))
    current_DMR_row <- DMR_table[i, ]
    
    
    subset_GrangesObject <- 
      GrangesObject[seqnames(GrangesObject) == current_DMR_row$seqnames & 
                      start(GrangesObject) >= current_DMR_row$start & 
                      end(GrangesObject) <= current_DMR_row$end]
    
    
    DMR_Genebody_list[[i]] <- subset_GrangesObject
    
  }
  return (DMR_Genebody_list)
}


CHH_DMRs_Selected2ndExpCvH_annabelle<- GetGrangesPositionFromDMRs(CHH_AA28vsAA42_H,
                                                                  methylationDataList_control_CHH_Annabelle_single_chr_join)


CalculatePercentageCpGPositions<- function(DMRs_list,SampleNames){
  
  GR_GeneBodyDMRslist_percentage <- list()
  for (DMR in seq(DMRs_list)) {
    ##Select the DMRs
    
    print(paste0("Calculating for DMR ", DMR))
    
    DMR_tables <- as.data.frame(DMRs_list[[DMR]])
    
    #initialize the percentage df
    DMR_tables_percentage_df <- DMR_tables[,c(1:6,ncol(DMR_tables))]
    for (i in SampleNames) {
      
      ##Select the columns....
      CpG_methylation_Single <- DMR_tables %>% 
        select(contains(i))
      
      ##Remove X from the column names if present
      colnames(CpG_methylation_Single) <- gsub("X","",colnames(CpG_methylation_Single))
      
      ##Calculate the methylation percentage
      CpG_methylation_Single_perc <- CpG_methylation_Single %>% transmute(
        !!paste0(i):= CpG_methylation_Single[,paste0(i,"_meth")] / CpG_methylation_Single[,paste0(i,"_cov")]
      )
      
      DMR_tables_percentage_df <- bind_cols(DMR_tables_percentage_df,
                                            CpG_methylation_Single_perc)
    }
    
    
    GR_GeneBodyDMRslist_percentage[[DMR]] <- DMR_tables_percentage_df
    
  }
  return(GR_GeneBodyDMRslist_percentage)
}

CHH_DMRs_Selected2ndExpCvH_annabelle_perc<- CalculatePercentageCpGPositions(
  CHH_DMRs_Selected2ndExpCvH_annabelle,
  names(methylationDataList_control_CHH_Annabelle_single_chr))

CalculateMeanMethylationDMRs<- function(DMRPercentageList){
  
  ##Initialize the empty list to store results mean methylation Results
  DMR_DF <- data.frame()
  for (DMRs in seq(DMRPercentageList)) {
    
    DMR_table_percentage_CpG<- DMRPercentageList[[DMRs]]
    
    #Get the chr for the DMRs, is in the first column, row number doesnt matter here
    DMR_table_percentage_CpG_Chr <- data.frame(chr=as.character(DMR_table_percentage_CpG[1,1]))
    
    ##Calculate the Mean methylation of all samples
    DMR_table_percentage_mean <- t(as.data.frame(colMeans(DMR_table_percentage_CpG[,8:ncol(DMR_table_percentage_CpG)])))
    rownames(DMR_table_percentage_mean) <- DMRs
    
    ##Get the DMR length that was based on the parents
    DMR_length <- data.frame(start = DMR_table_percentage_CpG[1,2],
                             end = DMR_table_percentage_CpG[nrow(DMR_table_percentage_CpG),2])
    
    #Combine the information together to a single row DF for that DMR
    DMR_mean_methylation <- cbind(DMR_table_percentage_CpG_Chr,DMR_length,DMR_table_percentage_mean)
    
    ##Append to the empty dataframe to store mean methylation results for all DMR regions
    DMR_DF <- rbind(DMR_DF,DMR_mean_methylation)
    
  }
  
  return(DMR_DF)
  
}

CHH_DMRs_Selected2ndExpCvH_annabelle_perc_mean<- CalculateMeanMethylationDMRs(
  CHH_DMRs_Selected2ndExpCvH_annabelle_perc
)


##Create metadata for labels
MetaInfo <- data.frame(
  SampleID= colnames(CHH_DMRs_Selected2ndExpCvH_annabelle_perc_mean)[4:ncol(CHH_DMRs_Selected2ndExpCvH_annabelle_perc_mean)]
)

MetaInfo <- MetaInfo %>%
  as.data.frame() %>%
  mutate(Genotypes=case_when(str_detect(SampleID,"AA") ~ "AA",
                             str_detect(SampleID,"CA") ~ "CA"),
         TimePoint= case_when(str_detect(SampleID,"28") ~ "28dap",
                              str_detect(SampleID,"42") ~ "42dap",
                              str_detect(SampleID,"48") ~ "48dap"),
         Conditions=case_when(str_detect(SampleID,"C") ~ "Control",
                              str_detect(SampleID,"H") ~ "Heat"))

library(circlize)
col_fun <- colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
ComplexHeatmap::Heatmap(
  as.matrix(na.omit(CHH_DMRs_Selected2ndExpCvH_annabelle_perc_mean[4:ncol(CHH_DMRs_Selected2ndExpCvH_annabelle_perc_mean)])),
  cluster_columns = F,
  cluster_rows =T,
  heatmap_legend_param = list(
    title = "Methylation"
  ), col = colorRamp2(c(0, 0.3, 1), c("blue", "white", "red")),
  show_row_names = F
)


PCA_DMR_Annabelle<- PCA(t(na.omit(CHH_DMRs_Selected2ndExpCvH_annabelle_perc_mean[4:ncol(CHH_DMRs_Selected2ndExpCvH_annabelle_perc_mean)])),
                          graph = FALSE,
                          scale.unit = TRUE)

fviz_pca_ind(PCA_DMR_Annabelle,
             # Fill individuals by groups
             geom.ind = "point", addEllipses = FALSE,
             col.var = "black",
             repel = TRUE) + geom_point(aes(shape = factor(MetaInfo$Conditions),
                                            colour = factor(MetaInfo$Genotypes)),size = 6) +
  guides(colour=guide_legend("Genotypes"),
         shape = guide_legend("Conditions")) +
  ggtitle("CHH DMRs") +
  geom_text(label = MetaInfo$SampleID, vjust=1.7)+
  scale_colour_manual(values=c("AA" = "#CC0000",
                               "CA" ="lightblue")) + theme(
                                 legend.title = element_text(face = "bold",size = 15),
                                 legend.text = element_text(face = "bold",size = 15),
                                 axis.title = element_text(face = "bold",size = 15),
                                 axis.text = element_text(face = "bold",size=12))


###### Inspect ABS difference of the DMRs identified in the 2nd exp for 1st exp #####


##### Verify candidate genes DMRs ####


GetGrangesRegionMethylation <-  function(DMRcallerObject, chr, START, END, number_samples){
  
  Region_MethylationCGlist <- lapply(DMRcallerObject, function(gr)
  {
    subset_gr <- gr[seqnames(gr) == chr & start(gr) >= START & end(gr) <= END]                
    
    return(subset_gr)
  })
  
  loops <- number_samples-1
  ##join each sample together
  Region_MethylationCGlist_GR<-joinReplicates(Region_MethylationCGlist[[1]], Region_MethylationCGlist[[2]]) 
  for (h in 2:loops) { #10 samples hence 10-1
    print(h)
    Region_MethylationCGlist_GR<-joinReplicates(Region_MethylationCGlist_GR, Region_MethylationCGlist[[h+1]])
  }
  return(Region_MethylationCGlist_GR)
}

##verify for one candidate GENE DNAJ HSP 02G017060
methylationDataList_control_CHH_Annabelle_DNAJ_HSP <- GetGrangesRegionMethylation(methylationDataList_control_CHH_Annabelle_single,
                                                                                  chr="chr02",
                                                                                  START = 31553584,
                                                                                  END = 31554044,
                                                                                  number_samples=10)

methylationDataList_control_CHH_Camel_DNAJ_HSP <- GetGrangesRegionMethylation(methylationDataList_CHH_Camel_single,
                                                                              chr="chr02",
                                                                              START = 31553584,
                                                                              END = 31554044,
                                                                              number_samples=10)

##edit such the the meth and cov of each sample will be labelled
column_names <- rep(names(methylationDataList_control_CHH_Annabelle_single), each=2)
extra_names <- rep(c("meth", "cov"),times=10)
column_extra_names <- paste(column_names,extra_names, sep ="_")

##rename the column names, not sure how to automate this...
names(mcols(methylationDataList_control_CHH_Annabelle_DNAJ_HSP))[2:21]<- column_extra_names


DMR_list <- list(
  AA_CvH = CHH_AA28vsAA42_H
)

plotLocalMethylationProfile(methylationDataList_CHH_Camel_single$CA_C_28_1,
                            methylationDataList_CHH_Camel_single$CA_H_48_1,
                            GRanges(seqnames = Rle("chr02"), ranges = IRanges(31553584,31554044)),
                            DMR_list,
                            context = "CHH",
                            conditionsNames = c("CA_C_28_1", "CA_H_48_1"),
                            windowSize = 200,
                            main="CHH methylation")

plotLocalMethylationProfile(methylationDataList_CHH_Camel_single$CA_C_42_2,
                            methylationDataList_CHH_Camel_single$CA_H_42_2,
                            GRanges(seqnames = Rle("chr02"), ranges = IRanges(31553584,31554044)),
                            DMR_list,
                            context = "CHH",
                            conditionsNames = c("CA_C_42_2", "CA_H_42_2"),
                            windowSize = 200,
                            main="CHH methylation")

###calculate the methylation
###calculate methylation perc
CalculateCpGMethylationSites<- function(GR_DF,SampleNames){
  
  GR_DF_matrix<- GR_DF %>% 
    select(-seqnames,-end,-width,-strand,-context,-trinucleotide_context) %>% 
    column_to_rownames("start") %>% as.matrix()
  
  Calculated_CpG_final <- data.frame(start=rownames(GR_DF_matrix))
  
  for (samples in SampleNames) {
    
    Selected_samples<- GR_DF_matrix %>% as.data.frame() %>%
      select(contains(samples))
    
    Calculated_CpG<- Selected_samples %>% transmute(
      !!paste0(samples):= Selected_samples[,paste0(samples,"_meth")] / Selected_samples[,paste0(samples,"_cov")]
    )
    
    Calculated_CpG_final <- cbind(Calculated_CpG_final,Calculated_CpG)
    
  }
  
  return(Calculated_CpG_final)
}

##Now have to calculate mean methylation of all CpG positions for all DMRs using colMeans
CalculateMeanMethylationDMRs<- function(DMRPercentageList){
  
  ##Initialize the empty list to store results mean methylation Results
  DMR_DF <- data.frame()
  for (DMRs in seq(DMRPercentageList)) {
    
    DMR_table_percentage_CpG<- DMRPercentageList[[DMRs]]
    
    #Get the chr for the DMRs, is in the first column, row number doesnt matter here
    DMR_table_percentage_CpG_Chr <- data.frame(chr=as.character(DMR_table_percentage_CpG[1,1]))
    
    ##Calculate the Mean methylation of all samples
    DMR_table_percentage_mean <- t(as.data.frame(colMeans(DMR_table_percentage_CpG[,8:ncol(DMR_table_percentage_CpG)])))
    rownames(DMR_table_percentage_mean) <- DMRs
    
    ##Get the DMR length that was based on the parents
    DMR_length <- data.frame(start = DMR_table_percentage_CpG[1,2],
                             end = DMR_table_percentage_CpG[nrow(DMR_table_percentage_CpG),2])
    
    #Combine the information together to a single row DF for that DMR
    DMR_mean_methylation <- cbind(DMR_table_percentage_CpG_Chr,DMR_length,DMR_table_percentage_mean)
    
    ##Append to the empty dataframe to store mean methylation results for all DMR regions
    DMR_DF <- rbind(DMR_DF,DMR_mean_methylation)
    
  }
  
  return(DMR_DF)
  
}



methylationDataList_control_CHH_Annabelle_DNAJ_HSP_perc<- CalculateCpGMethylationSites(as.data.frame(methylationDataList_control_CHH_Annabelle_DNAJ_HSP),
                             names(methylationDataList_control_CHH_Annabelle_single))

methylationDataList_control_CHH_Camel_DNAJ_HSP_perc<- CalculateCpGMethylationSites(as.data.frame(methylationDataList_control_CHH_Camel_DNAJ_HSP),
                                                                                       names(methylationDataList_CHH_Camel_single))

test <- CHHstepsDNAJ %>% as.data.frame() %>% 
  rownames_to_column("start") %>% filter(
    start > 31553684 & start < 31554044
  ) %>% column_to_rownames("start")

###try to calculate rolling mean
library(zoo)
window_size <- 10

CHHstepsDNAJ <- methylationDataList_control_CHH_Annabelle_DNAJ_HSP_perc %>%
  arrange(desc(start)) %>%  # Ensure sorted order
  transmute(
    AA_C_28_1 = rollmean(AA_C_28_1, window_size, fill = NA, align = "center"),
    AA_C_28_2 = rollmean(AA_C_28_2, window_size, fill = NA, align = "center"),
    AA_C_28_3 = rollmean(AA_C_28_3, window_size, fill = NA, align = "center"),
    AA_H_48_1 = rollmean(AA_H_48_1, window_size, fill = NA, align = "center"),
    AA_H_48_2 = rollmean(AA_H_48_2, window_size, fill = NA, align = "center"),
    AA_H_48_3 = rollmean(AA_H_48_3, window_size, fill = NA, align = "center"),
    
    AA_C_42_1 = rollmean(AA_C_42_1, window_size, fill = NA, align = "center"),
    AA_C_42_2 = rollmean(AA_C_42_2, window_size, fill = NA, align = "center"),
    
    AA_H_42_1 = rollmean(AA_H_42_1, window_size, fill = NA, align = "center"),
    AA_H_42_2 = rollmean(AA_H_42_2, window_size, fill = NA, align = "center")
  )

CHHstepsDNAJ_CA <- methylationDataList_control_CHH_Camel_DNAJ_HSP_perc %>%
  arrange(desc(start)) %>%  # Ensure sorted order
  transmute(
    CA_C_28_1 = rollmean(CA_C_28_1, window_size, fill = NA, align = "center"),
    CA_C_28_2 = rollmean(CA_C_28_2, window_size, fill = NA, align = "center"),
    CA_C_28_3 = rollmean(CA_C_28_3, window_size, fill = NA, align = "center"),
    CA_H_48_1 = rollmean(CA_H_48_1, window_size, fill = NA, align = "center"),
    CA_H_48_2 = rollmean(CA_H_48_2, window_size, fill = NA, align = "center"),
    CA_H_48_3 = rollmean(CA_H_48_3, window_size, fill = NA, align = "center"),
    
    CA_C_42_1 = rollmean(CA_C_42_1, window_size, fill = NA, align = "center"),
    CA_C_42_2 = rollmean(CA_C_42_2, window_size, fill = NA, align = "center"),
    
    CA_H_42_1 = rollmean(CA_H_42_1, window_size, fill = NA, align = "center"),
    CA_H_42_2 = rollmean(CA_H_42_2, window_size, fill = NA, align = "center")
  )

library(ComplexHeatmap)
library(circlize)
Heatmap(
  CHHstepsDNAJ_CA %>% as.matrix() %>% na.omit() %>%t(),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 16,fontface="bold"),
  heatmap_legend_param = list(title = "Methylation",
                              title_gp = gpar(fontsize = 14, fontface = "bold"),  # Increase title size & bold
                              labels_gp = gpar(fontsize = 12, fontface = "bold")),
  width = unit(45, "cm"),
  height = unit(5, "cm"),
  col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
)


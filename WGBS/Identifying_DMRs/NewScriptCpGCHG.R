library(DMRcaller)
library(GenomicRanges)
library(tidyverse)
library(genomation)
library(zoo)
path <- "/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/WGBS_HIGH_HEAT/03_BISMARK/CX_FILES/Chr_Context/"


CG_CHR_Qualtrim_AA4K_1_bismark_bt2_pe.deduplicated.CX_report.txt
CG_CHR_Qualtrim_AA5K_1_bismark_bt2_pe.deduplicated.CX_report.txt
CG_CHR_Qualtrim_AA15K_1_bismark_bt2_pe.deduplicated.CX_report.txt

CG_CHR_Qualtrim_AA43W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt
CG_CHR_Qualtrim_AA53W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt
CG_CHR_Qualtrim_AA153W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt

CG_CHR_Qualtrim_CA2K_1_bismark_bt2_pe.deduplicated.CX_report.txt
CG_CHR_Qualtrim_CA3_1_bismark_bt2_pe.deduplicated.CX_report.txt
CG_CHR_Qualtrim_CA11_1_bismark_bt2_pe.deduplicated.CX_report.txt

CG_CHR_Qualtrim_CA23W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt
CG_CHR_Qualtrim_CA33W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt
CG_CHR_Qualtrim_CA113W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt

methylationDataList_pooled <- GRangesList(
  ###PArents
  "AA_28dap_pooled" = readBismarkPool(c(paste0(path,"CG_CHR_Qualtrim_AA4K_1_bismark_bt2_pe.deduplicated.CX_report.txt"), 
                                    paste0(path,"CG_CHR_Qualtrim_AA5K_1_bismark_bt2_pe.deduplicated.CX_report.txt"),
                                    paste0(path,"CG_CHR_Qualtrim_AA15K_1_bismark_bt2_pe.deduplicated.CX_report.txt"))),
  
  "AA_42dap_H_pooled" = readBismarkPool(c(paste0(path,"CG_CHR_Qualtrim_AA43W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt"), 
                                          paste0(path,"CG_CHR_Qualtrim_AA53W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt"),
                                          paste0(path,"CG_CHR_Qualtrim_AA153W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt"))),
  
  "CA_28dap_pooled" = readBismarkPool(c(paste0(path,"CG_CHR_Qualtrim_CA2K_1_bismark_bt2_pe.deduplicated.CX_report.txt"), 
                                        paste0(path,"CG_CHR_Qualtrim_CA3_1_bismark_bt2_pe.deduplicated.CX_report.txt"),
                                        paste0(path,"CG_CHR_Qualtrim_CA11_1_bismark_bt2_pe.deduplicated.CX_report.txt"))),
  
  "CA_42dap_H_pooled" = readBismarkPool(c(paste0(path,"CG_CHR_Qualtrim_CA23W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt"), 
                                          paste0(path,"CG_CHR_Qualtrim_CA33W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt"),
                                          paste0(path,"CG_CHR_Qualtrim_CA113W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt")))
  
)

methylationDataList_CHG_pooled <- GRangesList(
  ###PArents
  "AA_28dap_pooled" = readBismarkPool(c(paste0(path,"CHG_CHR_Qualtrim_AA4K_1_bismark_bt2_pe.deduplicated.CX_report.txt"), 
                                        paste0(path,"CHG_CHR_Qualtrim_AA5K_1_bismark_bt2_pe.deduplicated.CX_report.txt"),
                                        paste0(path,"CHG_CHR_Qualtrim_AA15K_1_bismark_bt2_pe.deduplicated.CX_report.txt"))),
  
  "AA_42dap_H_pooled" = readBismarkPool(c(paste0(path,"CHG_CHR_Qualtrim_AA43W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt"), 
                                          paste0(path,"CHG_CHR_Qualtrim_AA53W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt"),
                                          paste0(path,"CHG_CHR_Qualtrim_AA153W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt"))),
  
  "CA_28dap_pooled" = readBismarkPool(c(paste0(path,"CHG_CHR_Qualtrim_CA2K_1_bismark_bt2_pe.deduplicated.CX_report.txt"), 
                                        paste0(path,"CHG_CHR_Qualtrim_CA3_1_bismark_bt2_pe.deduplicated.CX_report.txt"),
                                        paste0(path,"CHG_CHR_Qualtrim_CA11_1_bismark_bt2_pe.deduplicated.CX_report.txt"))),
  
  "CA_42dap_H_pooled" = readBismarkPool(c(paste0(path,"CHG_CHR_Qualtrim_CA23W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt"), 
                                          paste0(path,"CHG_CHR_Qualtrim_CA33W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt"),
                                          paste0(path,"CHG_CHR_Qualtrim_CA113W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt")))
  
)

filter_scaffolds <- function(gr) {
  return(subset(gr, !grepl("scaffold",seqnames(gr))))
}
##CG
methylationDataList_pooled_chr <-lapply(methylationDataList_pooled, 
                                                     filter_scaffolds) 
##CHG
methylationDataList_pooled_CHG_chr <-lapply(methylationDataList_CHG_pooled, 
                                        filter_scaffolds) 

computeDMRGaussian <- function(GrangesReference,GrangesComparison,Context,minProportionDiff,DMRsize,pvalue){
  
  DMRGaussian  <- computeDMRs(GrangesReference,
                              GrangesComparison,
                              regions = NULL,
                              context = Context,
                              method = "noise_filter",
                              windowSize = 100,
                              kernelFunction = "gaussian",
                              test = "fisher", ##change back to fisher
                              pValueThreshold = pvalue,
                              minCytosinesCount = 4, # lower parameter, 
                              minProportionDifference = minProportionDiff, #0.4
                              minGap = 200,
                              minSize = DMRsize, #50 
                              minReadsPerCytosine = 8,
                              cores = 12)
  
  return(DMRGaussian)
}

###RUN CG DMR
##AAvsAA_H
AA28vsAA42_H<- computeDMRGaussian(methylationDataList_pooled_chr[["AA_28dap_pooled"]],
                                  methylationDataList_pooled_chr[["AA_42dap_H_pooled"]],
                                  Context="CG",
                                  minProportionDiff = 0.4,
                                  DMRsize = 50,
                                  pvalue=0.05)
write.table(as.data.frame(AA28vsAA42_H),"/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/WGBS_HIGH_HEAT/ANALYSIS/CpG_DMRs/FISHER/AA28vsAA42_H.txt",
            col.names = TRUE)

##CAvsCA_H score....
CA28vsCA42_H <- computeDMRGaussian(methylationDataList_pooled_chr[["CA_28dap_pooled"]],
                                   methylationDataList_pooled_chr[["CA_42dap_H_pooled"]],
                                  Context="CG",
                                  minProportionDiff = 0.4,
                                  DMRsize = 50,
                                  pvalue=0.05)
write.table(as.data.frame(CA28vsCA42_H),"/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/WGBS_HIGH_HEAT/ANALYSIS/CpG_DMRs/FISHER/CA28vsCA42_H.txt",
            col.names = TRUE)

##AA vs CA C
AA28_CvsCA28_C<- computeDMRGaussian(methylationDataList_pooled_chr[["AA_28dap_pooled"]],
                                    methylationDataList_pooled_chr[["CA_28dap_pooled"]],
                                  Context="CG",
                                  minProportionDiff = 0.4,
                                  DMRsize = 50,
                                  pvalue=0.05)

write.table(as.data.frame(AA28_CvsCA28_C),"/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/WGBS_HIGH_HEAT/ANALYSIS/CpG_DMRs/FISHER/AA28_CvsCA28_C.txt",
            col.names = TRUE)

##AA vs CA H
AA42_HvsCA42_H<- computeDMRGaussian(methylationDataList_pooled_chr[["AA_42dap_H_pooled"]],
                                    methylationDataList_pooled_chr[["CA_42dap_H_pooled"]],
                                  Context="CG",
                                  minProportionDiff = 0.4,
                                  DMRsize = 50,
                                  pvalue=0.05)

write.table(as.data.frame(AA42_HvsCA42_H),"/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/WGBS_HIGH_HEAT/ANALYSIS/CpG_DMRs/FISHER/AA42_HvsCA42_H.txt",
            col.names = TRUE)




###RUN CHG DMR
##AAvsAA_H
CHG_AA28vsAA42_H<- computeDMRGaussian(methylationDataList_pooled_CHG_chr[["AA_28dap_pooled"]],
                                      methylationDataList_pooled_CHG_chr[["AA_42dap_H_pooled"]],
                                  Context="CHG",
                                  minProportionDiff = 0.2,
                                  DMRsize = 50,
                                  pvalue=0.05)
write.table(as.data.frame(CHG_AA28vsAA42_H),"/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/WGBS_HIGH_HEAT/ANALYSIS/CHG_DMRs/FISHER/CHG_AA28vsAA42_H.txt",
            col.names = TRUE)

##CAvsCA_H
CHG_CA28vsCA42_H <- computeDMRGaussian(methylationDataList_pooled_CHG_chr[["CA_28dap_pooled"]],
                                       methylationDataList_pooled_CHG_chr[["CA_42dap_H_pooled"]],
                                   Context="CHG",
                                   minProportionDiff = 0.2,
                                   DMRsize = 50,
                                   pvalue=0.05)
write.table(as.data.frame(CHG_CA28vsCA42_H),"/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/WGBS_HIGH_HEAT/ANALYSIS/CHG_DMRs/FISHER/CHG_CA28vsCA42_H.txt",
            col.names = TRUE)

###save the cvh part as bed file...
AA28vsAA42_H <- read.csv("/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/WGBS_HIGH_HEAT/ANALYSIS/CpG_DMRs/FISHER/AA28vsAA42_H.txt", sep="")
CA28vsCA42_H <- read.csv("/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/WGBS_HIGH_HEAT/ANALYSIS/CpG_DMRs/FISHER/CA28vsCA42_H.txt", sep="")

library(rtracklayer)
AA28vsAA42_H_GR <- GRanges(
  seqnames = AA28vsAA42_H$seqnames,
  ranges = IRanges(start = AA28vsAA42_H$start, end = AA28vsAA42_H$end),
  strand = AA28vsAA42_H$strand,
  score = AA28vsAA42_H$direction
)

export(AA28vsAA42_H_GR, "/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/WGBS_HIGH_HEAT/ANALYSIS/CpG_DMRs/FISHER/AA28vsAA42_H.bed", format = "BED")

CA28vsCA42_H_GR <- GRanges(
  seqnames = CA28vsCA42_H$seqnames,
  ranges = IRanges(start = CA28vsCA42_H$start, end = CA28vsCA42_H$end),
  strand = CA28vsCA42_H$strand,
  score = CA28vsCA42_H$direction
)

export(CA28vsCA42_H_GR, "/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/WGBS_HIGH_HEAT/ANALYSIS/CpG_DMRs/FISHER/CA28vsCA42_H.bed", format = "BED")

##AA vs CA C
CHG_AA28_CvsCA28_C<- computeDMRGaussian(methylationDataList_pooled_CHG_chr[["AA_28dap_pooled"]],
                                        methylationDataList_pooled_CHG_chr[["CA_28dap_pooled"]],
                                    Context="CHG",
                                    minProportionDiff = 0.2,
                                    DMRsize = 50,
                                    pvalue=0.05)

write.table(as.data.frame(CHG_AA28_CvsCA28_C),"/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/WGBS_HIGH_HEAT/ANALYSIS/CHG_DMRs/FISHER/CHG_AA28_CvsCA28_C.txt",
            col.names = TRUE)

##AA vs CA H
CHG_AA42_HvsCA42_H<- computeDMRGaussian(methylationDataList_pooled_CHG_chr[["AA_42dap_H_pooled"]],
                                        methylationDataList_pooled_CHG_chr[["CA_42dap_H_pooled"]],
                                    Context="CHG",
                                    minProportionDiff = 0.2,
                                    DMRsize = 50,
                                    pvalue=0.05)

write.table(as.data.frame(CHG_AA42_HvsCA42_H),"/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/WGBS_HIGH_HEAT/ANALYSIS/CHG_DMRs/FISHER/CHG_AA42_HvsCA42_H.txt",
            col.names = TRUE)



#### reverse make AA the comparison just as confirmation
##AA vs CA C
AA28_CvsCA28_C_test<- computeDMRGaussian(methylationDataList_pooled_chr[["CA_28dap_pooled"]],
                                    methylationDataList_pooled_chr[["AA_28dap_pooled"]],
                                    Context="CG",
                                    minProportionDiff = 0.4,
                                    DMRsize = 50,
                                    pvalue=0.05)

write.table(as.data.frame(AA28_CvsCA28_C_test),"/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/WGBS_HIGH_HEAT/ANALYSIS/CpG_DMRs/FISHER/AA28_CvsCA28_C_test.txt",
            col.names = TRUE)

library(rtracklayer)
names(mcols(AA28_CvsCA28_C_test))[1] <- "score"
export(AA28_CvsCA28_C_test,
       "/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/WGBS_HIGH_HEAT/ANALYSIS/CpG_DMRs/FISHER/AA28_CvsCA28_C_test.bed",
       format = "BED")

###CG methylation annotation

##Function to find overlaps in genic regions
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

colnames(GENES_df)[1] <- "seqnames"
colnames(GENES_df)[12] <- "locusName"
GENES_df$locusName <- gsub("Name=","",GENES_df$locusName)

CG_DMR_list <- list(
  AA_CvH= read.csv("/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/WGBS_HIGH_HEAT/ANALYSIS/CpG_DMRs/FISHER/AA28vsAA42_H.txt", sep=""),
  CA_CvH= read.csv("/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/WGBS_HIGH_HEAT/ANALYSIS/CpG_DMRs/FISHER/CA28vsCA42_H.txt", sep=""),
  AAvCA_C= read.csv("/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/WGBS_HIGH_HEAT/ANALYSIS/CpG_DMRs/FISHER/AA28_CvsCA28_C.txt", sep=""),
  AAvCA_H= read.csv("/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/WGBS_HIGH_HEAT/ANALYSIS/CpG_DMRs/FISHER/AA42_HvsCA42_H.txt", sep="") 
)


for (DMR_comparison in names(CG_DMR_list)) {
  
  print(nrow(as.data.frame(CG_DMR_list[[DMR_comparison]])))
  
}
library(rtracklayer)
library(GenomicRanges)

CG_DMR_list_annot <- FindDMRoverlaps(CG_DMR_list,GENES_df)
### correct for strand information for upstream and downstream region
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

CG_DMR_list_annot_correct <- list()
for (DMR_comparison in names(CG_DMR_list_annot)) {
  
  CG_DMR_list_annot_correct[[DMR_comparison]] <- CorrectForUpstreamDownstream(
    DMR_Genebody = CG_DMR_list_annot[[DMR_comparison]]$GeneBody,
    DMR_putative_Upstream = CG_DMR_list_annot[[DMR_comparison]]$Upstream,
    DMR_putative_downstream = CG_DMR_list_annot[[DMR_comparison]]$Downstream
  )
  
}

###Controls COmparison
for (DMR_comparisons in names(CG_DMR_list_annot_correct$AAvCA_C)) {
  
  write.table(CG_DMR_list_annot_correct$AAvCA_C[[DMR_comparisons]],
              paste0("/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/WGBS_HIGH_HEAT/ANALYSIS/CpG_DMRs/FISHER/Controls/",
                     DMR_comparisons,".txt"),
              col.names = TRUE)
  
}

###Heat COmparison
for (DMR_comparisons in names(CG_DMR_list_annot_correct$AAvCA_H)) {
  
  write.table(CG_DMR_list_annot_correct$AAvCA_H[[DMR_comparisons]],
              paste0("/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/WGBS_HIGH_HEAT/ANALYSIS/CpG_DMRs/FISHER/Heat/",
                     DMR_comparisons,".txt"),
              col.names = TRUE)
  
}

SaveDMRAnnotationTable <- function(DMR_ComparisonTableList,
                                   path){
  
 for (DMR_Comparison in names(DMR_ComparisonTableList)) {
   
   write.table(DMR_ComparisonTableList[[DMR_Comparison]],
               paste0(path,DMR_Comparison,".txt"),
               col.names = TRUE)
   
 }
  
  
}


SaveDMRAnnotationTable(CG_DMR_list_annot_correct$AA_CvH,
  path = "/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/WGBS_HIGH_HEAT/ANALYSIS/CpG_DMRs/FISHER/Annabelle_CvH/"
)

SaveDMRAnnotationTable(CG_DMR_list_annot_correct$CA_CvH,
                       path = "/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/WGBS_HIGH_HEAT/ANALYSIS/CpG_DMRs/FISHER/Camel_CvH/"
)

##load DMR DEGs
AAvsCA_C_1st <- list(
  UpstreamCorrected= read.csv("/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/testDMRCaller/DMR_caller_CG_GenicRegions/NEWEST_1st_EXPERIMENT/CpG/DMR_CALL_PVALUE_005/Control/UpstreamCorrected.txt", sep=""),
  GeneBody = read.csv("/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/testDMRCaller/DMR_caller_CG_GenicRegions/NEWEST_1st_EXPERIMENT/CpG/DMR_CALL_PVALUE_005/Control/Genebody.txt", sep=""),
  DownstreamCorrected = read.csv("/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/testDMRCaller/DMR_caller_CG_GenicRegions/NEWEST_1st_EXPERIMENT/CpG/DMR_CALL_PVALUE_005/Control/DownstreamCorrected.txt", sep="")
)


CG_DMR_list$AAvCA_C
AAvCA_C_GR <- GRanges(
  seqnames = CG_DMR_list$AAvCA_C$seqnames,
  ranges = IRanges(start = CG_DMR_list$AAvCA_C$start, end = CG_DMR_list$AAvCA_C$end),
  strand = CG_DMR_list$AAvCA_C$strand
)

export(AAvCA_C_GR, "/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/WGBS_HIGH_HEAT/ANALYSIS/CpG_DMRs/FISHER/AA28_CvsCA28_C.bed", 
       format = "BED")



##Load DMR genes of 1st experiment (CONTROLS COMPARISON)
AAvsCA_C_1st <- list(
  UpstreamCorrected= read.csv("/media/rna/Epipotato16TB/Analyses/DMR_caller/DMR_caller_CG_test/DMR_CALL_PVALUE_005/Control/UpstreamCorrected.txt", sep=""),
  GeneBody = read.csv("/media/rna/Epipotato16TB/Analyses/DMR_caller/DMR_caller_CG_test/DMR_CALL_PVALUE_005/Control/Genebody.txt", sep=""),
  DownstreamCorrected = read.csv("//media/rna/Epipotato16TB/Analyses/DMR_caller/DMR_caller_CG_test/DMR_CALL_PVALUE_005/Control/DownstreamCorrected.txt", sep="")
)

AAvsCA_H_1st <- list(
  UpstreamCorrected= read.csv("/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/testDMRCaller/DMR_caller_CG_GenicRegions/NEWEST_1st_EXPERIMENT/CpG/DMR_CALL_PVALUE_005/Heat/UpstreamCorrected.txt", sep=""),
  GeneBody = read.csv("/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/testDMRCaller/DMR_caller_CG_GenicRegions/NEWEST_1st_EXPERIMENT/CpG/DMR_CALL_PVALUE_005/Heat/Genebody.txt", sep=""),
  DownstreamCorrected = read.csv("/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/testDMRCaller/DMR_caller_CG_GenicRegions/NEWEST_1st_EXPERIMENT/CpG/DMR_CALL_PVALUE_005/Heat/DownstreamCorrected.txt", sep="")
)

##count total genic DMRs 2nd heat stress exp
for (DMR_comparison in names(CG_DMR_list_annot_correct$AAvCA_C)) {
  
  print(nrow(as.data.frame(CG_DMR_list_annot_correct$AAvCA_C[[DMR_comparison]])))
  
}
##count total genic DMRs 1st heat stress exp
for (DMR_comparison in names(AAvsCA_C_1st)) {
  
  print(nrow(as.data.frame(AAvsCA_C_1st[[DMR_comparison]])))
  
}

###compare overlaps for controls
library(ggvenn)


test <- anti_join(CG_DMR_list_annot_correct$AAvCA_C$UpstreamCorrected,
                  AAvsCA_C_1st$UpstreamCorrected,
          by="locusName"
          )

### load the constitutive DMRs for second experiment 
CG_DMR_list_annot_correct <- list(
  UpstreamCorrected = read.csv("/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/WGBS_HIGH_HEAT/ANALYSIS/CpG_DMRs/FISHER/Controls/UpstreamCorrected.txt", sep=""),
  Genebody = read.csv("/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/WGBS_HIGH_HEAT/ANALYSIS/CpG_DMRs/FISHER/Controls/Genebody.txt", sep=""),
  DownstreamCorrected = read.csv("/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/WGBS_HIGH_HEAT/ANALYSIS/CpG_DMRs/FISHER/Controls/DownstreamCorrected.txt", sep="")
)

ggvenn(
  list(
    Second_28dap= CG_DMR_list_annot_correct$AAvCA_C$UpstreamCorrected$locusName %>% unique(),
    first_42dap = AAvsCA_C_1st$UpstreamCorrected$locusName %>%unique()
  )
)






test1 <- inner_join(unique(test1), S.tuberosum.v6.1_latest)

CG_DMR_list_annot_correct$AAvCA_C$UpstreamCorrected$locusName

###plot sp6a region,
###FIRST AND SECOND EXPERIMENT
path2 <- "/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/WGBS_HIGH_HEAT/03_BISMARK/CX_FILES/Chr_Context/"
path <- "/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/testDMRCaller/CG_CX/"
methylationDataList_control_single <- GRangesList(
  ###PArents
  "AA_28_1" = readBismark(paste0(path2,"CG_CHR_Qualtrim_AA4K_1_bismark_bt2_pe.deduplicated.CX_report.txt")), 
  "AA_28_2" = readBismark(paste0(path2,"CG_CHR_Qualtrim_AA5K_1_bismark_bt2_pe.deduplicated.CX_report.txt")),
  "AA_28_3" = readBismark(paste0(path2,"CG_CHR_Qualtrim_AA15K_1_bismark_bt2_pe.deduplicated.CX_report.txt")), 
  
  "CA_28_1" = readBismark(paste0(path2,"CG_CHR_Qualtrim_CA2K_1_bismark_bt2_pe.deduplicated.CX_report.txt")),
  "CA_28_2" = readBismark(paste0(path2,"CG_CHR_Qualtrim_CA3_1_bismark_bt2_pe.deduplicated.CX_report.txt")),
  "CA_28_3" = readBismark(paste0(path2,"CG_CHR_Qualtrim_CA11_1_bismark_bt2_pe.deduplicated.CX_report.txt")),
  
  "AA_42_1" = readBismark(paste0(path,"_AA_A_1_bismark_bt2_pe.ded_CG.txt")), 
  "AA_42_2" = readBismark(paste0(path,"_AA_B_1_bismark_bt2_pe.ded_CG.txt")),
  
  "CA_42_1" = readBismark(paste0(path,"_CA_A_1_bismark_bt2_pe.ded_CG.txt")), 
  "CA_42_2" = readBismark(paste0(path,"_CA_B_1_bismark_bt2_pe.ded_CG.txt")))

methylationDataList_control_single_chr <-lapply(methylationDataList_control_single, 
                                        filter_scaffolds) 

###GET SP6A region
###Create function to get Granges region and calculate methylation perc of that region

GetGrangesRegionMethylation <-  function(DMRcallerObject, chr,START, END, number_samples){
  
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

SP6A_Region_df_GR<- GetGrangesRegionMethylation(methylationDataList_control_single,
                                                chr="chr05",
                                                START = 54497343,
                                                END = 54503088,
                                                number_samples=10)

DRE_Region_df_GR<- GetGrangesRegionMethylation(methylationDataList_control_single,
                                                chr="chr06",
                                                START = 31370772,
                                                END = 31376864,
                                                number_samples=10)


AB_Hydrolases_Region_df_GR<- GetGrangesRegionMethylation(methylationDataList_control_single,
                                               chr="chr06",
                                               START = 41377854,
                                               END = 41389574,
                                               number_samples=10)

St01G045040_Hydrolases_Region_df_GR<- GetGrangesRegionMethylation(methylationDataList_control_single,
                                                         chr="chr01",
                                                         START = 82882814,
                                                         END = 82888058,
                                                         number_samples=10)


St12G001580_ARR_Region_df_GR<- GetGrangesRegionMethylation(methylationDataList_control_single,
                                                                  chr="chr12",
                                                                  START = 1338072,
                                                                  END = 1346160,
                                                                  number_samples=10)
RenameTheDFtoSampleName<- function(COmbinedGRDMRCaller,namestorename,number_of_samples){
  
  ##Convert to DF
  Region_df_GR<- as.data.frame(COmbinedGRDMRCaller)
  
  column_names <- namestorename
  
  ##edit such the the meth and cov of each sample will be labelled
  column_names <- rep(column_names, each=2)
  extra_names <- rep(c("meth", "cov"),times=number_of_samples)
  column_extra_names <- paste(column_names,extra_names, sep ="_")
  
  dummy_end_of_column_to_rename <- (length(column_extra_names) + 7)-1 ##7 is the start of the column to rename
  
  ##rename the column names, not sure how to automate this...
  colnames(Region_df_GR)[7:dummy_end_of_column_to_rename] <- column_extra_names
  
  return(Region_df_GR)
}
##rename the columns
DRE_Region_df_GR<- RenameTheDFtoSampleName(DRE_Region_df_GR,
                                           namestorename = names(methylationDataList_control_single),
                                           number_of_samples = 10)

AB_Hydrolases_Region_df_GR<- RenameTheDFtoSampleName(AB_Hydrolases_Region_df_GR,
                                           namestorename = names(methylationDataList_control_single),
                                           number_of_samples = 10)

St01G045040_Hydrolases_Region_df_GR<- RenameTheDFtoSampleName(St01G045040_Hydrolases_Region_df_GR,
                                                     namestorename = names(methylationDataList_control_single),
                                                     number_of_samples = 10)

St12G001580_ARR_Region_df_GR<- RenameTheDFtoSampleName(St12G001580_ARR_Region_df_GR,
                                                              namestorename = names(methylationDataList_control_single),
                                                              number_of_samples = 10)

##get SP6A GRANGES
SP6A_Region_df_GR <- as.data.frame(SP6A_Region_df_GR)
column_names <- names(methylationDataList_control_single)

column_names <- rep(column_names, each=2)
extra_names <- rep(c("meth", "cov"),times=10)
column_extra_names <- paste(column_names,extra_names, sep ="_")

colnames(SP6A_Region_df_GR)[7:26] <- column_extra_names

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

DRE_Region_CpG_methylation<- CalculateCpGMethylationSites(DRE_Region_df_GR,
                                names(methylationDataList_control_single))

AB_Hydrolases_Region_CpG_methylation<- CalculateCpGMethylationSites(AB_Hydrolases_Region_df_GR,
                                                          names(methylationDataList_control_single))

St01G045040_Region_CpG_methylation<- CalculateCpGMethylationSites(St01G045040_Hydrolases_Region_df_GR,
                                                                    names(methylationDataList_control_single))

St12G001580_Region_CpG_methylation<- CalculateCpGMethylationSites(St12G001580_ARR_Region_df_GR,
                                                                  names(methylationDataList_control_single))

##cALCULATE THE METHYLATION PERCENTAGE
SP6A_Region_df_GR_matrix<- SP6A_Region_df_GR %>% 
  select(-seqnames,-end,-width,-strand,-context,-trinucleotide_context) %>% 
  column_to_rownames("start") %>% as.matrix()


Calculated_CpG_SP6A_final <- data.frame(start=rownames(SP6A_Region_df_GR_matrix))
for (samples in names(SP6A_Region_MethylationCGlist)) {
  
  Selected_samples<- SP6A_Region_df_GR_matrix %>% as.data.frame() %>%
    select(contains(samples))
  
  Calculated_CpG_SP6A<- Selected_samples %>% transmute(
    !!paste0(samples):= Selected_samples[,paste0(samples,"_meth")] / Selected_samples[,paste0(samples,"_cov")]
  )
  
  Calculated_CpG_SP6A_final <- cbind(Calculated_CpG_SP6A_final,Calculated_CpG_SP6A)
  
}

##plot single regions sp6a
Heatmap(Calculated_CpG_SP6A_final %>% dplyr::select(contains("AA"),contains("CA"),start) %>% filter(start >= 54501100 & start <= 54501380) %>% arrange(desc(start)) %>% select(-start) %>%t(),
        show_row_names=TRUE,cluster_columns  = FALSE, cluster_rows =FALSE,
        heatmap_legend_param = list(title = "Methylation"),
        row_names_gp = gpar(fontsize = 16,fontface="bold"),
        column_names_gp = gpar(fontsize = 16,fontface="bold"),
        width = unit(5, "cm"),
        height = unit(5, "cm"))

###try to calculate rolling mean
Calculated_CpG_SP6A_final
library(zoo)

window_size <- 10

CpGstepsSP6A <- Calculated_CpG_SP6A_final %>%
  arrange(desc(start)) %>%  # Ensure sorted order
  transmute(
    AA_28_1 = rollmean(AA_28_1, window_size, fill = NA, align = "center"),
    AA_28_2 = rollmean(AA_28_2, window_size, fill = NA, align = "center"),
    AA_28_3 = rollmean(AA_28_3, window_size, fill = NA, align = "center"),
    CA_28_1 = rollmean(CA_28_1, window_size, fill = NA, align = "center"),
    CA_28_2 = rollmean(CA_28_2, window_size, fill = NA, align = "center"),
    CA_28_3 = rollmean(CA_28_3, window_size, fill = NA, align = "center"),
    
    AA_42_1 = rollmean(AA_42_1, window_size, fill = NA, align = "center"),
    AA_42_2 = rollmean(AA_42_2, window_size, fill = NA, align = "center"),
    
    CA_42_1 = rollmean(CA_42_1, window_size, fill = NA, align = "center"),
    CA_42_2 = rollmean(CA_42_2, window_size, fill = NA, align = "center")
    )
  

CpGstepsDRE <- DRE_Region_CpG_methylation %>%
  arrange(desc(start)) %>%  # Ensure sorted order
  transmute(
    AA_28_1 = rollmean(AA_28_1, window_size, fill = NA, align = "center"),
    AA_28_2 = rollmean(AA_28_2, window_size, fill = NA, align = "center"),
    AA_28_3 = rollmean(AA_28_3, window_size, fill = NA, align = "center"),
    CA_28_1 = rollmean(CA_28_1, window_size, fill = NA, align = "center"),
    CA_28_2 = rollmean(CA_28_2, window_size, fill = NA, align = "center"),
    CA_28_3 = rollmean(CA_28_3, window_size, fill = NA, align = "center"),
    
    AA_42_1 = rollmean(AA_42_1, window_size, fill = NA, align = "center"),
    AA_42_2 = rollmean(AA_42_2, window_size, fill = NA, align = "center"),
    
    CA_42_1 = rollmean(CA_42_1, window_size, fill = NA, align = "center"),
    CA_42_2 = rollmean(CA_42_2, window_size, fill = NA, align = "center")
  )

CpGstepsSt01G045040 <- St01G045040_Region_CpG_methylation %>%
  arrange(desc(start)) %>%  # Ensure sorted order
  transmute(
    AA_28_1 = rollmean(AA_28_1, window_size, fill = NA, align = "center"),
    AA_28_2 = rollmean(AA_28_2, window_size, fill = NA, align = "center"),
    AA_28_3 = rollmean(AA_28_3, window_size, fill = NA, align = "center"),
    CA_28_1 = rollmean(CA_28_1, window_size, fill = NA, align = "center"),
    CA_28_2 = rollmean(CA_28_2, window_size, fill = NA, align = "center"),
    CA_28_3 = rollmean(CA_28_3, window_size, fill = NA, align = "center"),
    
    AA_42_1 = rollmean(AA_42_1, window_size, fill = NA, align = "center"),
    AA_42_2 = rollmean(AA_42_2, window_size, fill = NA, align = "center"),
    
    CA_42_1 = rollmean(CA_42_1, window_size, fill = NA, align = "center"),
    CA_42_2 = rollmean(CA_42_2, window_size, fill = NA, align = "center")
  )

CpGstepsAB_Hydrolases <- AB_Hydrolases_Region_CpG_methylation %>%
  arrange(desc(start)) %>%  # Ensure sorted order
  transmute(
    AA_28_1 = rollmean(AA_28_1, window_size, fill = NA, align = "center"),
    AA_28_2 = rollmean(AA_28_2, window_size, fill = NA, align = "center"),
    AA_28_3 = rollmean(AA_28_3, window_size, fill = NA, align = "center"),
    CA_28_1 = rollmean(CA_28_1, window_size, fill = NA, align = "center"),
    CA_28_2 = rollmean(CA_28_2, window_size, fill = NA, align = "center"),
    CA_28_3 = rollmean(CA_28_3, window_size, fill = NA, align = "center"),
    
    AA_42_1 = rollmean(AA_42_1, window_size, fill = NA, align = "center"),
    AA_42_2 = rollmean(AA_42_2, window_size, fill = NA, align = "center"),
    
    CA_42_1 = rollmean(CA_42_1, window_size, fill = NA, align = "center"),
    CA_42_2 = rollmean(CA_42_2, window_size, fill = NA, align = "center")
  )


CpGstepsARR <- St12G001580_Region_CpG_methylation %>%
  arrange(desc(start)) %>%  # Ensure sorted order
  transmute(
    AA_28_1 = rollmean(AA_28_1, window_size, fill = NA, align = "center"),
    AA_28_2 = rollmean(AA_28_2, window_size, fill = NA, align = "center"),
    AA_28_3 = rollmean(AA_28_3, window_size, fill = NA, align = "center"),
    CA_28_1 = rollmean(CA_28_1, window_size, fill = NA, align = "center"),
    CA_28_2 = rollmean(CA_28_2, window_size, fill = NA, align = "center"),
    CA_28_3 = rollmean(CA_28_3, window_size, fill = NA, align = "center"),
    
    AA_42_1 = rollmean(AA_42_1, window_size, fill = NA, align = "center"),
    AA_42_2 = rollmean(AA_42_2, window_size, fill = NA, align = "center"),
    
    CA_42_1 = rollmean(CA_42_1, window_size, fill = NA, align = "center"),
    CA_42_2 = rollmean(CA_42_2, window_size, fill = NA, align = "center")
  )

test1 <- St01G045040_Region_CpG_methylation %>% filter(start> 82887208 & start < 82887869)

test<- CpGstepsARR %>% rownames_to_column("start") %>%
  filter(start > 1338300 & start < 1338600) %>% column_to_rownames("start")

1338365 - 1338526

Heatmap(
  CpGstepsARR %>% as.matrix() %>%t(),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 16,fontface="bold"),
  heatmap_legend_param = list(title = "Methylation",
                              title_gp = gpar(fontsize = 14, fontface = "bold"),  # Increase title size & bold
                              labels_gp = gpar(fontsize = 12, fontface = "bold")),
  width = unit(30, "cm"),
  height = unit(5, "cm")
)


####verify DMRs consistency from 1st and second experiments

ggvenn(
  list(
    Second_28dap= CG_DMR_list_annot_correct$AAvCA_C$UpstreamCorrected$locusName %>% unique(),
    first_42dap = AAvsCA_C_1st$UpstreamCorrected$locusName %>%unique()
  )
)
##Upstream genes.. unique for
anti_join(CG_DMR_list_annot_correct$AAvCA_C$UpstreamCorrected %>% select(locusName) %>% unique(),
          AAvsCA_C_1st$UpstreamCorrected %>% select(locusName) %>% unique())


GetDMRGenesFromAntiJoin<- function(DMR_GeneTable1, DMR_GeneTable2, NameSet1, NameSet2){
  
  ##Get just the gene names, ensure unique, cause some genes may have two dmrs in said region
  GeneID1 <- DMR_GeneTable1 %>% select(locusName) %>% unique()
  GeneID2 <- DMR_GeneTable2 %>% select(locusName) %>% unique()
  
  ##Antijoin to get unqiue sets
  Unique_GeneID1<- anti_join(GeneID1,GeneID2)
  Unique_GeneID2<- anti_join(GeneID2,GeneID1)
  
  ###Get only the unique gene DMRs in each set
  
  DMR_GeneTableUniqueAnti1<- DMR_GeneTable1[DMR_GeneTable1$locusName %in% Unique_GeneID1$locusName,]
  DMR_GeneTableUniqueAnti2<- DMR_GeneTable2[DMR_GeneTable2$locusName %in% Unique_GeneID2$locusName,]
  
  return(setNames(list(DMR_GeneTableUniqueAnti1, 
                       DMR_GeneTableUniqueAnti2), 
                  c(NameSet1, NameSet2)))
  
}
library(ggvenn)
ggvenn(list(
 second= CG_DMR_list_annot_correct$UpstreamCorrected$locusName,
  first= AAvsCA_C_1st$UpstreamCorrected$locusName
))

Antijoin_upstreamDMRGenes<- GetDMRGenesFromAntiJoin(DMR_GeneTable1 = CG_DMR_list_annot_correct$UpstreamCorrected,
                        DMR_GeneTable2 = AAvsCA_C_1st$UpstreamCorrected,
                        NameSet1 = "Second_28dap",
                        NameSet2 = "first_42dap")

###Function to get positions of GrangesObject based on DMRs
GetGrangesPositionFromDMRs <-function(DMR_table,GrangesObject){ 
  
  DMR_Genebody_list <- list()
  for (i in 1:nrow(DMR_table)) {
    
    print(paste0("Getting Granges for DMR ",i))
    current_DMR_row <- DMR_table[i, ]
    
    
    subset_GrangesObject <- 
      GrangesObject[seqnames(GrangesObject) == current_DMR_row$seqnames & 
                      start(GrangesObject) >= current_DMR_row$start.x & 
                      end(GrangesObject) <= current_DMR_row$end.x]
    
    
    DMR_Genebody_list[[i]] <- subset_GrangesObject
    
  }
  return (DMR_Genebody_list)
}


JoinMethylationSamples<- function(DMRcallerObject, number_samples, namestorename){
  
  loops <- number_samples-1
  ##join each sample together
  Combined_DMRcallerObject_GR<-joinReplicates(DMRcallerObject[[1]], DMRcallerObject[[2]]) 
  for (h in 2:loops) { #10 samples hence 10-1
    print(h)
    Combined_DMRcallerObject_GR<-joinReplicates(Combined_DMRcallerObject_GR, DMRcallerObject[[h+1]])
  }
  
  ###Define the column names
  column_names <- namestorename
  
  ##edit such the the meth and cov of each sample will be labelled
  column_names <- rep(column_names, each=2)
  extra_names <- rep(c("meth", "cov"),times=number_samples)
  column_extra_names <- paste(column_names,extra_names, sep ="_")
  
  dummy_end_of_column_to_rename <- (length(column_extra_names) + 2)-1 ##7 is the start of the column to rename
  
  ##rename the column names, not sure how to automate this...
  names(mcols(Combined_DMRcallerObject_GR))[2:dummy_end_of_column_to_rename]<- column_extra_names
  
  
  return(Combined_DMRcallerObject_GR)
}

###join the samples together and rename the samples
methylationDataList_control_single_chr_combined<- JoinMethylationSamples(methylationDataList_control_single_chr,
                       10,
                       namestorename = names(methylationDataList_control_single_chr))

AntijoinUniqueSetsDMRsGenes<- lapply(Antijoin_upstreamDMRGenes, 
       GetGrangesPositionFromDMRs,
       methylationDataList_control_single_chr_combined)

##Calculate percentage fucntion for every CpG positions in the DMRs of every sample
CalculatePercentageCpGPositions<- function(DMRs_list,SampleNames){
  
  GR_DMRslist_percentage <- data.frame()
  for (DMR in seq(DMRs_list)) {
    ##Select the DMRs
    DMR_tables <- as.data.frame(DMRs_list[[DMR]])
    
    region_id <- paste0(DMR_tables$seqnames[1], "_",
                        DMR_tables$start[1], "_",
                        DMR_tables$end[nrow(DMR_tables)])
    
    print(paste0("calculating DMR ", DMR))
    #initialize the percentage df
    DMR_tables_weighted_meth_df <- data.frame(region = region_id)
    for (i in SampleNames) {
      
      ##Select the columns....
      CpG_methylation_Single <- DMR_tables %>% 
        select(contains(i))
      
      ##Remove X from the column names if present
      colnames(CpG_methylation_Single) <- gsub("X","",colnames(CpG_methylation_Single))
      
      
      ####aggregate all regions.. to calculate weighted meth reads and weighted cov
      CpG_methylation_single_Weighted <- colSums(CpG_methylation_Single)
      
      # Rename them with sample prefix
      names(CpG_methylation_single_Weighted) <- c("meth","cov")
      
      ##calculate weighted methylation
      CpG_methylation_single_Weighted_DMR_meth <-  (CpG_methylation_single_Weighted[["meth"]] / CpG_methylation_single_Weighted[["cov"]] ) * 100
      
      # Add as new column with sample-specific name (combing with initialized df)
      DMR_tables_weighted_meth_df[[paste0(i)]] <- CpG_methylation_single_Weighted_DMR_meth
      
    }
    
    
    GR_DMRslist_percentage <- rbind(GR_DMRslist_percentage,
                                    DMR_tables_weighted_meth_df)
    
  }
  return(GR_DMRslist_percentage)
}


test <- AntijoinUniqueSetsDMRsGenes$Second_28dap[[1]]

###get methylaion percentage
AntijoinUniqueSetsDMRsGenes_CpG_perc <- lapply(AntijoinUniqueSetsDMRsGenes, 
                                               CalculatePercentageCpGPositions,
                                               SampleNames=names(methylationDataList_control_single_chr))

test <- AntijoinUniqueSetsDMRsGenes_CpG_perc$Second_28dap[[1]]

###Now have to calculate mean methylation of all CpG positions for all DMRs using colMeans
#CalculateMeanMethylationDMRs<- function(DMRPercentageList){
#  
#  ##Initialize the empty list to store results mean methylation Results
#  DMR_DF <- data.frame()
#  for (DMRs in seq(DMRPercentageList)) {
#    
#    DMR_table_percentage_CpG<- DMRPercentageList[[DMRs]]
#    
#    #Get the chr for the DMRs, is in the first column, row number doesnt matter here
#    DMR_table_percentage_CpG_Chr <- data.frame(chr=as.character(DMR_table_percentage_CpG[1,1]))
#    
#    ##Calculate the Mean methylation of all samples
#    DMR_table_percentage_mean <- t(as.data.frame(colMeans(DMR_table_percentage_CpG[,8:ncol(DMR_table_percentage_CpG)])))
#    rownames(DMR_table_percentage_mean) <- DMRs
#    
#    ##Get the DMR length that was based on the parents
#    DMR_length <- data.frame(start = DMR_table_percentage_CpG[1,2],
#                             end = DMR_table_percentage_CpG[nrow(DMR_table_percentage_CpG),2])
#    
#    #Combine the information together to a single row DF for that DMR
#    DMR_mean_methylation <- cbind(DMR_table_percentage_CpG_Chr,DMR_length,DMR_table_percentage_mean)
#    
#    ##Append to the empty dataframe to store mean methylation results for all DMR regions
#    DMR_DF <- rbind(DMR_DF,DMR_mean_methylation)
#    
#  }
#  
#  return(DMR_DF)
#  
#}


#AntijoinUniqueSetsDMRsGenes_CpG_perc_mean<- lapply(AntijoinUniqueSetsDMRsGenes_CpG_perc, 
#       CalculateMeanMethylationDMRs)

##Create metadata for labels
MetaInfo <- data.frame(
  SampleID= colnames(AntijoinUniqueSetsDMRsGenes_CpG_perc$Second_28dap)[4:ncol(AntijoinUniqueSetsDMRsGenes_CpG_perc$Second_28dap)]
)

MetaInfo <- MetaInfo %>%
  as.data.frame() %>%
  mutate(Genotypes=case_when(str_detect(SampleID,"AA") ~ "AA",
                             str_detect(SampleID,"CA") ~ "CA"),
         TimePoint= case_when(str_detect(SampleID,"28") ~ "Control_28dap",
                               str_detect(SampleID,"42") ~ "Control_42dap"))


library(factoextra)
library(FactoMineR)
PCA_DMR_second28dap<- PCA(t(na.omit(AntijoinUniqueSetsDMRsGenes_CpG_perc_mean$Second_28dap[4:ncol(AntijoinUniqueSetsDMRsGenes_CpG_perc_mean$Second_28dap)])),
    graph = FALSE,
    scale.unit = TRUE)

fviz_pca_ind(PCA_DMR_second28dap,
             # Fill individuals by groups
             geom.ind = "point", addEllipses = FALSE,
             col.var = "black",
             repel = TRUE) + geom_point(aes(shape = factor(MetaInfo$TimePoint),
                                            colour = factor(MetaInfo$Genotypes)),size = 6) +
  guides(colour=guide_legend("Genotypes"),
         shape = guide_legend("Conditions")) +
  ggtitle("CpG methylation percentage (1st and 2nd experiment)") +
  geom_text(label = MetaInfo$SampleID, vjust=1.7)+
  scale_colour_manual(values=c("AA" = "#CC0000",
                               "CA" ="lightblue")) + theme(
                                 legend.title = element_text(face = "bold",size = 15),
                                 legend.text = element_text(face = "bold",size = 15),
                                 axis.title = element_text(face = "bold",size = 15),
                                 axis.text = element_text(face = "bold",size=12))

library(Rtsne)
tsne_out <- Rtsne(PCA_DMR_second28dap$ind$coord, perplexity = 3)

ggplot(as.data.frame(tsne_out$Y), aes(x=V1, y=V2)) +
  geom_point(aes(shape = factor(MetaInfo$TimePoint),
                 colour = factor(MetaInfo$Genotypes)),size = 6) +
  guides(colour=guide_legend("Genotypes"),
         shape = guide_legend("Conditions")) +
  scale_colour_manual(values=c("AA" = "#CC0000",
                               "CA" ="lightblue")) +
  xlab("tSNE1")+
  ylab("tSNE2")+
  theme_bw()+
  theme(
    legend.title = element_text(face = "bold",size = 12),
    legend.text = element_text(face = "bold",size = 10),
    axis.title = element_text(face = "bold",size = 12),
    axis.text = element_text(face = "bold"))


Second_28dap_DMR_mean<- AntijoinUniqueSetsDMRsGenes_CpG_perc_mean$Second_28dap[4:ncol(AntijoinUniqueSetsDMRsGenes_CpG_perc_mean$Second_28dap)]
First_42dap_DMR_mean<- AntijoinUniqueSetsDMRsGenes_CpG_perc_mean$first_42dap[4:ncol(AntijoinUniqueSetsDMRsGenes_CpG_perc_mean$first_42dap)]


rownames(AntijoinUniqueSetsDMRsGenes_CpG_perc$Second_28dap) <- NULL
rownames(AntijoinUniqueSetsDMRsGenes_CpG_perc$first_42dap) <- NULL
ComplexHeatmap::Heatmap(
  as.matrix(AntijoinUniqueSetsDMRsGenes_CpG_perc$first_42dap %>% unique() %>% select(-region)),
  cluster_columns = F,
  cluster_rows =T,
  heatmap_legend_param = list(
    title = "Methylation"
  ),
  show_row_names = F
)

ComplexHeatmap::Heatmap(
  as.matrix(na.omit(First_42dap_DMR_mean)),
  cluster_columns = FALSE,
  cluster_rows =T,
  heatmap_legend_param = list(
    title = "Methylation"
  ),
  show_row_names = F
)





########Compare methylation differences between two experiments######






methylationDataList_1st2nd_pooled <- GRangesList(
  ###PArents
  "AA_2nd_pooled" = readBismarkPool(c(paste0(path2,"CG_CHR_Qualtrim_AA4K_1_bismark_bt2_pe.deduplicated.CX_report.txt"), 
                                      paste0(path2,"CG_CHR_Qualtrim_AA5K_1_bismark_bt2_pe.deduplicated.CX_report.txt"),
                                      paste0(path2,"CG_CHR_Qualtrim_AA15K_1_bismark_bt2_pe.deduplicated.CX_report.txt"))),
  
  "CA_2nd_pooled" = readBismarkPool(c(paste0(path2,"CG_CHR_Qualtrim_CA2K_1_bismark_bt2_pe.deduplicated.CX_report.txt"), 
                                      paste0(path2,"CG_CHR_Qualtrim_CA3_1_bismark_bt2_pe.deduplicated.CX_report.txt"),
                                      paste0(path2,"CG_CHR_Qualtrim_CA11_1_bismark_bt2_pe.deduplicated.CX_report.txt"))),
  
  "AA_1st_pooled" = readBismarkPool(c(paste0(path,"_AA_A_1_bismark_bt2_pe.ded_CG.txt"), 
                                      paste0(path,"_AA_B_1_bismark_bt2_pe.ded_CG.txt"))),
  
  "CA_1st_pooled" = readBismarkPool(c(paste0(path,"_CA_A_1_bismark_bt2_pe.ded_CG.txt"), 
                                      paste0(path,"_CA_B_1_bismark_bt2_pe.ded_CG.txt")))

)

methylationDataList_1st2nd_pooled_chr <-lapply(methylationDataList_1st2nd_pooled, 
                                        filter_scaffolds) 

methylationDataList_1st2nd_pooled_chr_joined<- JoinMethylationSamples(methylationDataList_1st2nd_pooled_chr,
                                                               number_samples=4,
                                                               namestorename = names(methylationDataList_1st2nd_pooled_chr))


AntijoinUniqueSetsDMRsGenes_pooled<- lapply(Antijoin_upstreamDMRGenes, 
                                     GetGrangesPositionFromDMRs,
                                     methylationDataList_1st2nd_pooled_chr_joined)


AntijoinUniqueSetsDMRsGenes_pooled_CpG<- lapply(AntijoinUniqueSetsDMRsGenes_pooled, 
       CalculatePercentageCpGPositionsForEachDMR,
       names(methylationDataList_1st2nd_pooled_chr))

AntijoinUniqueSetsDMRsGenes_pooled_CpG_mean<- lapply(AntijoinUniqueSetsDMRsGenes_pooled_CpG, 
                                                     CalculateMeanMethylationDMRs)

AntijoinUniqueSetsDMRsGenes_ProportionDifferenceCpG <- list()
for (SecondFirst in names(AntijoinUniqueSetsDMRsGenes_pooled_CpG_mean)) {
  
  AntijoinUniqueSetsDMRsGenes_ProportionDifferenceCpG[[SecondFirst]]<- 
    AntijoinUniqueSetsDMRsGenes_pooled_CpG_mean[[SecondFirst]] %>%
    transmute(
      chr=chr,
      start=start,
      end=end,
      AA_CA_2nd = CA_2nd_pooled/AA_2nd_pooled,
      AA_CA_1st= CA_1st_pooled/AA_1st_pooled
    )
}

AntijoinUniqueSetsDMRsGenes_ProportionDifferenceCpG_ABS_DIFF <- list()
for (SecondFirst in names(AntijoinUniqueSetsDMRsGenes_pooled_CpG_mean)) {
  
  AntijoinUniqueSetsDMRsGenes_ProportionDifferenceCpG_ABS_DIFF[[SecondFirst]]<- 
    AntijoinUniqueSetsDMRsGenes_pooled_CpG_mean[[SecondFirst]] %>%
    transmute(
      chr=chr,
      start=start,
      end=end,
      AA_CA_2nd = abs(CA_2nd_pooled - AA_2nd_pooled),
      AA_CA_1st= abs(CA_1st_pooled - AA_1st_pooled)
    ) %>% pivot_longer(
                       cols = c("AA_CA_2nd", "AA_CA_1st"), 
                       names_to = "Experiment", 
                       values_to = "Abs_Difference")
  
}

ggplot(AntijoinUniqueSetsDMRsGenes_ProportionDifferenceCpG_ABS_DIFF$Second_28dap, 
       aes(x = Experiment, y = Abs_Difference, fill = Experiment)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.4, linetype = "dashed", color = "red") +
  labs(title = "Comparison of Absolute Differences Between AA and CA (DMRs detected in 2nd Exp)",
       x = "Experiment", y = "Absolute Difference (AA - CA)") +
  scale_fill_manual(values = c("blue4", "red4")) +
  theme_bw() +theme(
    axis.title = element_text(size = 12,face = "bold"),
    axis.text = element_text(size = 12,face = "bold"),
    legend.title = element_text(size = 12,face = "bold"),
    legend.text = element_text(size = 12,face = "bold")
  )



# Scatter plot of fold changes
test <- na.omit(AntijoinUniqueSetsDMRsGenes_ProportionDifferenceCpG$first_42dap)
test$AA_CA_2nd <- log10(test$AA_CA_2nd+1)
test$AA_CA_1st <- log10(test$AA_CA_1st+1)

ggplot(test, 
       aes(x = AA_CA_1st, y = AA_CA_2nd)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", linetype = "solid") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +  # Add identity line
  labs(title = "Fold Change Comparison between Two Experiments",
       x = "AA_CA_1st Experiment (log10 FC)",
       y = "AA_CA_2nd Experiment (log10 FC)") +
  theme_minimal()

###Combined analysis for coverage
path2 <- "/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/WGBS_HIGH_HEAT/03_BISMARK/CX_FILES/Chr_Context/"
path <- "/media/rna/Epipotato16TB/Epipotato_methylome/04_BISMARK/METHYLATION_CX_REPORT/"

methylationDataList_control_single <- GRangesList(
  ###PArents
  "AA_28_1" = readBismark(paste0(path2,"CG_CHR_Qualtrim_AA4K_1_bismark_bt2_pe.deduplicated.CX_report.txt")), 
  "AA_28_2" = readBismark(paste0(path2,"CG_CHR_Qualtrim_AA5K_1_bismark_bt2_pe.deduplicated.CX_report.txt")),
  "AA_28_3" = readBismark(paste0(path2,"CG_CHR_Qualtrim_AA15K_1_bismark_bt2_pe.deduplicated.CX_report.txt")), 
  
  "CA_28_1" = readBismark(paste0(path2,"CG_CHR_Qualtrim_CA2K_1_bismark_bt2_pe.deduplicated.CX_report.txt")),
  "CA_28_2" = readBismark(paste0(path2,"CG_CHR_Qualtrim_CA3_1_bismark_bt2_pe.deduplicated.CX_report.txt")),
  "CA_28_3" = readBismark(paste0(path2,"CG_CHR_Qualtrim_CA11_1_bismark_bt2_pe.deduplicated.CX_report.txt")),
  
  "AA_42_1" = readBismark(paste0(path,"_AA_A_1_bismark_bt2_pe.ded_CG.txt")), 
  "AA_42_2" = readBismark(paste0(path,"_AA_B_1_bismark_bt2_pe.ded_CG.txt")),
  
  "CA_42_1" = readBismark(paste0(path,"_CA_A_1_bismark_bt2_pe.ded_CG.txt")), 
  "CA_42_2" = readBismark(paste0(path,"_CA_B_1_bismark_bt2_pe.ded_CG.txt")))

####FOR CHG
methylationDataList__CHG_control_single <- GRangesList(
  ###PArents
  "AA_28_1" = readBismark(paste0(path2,"CHG_CHR_Qualtrim_AA4K_1_bismark_bt2_pe.deduplicated.CX_report.txt")), 
  "AA_28_2" = readBismark(paste0(path2,"CHG_CHR_Qualtrim_AA5K_1_bismark_bt2_pe.deduplicated.CX_report.txt")),
  "AA_28_3" = readBismark(paste0(path2,"CHG_CHR_Qualtrim_AA15K_1_bismark_bt2_pe.deduplicated.CX_report.txt")), 
  
  "CA_28_1" = readBismark(paste0(path2,"CHG_CHR_Qualtrim_CA2K_1_bismark_bt2_pe.deduplicated.CX_report.txt")),
  "CA_28_2" = readBismark(paste0(path2,"CHG_CHR_Qualtrim_CA3_1_bismark_bt2_pe.deduplicated.CX_report.txt")),
  "CA_28_3" = readBismark(paste0(path2,"CHG_CHR_Qualtrim_CA11_1_bismark_bt2_pe.deduplicated.CX_report.txt")),
  
  "AA_42_1" = readBismark(paste0(path,"_AA_A_1_bismark_bt2_pe.ded_CHG.txt")), 
  "AA_42_2" = readBismark(paste0(path,"_AA_B_1_bismark_bt2_pe.ded_CHG.txt")),
  
  "CA_42_1" = readBismark(paste0(path,"_CA_A_1_bismark_bt2_pe.ded_CHG.txt")), 
  "CA_42_2" = readBismark(paste0(path,"_CA_B_1_bismark_bt2_pe.ded_CHG.txt")))


methylationDataList_pooled <- GRangesList(
  ###PArents
  "AA_2nd_pooled" = readBismarkPool(c(paste0(path2,"CG_CHR_Qualtrim_AA4K_1_bismark_bt2_pe.deduplicated.CX_report.txt"), 
                                        paste0(path2,"CG_CHR_Qualtrim_AA5K_1_bismark_bt2_pe.deduplicated.CX_report.txt"),
                                        paste0(path2,"CG_CHR_Qualtrim_AA15K_1_bismark_bt2_pe.deduplicated.CX_report.txt"))),
  
  "CA_2nd_pooled" = readBismarkPool(c(paste0(path2,"CG_CHR_Qualtrim_CA2K_1_bismark_bt2_pe.deduplicated.CX_report.txt"), 
                                        paste0(path2,"CG_CHR_Qualtrim_CA3_1_bismark_bt2_pe.deduplicated.CX_report.txt"),
                                        paste0(path2,"CG_CHR_Qualtrim_CA11_1_bismark_bt2_pe.deduplicated.CX_report.txt"))),
  
  "AA_2nd_H_pooled" = readBismarkPool(c(paste0(path2,"CG_CHR_Qualtrim_AA43W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt"), 
                                          paste0(path2,"CG_CHR_Qualtrim_AA53W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt"),
                                          paste0(path2,"CG_CHR_Qualtrim_AA153W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt"))),
  
  "CA_2nd_H_pooled" = readBismarkPool(c(paste0(path2,"CG_CHR_Qualtrim_CA23W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt"), 
                                          paste0(path2,"CG_CHR_Qualtrim_CA33W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt"),
                                          paste0(path2,"CG_CHR_Qualtrim_CA113W37C_1_bismark_bt2_pe.deduplicated.CX_report.txt"))),
  
  
  "AA_1st_pooled" = readBismarkPool(c(paste0(path,"_AA_A_1_bismark_bt2_pe.ded_CG.txt"), 
                                        paste0(path,"_AA_B_1_bismark_bt2_pe.ded_CG.txt"))),
  
  "CA_1st_pooled" = readBismarkPool(c(paste0(path,"_CA_A_1_bismark_bt2_pe.ded_CG.txt"), 
                                          paste0(path,"_CA_B_1_bismark_bt2_pe.ded_CG.txt"))),
  
  "AA_1st__H_pooled" = readBismarkPool(c(paste0(path,"_AA_F_1_bismark_bt2_pe.ded_CG.txt"), 
                                        paste0(path,"_AA_G_1_bismark_bt2_pe.ded_CG.txt"))),
  
  "CA_1st__H_pooled" = readBismarkPool(c(paste0(path,"_CA_F_1_bismark_bt2_pe.ded_CG.txt"), 
                                      paste0(path,"_CA_G_1_bismark_bt2_pe.ded_CG.txt"))),
  
  
)

methylationDataList_pooled_chr <-lapply(methylationDataList_pooled, 
                                        filter_scaffolds) 

methylationDataList_control_single_chr <-lapply(methylationDataList_control_single, 
                                        filter_scaffolds) 
##CHG
methylationDataList__CHG_control_single_chr<-lapply(methylationDataList__CHG_control_single, 
                                                    filter_scaffolds) 


#####GET COVERAGE...
methylationDataList_pooled_chr_joined<- JoinMethylationSamples(methylationDataList_pooled_chr,
                       number_samples=8,
                       namestorename = names(methylationDataList_pooled_chr))

methylationDataList_pooled_chr_joined_COV_DF <-as.data.frame(methylationDataList_pooled_chr_joined) %>% 
  select(seqnames,start,end,context,contains("cov"))



tail(methylationDataList_pooled_chr_joined_COV_DF)[,5:8]


methylationDataList_pooled_chr_joined_COV_DF_long <- methylationDataList_pooled_chr_joined_COV_DF %>%
  pivot_longer(cols = starts_with("AA") | starts_with("CA"), 
               names_to = "sample", 
               values_to = "coverage")

methylationDataList_pooled_chr_joined_COV_DF_long$coverage <- log10(methylationDataList_pooled_chr_joined_COV_DF_long$coverage +1)

ggplot(methylationDataList_pooled_chr_joined_COV_DF_long, 
       aes(x = coverage, fill = sample)) +
  geom_histogram(position = "dodge", binwidth = 0.4, color = "black", alpha = 0.3,
                 bins = 30) +
  labs(title = "Coverage Distribution", x = "Coverage", y = "Frequency") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1")


coverageCG_AA_28 <- computeMethylationDataCoverage(methylationDataList_pooled_chr[["AA_28dap_pooled"]],
                                               context="CG",
                                               breaks = c(1,5,10,15,20,25,30,35,40))

Coverage_meth_GenomeWidelist<- list(
  AA_28dap_pooled =computeMethylationDataCoverage(methylationDataList_pooled_chr[["CA_28dap_pooled"]],
                                                     context="CG",
                                                     breaks = c(1,5,10,15,20,25,30,35,40,45,50,55,60)),
  CA_28dap_pooled =computeMethylationDataCoverage(methylationDataList_pooled_chr[["AA_28dap_pooled"]],
                                                  context="CG",
                                                  breaks = c(1,5,10,15,20,25,30,35,40,45,50,55,60)),
  AA_42dap_pooled =computeMethylationDataCoverage(methylationDataList_pooled_chr[["AA_42dap_pooled"]],
                                                  context="CG",
                                                  breaks = c(1,5,10,15,20,25,30,35,40,45,50,55,60)),
  CA_42dap_pooled =computeMethylationDataCoverage(methylationDataList_pooled_chr[["CA_42dap_pooled"]],
                                                  context="CG",
                                                  breaks = c(1,5,10,15,20,25,30,35,40,45,50,55,60))
)
methylationDataList_control_single_chr

CalculateCoverageMethylationData <- function(methylationData,
                                             CONTEXT,
                                             BreakPoint=c(1,5,10,15,20,25)){
  
  ##calculate proportion of cytosines for a given read count threshold 
 CoverageVector <- computeMethylationDataCoverage(
    methylationData,
    context=CONTEXT,
    breaks=BreakPoint)
 
 
 return(CoverageVector)
  
}

Coverage_meth_GenomeWide_Single_list <- list()
for (samples in names(methylationDataList_control_single_chr)) {
  
  print(paste0("Calculating coverage for ", samples))
  
  Coverage_meth_GenomeWide_Single_list[[samples]]<- computeMethylationDataCoverage(
    methylationDataList_control_single_chr[[samples]], 
    context="CG", breaks=c(1,5,10,15,20,25,30,35,40,45,50))
  
}

methylationDataList__CHG_control_single_chr
Coverage_meth_CHG_GenomeWide_Single_list <- list()
for (samples in names(methylationDataList__CHG_control_single_chr)) {
  
  print(paste0("Calculating coverage for ", samples))
  
  Coverage_meth_CHG_GenomeWide_Single_list[[samples]]<- computeMethylationDataCoverage(
    methylationDataList__CHG_control_single_chr[[samples]], 
    context="CHG", breaks=c(1,5,10,15,20,25,30,35,40,45,50))
  
}

#Coverage_meth_GenomeWide_Single_list<- lapply(methylationDataList_control_single_chr, 
#       CalculateCoverageMethylationData,
#       CONTEXT="CG",
#       BreakPoint=c(1,5,10,15,20,25,30,35,40,45,50))

Coverage_meth_GenomeWide_DF<- as.data.frame(Coverage_meth_GenomeWidelist)
Coverage_meth_GenomeWide_DF$MinimumCounts <- c(1,5,10,15,20,25,30,35,40,45,50,55,60)

Coverage_meth_GenomeWide_DF_long <- Coverage_meth_GenomeWide_DF %>%
  pivot_longer(cols = starts_with("AA") | starts_with("CA"), 
               names_to = "sample", 
               values_to = "proportion")

# Plot using ggplot2
ggplot(Coverage_meth_GenomeWide_DF_long, 
       aes(x = MinimumCounts, y = proportion, color = sample, group = sample)) +
  geom_line(size = 1) + 
  geom_point(size = 2) +
  labs(title = "Proportion of Cytosines Covered at Different Read Thresholds",
       x = "Coverage Threshold (Reads)", 
       y = "Proportion of Cytosines Covered") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")


Coverage_meth_GenomeWide_Single_DF<- as.data.frame(Coverage_meth_GenomeWide_Single_list)
Coverage_meth_GenomeWide_Single_DF$MinimumCounts <- c(1,5,10,15,20,25,30,35,40,45,50)

##CHG
Coverage_meth_CHG_GenomeWide_Single_DF <- as.data.frame(Coverage_meth_CHG_GenomeWide_Single_list)
Coverage_meth_CHG_GenomeWide_Single_DF$MinimumCounts <- c(1,5,10,15,20,25,30,35,40,45,50)

Coverage_meth_GenomeWide_Single_DF_long <- pivot_longer(Coverage_meth_GenomeWide_Single_DF, cols = starts_with("AA") | starts_with("CA"),
                                   names_to = "Sample", values_to = "Proportion") %>%
                                      mutate(
                                     TimePoint= case_when(str_detect(Sample,"28") ~ "Control_28dap",
                                                          str_detect(Sample,"42") ~ "Control_42dap"))

##CHG
Coverage_meth_CHG_GenomeWide_Single_DF_long <- pivot_longer(Coverage_meth_CHG_GenomeWide_Single_DF, cols = starts_with("AA") | starts_with("CA"),
                                                        names_to = "Sample", values_to = "Proportion") %>%
  mutate(
    TimePoint= case_when(str_detect(Sample,"28") ~ "Control_28dap",
                         str_detect(Sample,"42") ~ "Control_42dap"))

Coverage_meth_GenomeWide_Single_DF_long$MinimumCounts<- as.factor(Coverage_meth_GenomeWide_Single_DF_long$MinimumCounts)

#CHG
Coverage_meth_CHG_GenomeWide_Single_DF_long$MinimumCounts <- as.factor(Coverage_meth_CHG_GenomeWide_Single_DF_long$MinimumCounts)

# Reshape data to long format for plotting
ggplot(Coverage_meth_CHG_GenomeWide_Single_DF_long, 
       aes(x = MinimumCounts, y = Proportion, color = TimePoint)) +
geom_boxplot()+
  labs(title = "Proportion of Cytosines Covered at Different Read Thresholds",
       x = "Coverage Threshold (Reads)", 
       y = "Proportion of Cytosines Covered") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")+theme(
    axis.title = element_text(size = 14,face = "bold"),
    axis.text = element_text(size = 12,face = "bold"),
    legend.text = element_text(size = 12,face = "bold"),
    legend.title  = element_text(size = 14,face = "bold"),
  )


#####TRY to COMPUTE COVERAGE FOR only identified DMRS sites...
Antijoin_upstreamDMRGenes$Second_28dap
methylationDataList_control_single_chr$AA_28_1

###For DMRs detected in the second experiment  28dap (controls conditions)
DMR_regions <- GRanges(
  seqnames = Antijoin_upstreamDMRGenes$first_42dap$seqnames,
  ranges = IRanges(start = Antijoin_upstreamDMRGenes$first_42dap$start.x,
                   end = Antijoin_upstreamDMRGenes$first_42dap$end.x)
)

####subset the regions by the dmr regions identified in second experiment###
MethSubset_control_single_DMRRegions<- lapply(methylationDataList_control_single_chr, 
       subsetByOverlaps,
       DMR_regions)


MethSubset_control_single_DMRRegions_cov<- lapply(MethSubset_control_single_DMRRegions, 
                                              computeMethylationDataCoverage,
                                              context = "CG",
                                              breaks = c(1,5,10,15,20,25,30,35,40,45,50))

MethSubset_control_single_DMRRegions_covDF <- as.data.frame(MethSubset_control_single_DMRRegions_cov)
MethSubset_control_single_DMRRegions_covDF$MinimumCounts <- c(1,5,10,15,20,25,30,35,40,45,50)


MethSubset_control_single_DMRRegions_covDF_long <- pivot_longer(MethSubset_control_single_DMRRegions_covDF, cols = starts_with("AA") | starts_with("CA"),
                                                        names_to = "Sample", values_to = "Proportion") %>%
  mutate(
    TimePoint= case_when(str_detect(Sample,"28") ~ "Control_28dap",
                         str_detect(Sample,"42") ~ "Control_42dap"))



# Reshape data to long format for plotting
ggplot(MethSubset_control_single_DMRRegions_covDF_long, 
       aes(x = as.factor(MinimumCounts), y = Proportion, color = TimePoint)) +
  geom_boxplot()+
  labs(title = "Proportion of Cytosines Covered at Different Read Thresholds",
       x = "Coverage Threshold (Reads)", 
       y = "Proportion of Cytosines Covered") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")+theme(
    axis.title = element_text(size = 14,face = "bold"),
    axis.text = element_text(size = 12,face = "bold"),
    legend.text = element_text(size = 12,face = "bold"),
    legend.title  = element_text(size = 14,face = "bold"),
  )

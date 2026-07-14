###Do enrichment analyses (both GO and KO) for lists of genes directly from dataframe of DESeq2
GO_KO_enrichment_with_enricher_CP <- function(DEGlists, term2gene_GO, term2name_GO, pvalue = 0.05){
  
  Tolerance_DEG_enrichment_list <- list()
  
  for (i in names(DEGlists)) {
    print(paste("Processing gene vector:", i))
    
    enrichment_intersections_GO<- enricher(
      gene = DEGlists[[i]],
      pvalueCutoff = pvalue,
      pAdjustMethod = "fdr",
      universe = NULL,
      minGSSize = 10,
      maxGSSize = 500,
      qvalueCutoff = 0.2,
      gson = NULL,
      TERM2GENE = term2gene_GO,
      TERM2NAME = term2name_GO
    )
    
#  enrichment_intersections_KO<- enricher(
#    gene =  DEGlists[[i]]$locusName,
#    pvalueCutoff = pvalue,
#    pAdjustMethod = "BH",
#    universe = NULL,
#    minGSSize = 10,
#    maxGSSize = 500,
#    qvalueCutoff = 0.2,
#    gson = NULL,
#    TERM2GENE = term2gene_KO,
#    TERM2NAME = term2name_KO
#  )
#  
    Tolerance_DEG_enrichment_list[[i]] <- enrichment_intersections_GO
    
  }
  
  return(Tolerance_DEG_enrichment_list)
  
}


CombineUpDownRegulation_SelectTopEnrichmentTermsByCounts_Plot <- function(TopCounts=20, 
                                                                     Upregulated_DF,
                                                                     Downregulated_DF){
  CombineRegulation_DF <- rbind(Upregulated_DF %>% 
                   mutate(Regulation="Upregulated") %>% 
                   select(Description,GeneRatio,Count,p.adjust,Regulation,geneID) %>% 
                   arrange(desc(Count), .by_group =TRUE) %>% dplyr::filter(p.adjust < 0.05) %>%
                   slice(1:TopCounts),
                 Downregulated_DF %>% 
                   mutate(Regulation="Downregulated") %>% 
                   select(Description,GeneRatio,Count,p.adjust,Regulation,geneID) %>% 
                   arrange(desc(Count), .by_group =TRUE) %>% dplyr::filter(p.adjust < 0.05) %>%
                   slice(1:TopCounts))
  
 split_elements <- strsplit(CombineRegulation_DF$GeneRatio, "/")
  
 numerators <- sapply(split_elements, function(x) as.numeric(x[1]))
 denominators <- sapply(split_elements, function(x) as.numeric(x[2]))
 
 
 CombineRegulation_DF <- CombineRegulation_DF %>% 
   mutate(GeneRatioNumeric = as.numeric(numerators/denominators)) %>%
   select(Description,GeneRatio,GeneRatioNumeric,Count,p.adjust,Regulation,geneID)
  
  return(CombineRegulation_DF)
  
}

Check_for_overlaps_between_two_GENE_list <- function(DF1,DF2){
  
  #Overlaps_between_two_DFs
  Overlaps_between_two_DFs <- inner_join(
    data.frame(locusName=DF1$locusName),
    data.frame(locusName=DF2$locusName))
  
  #DF1 unique genes with overlaps from DF2
  DF1withOverlapsFromDF2 <- left_join(data.frame(locusName=DF1$locusName),
                                      data.frame(locusName=DF2$locusName))
                                                         
  #DF2 unique genes with overlaps from DF1
  DF2withOverlapsFromDF1 <- left_join(data.frame(locusName=DF2$locusName),
                                      data.frame(locusName=DF1$locusName))
  
  #Create binaries 
  DF1withOverlapsFromDF2 <- DF1withOverlapsFromDF2 %>% 
    mutate(
      overlaps=case_when(str_detect(locusName,
                                    paste(Overlaps_between_two_DFs$locusName, 
                                          collapse = "|")) ~ 1, TRUE ~ 0)) %>% arrange(desc(overlaps))
  
  #Create binaries 
  DF2withOverlapsFromDF1 <- DF2withOverlapsFromDF1 %>% 
    mutate(
      overlaps=case_when(str_detect(locusName,
                                    paste(Overlaps_between_two_DFs$locusName, 
                                          collapse = "|")) ~ 1, TRUE ~ 0)) %>% arrange(desc(overlaps))
  
  #Only works for me, :)
  DF1withOverlapsFromDF2 <- unique(left_join(DF1withOverlapsFromDF2,S.tuberosum.v6.1_latest))
  DF2withOverlapsFromDF1 <- unique(left_join(DF2withOverlapsFromDF1,S.tuberosum.v6.1_latest))
  
  
  OverlapsDF1_DF2 <- list(
    DF1withOverlapsFromDF2,
    DF2withOverlapsFromDF1 
  )
  
  return(OverlapsDF1_DF2)
  
}


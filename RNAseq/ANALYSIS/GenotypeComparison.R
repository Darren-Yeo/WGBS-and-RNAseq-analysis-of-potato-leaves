library(tidyverse)
library(DESeq2)
library(umap)
library(EnhancedVolcano)
library(gridExtra)
library(clusterProfiler)
library(matrixStats)
library(dendextend)
library(factoextra)
library(FactoMineR)
library(cluster)
library(circlize)
library(ggvenn)


GeneID2GeneLabel<- S.tuberosum.v6.1_latest %>% dplyr::select(locusName,v6.1.Description,aktualisierte.annotation..bitte.ergänzen.,Mapman.description)

#Import unique raw counts....
UniqueCountsDuplicatedHighStress <- read.delim("/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/05_FEATURECOUNTS/UniqueCountsDuplicatedHighStress.txt", comment.char="#")
#Import Meta Data

GeneInfo <- UniqueCountsDuplicatedHighStress[,1:6]

#Get only the raw counts
UniqueCountsDuplicatedHighStress_counts <- UniqueCountsDuplicatedHighStress %>%
  column_to_rownames("Geneid") %>% 
  dplyr::select(-colnames(GeneInfo)[2:6])

colnames(UniqueCountsDuplicatedHighStress_counts) <- substring(colnames(UniqueCountsDuplicatedHighStress_counts),
                                                               1, 
                                                               nchar(colnames(UniqueCountsDuplicatedHighStress_counts)) - 35)



MetaData<- data.frame(SampleID=colnames(UniqueCountsDuplicatedHighStress_counts))

MetaData<- MetaData%>% 
  mutate(Conditions=case_when(str_detect(SampleID,"K") ~ "Control",
                              str_detect(SampleID,"H") ~ "Heat"),
         Genotype=case_when(str_detect(SampleID,"2K") ~ "CA",
                            str_detect(SampleID,"3K") ~ "CA",
                            str_detect(SampleID,"11K") ~ "CA",
                            str_detect(SampleID,"15K") ~ "AA",
                            str_detect(SampleID,"4K") ~ "AA",
                            str_detect(SampleID,"5K") ~ "AA",
                            str_detect(SampleID,"3H") ~ "CA",
                            str_detect(SampleID,"11H") ~ "CA",
                            str_detect(SampleID,"2H") ~ "CA",
                            str_detect(SampleID,"4H") ~ "AA",
                            str_detect(SampleID,"5H") ~ "AA",
                            str_detect(SampleID,"15H") ~ "AA"),
         Timepoint=case_when(str_detect(SampleID,"H") ~ "48dap",
                             str_detect(SampleID,"K") ~ "28dap"))

###get only respective samples for the comparison
UniqueCountsDuplicatedHighStress_counts_Control<- UniqueCountsDuplicatedHighStress_counts %>% dplyr::select(contains("K"))
UniqueCountsDuplicatedHighStress_counts_Heat<- UniqueCountsDuplicatedHighStress_counts %>% dplyr::select(contains("H"))

###Get respective data for the metadata
MetaData_Control <- MetaData %>% dplyr::filter(str_detect(Conditions,"Control"))
MetaData_Heat <- MetaData %>% dplyr::filter(str_detect(Conditions,"Heat"))
###Analysis with interacting terms.. conditions as the main comparison

dds_Genotype_list<- list(
    control = DESeqDataSetFromMatrix(countData=UniqueCountsDuplicatedHighStress_counts_Control, 
                               colData=MetaData_Control, 
                               design= ~ Genotype),

    Heat = DESeqDataSetFromMatrix(countData=UniqueCountsDuplicatedHighStress_counts_Heat, 
                                                 colData=MetaData_Heat, 
                                                 design= ~ Genotype) )


##Get results from the comparisons
dds_Genotype_res_list <- list()
for (dds_name in names(dds_Genotype_list)) {
  
  #run the DEseq function
  dds_Genotype <- DESeq(dds_Genotype_list[[dds_name]])
  
  #get the results
  
  res <- results(dds_Genotype, contrast = c("Genotype","AA","CA"), pAdjustMethod = "fdr")
  
  rownames(res) <- gsub(".v6.1","",rownames(res))
  
  ##label the genes
  dds_Genotype_res_list[[dds_name]] <- res %>% 
    as.data.frame() %>% 
    rownames_to_column("locusName") %>%
    inner_join(GeneID2GeneLabel)
  
}

ggsave("VolcanoAAvsCA_Control.svg",
       EnhancedVolcano::EnhancedVolcano(dds_Genotype_res_list$control,
                                        lab = dds_Genotype_res_list$control$aktualisierte.annotation..bitte.ergänzen.,
                                        selectLab = c("StPEBP9, SP6A","StPEBP4, SP3D"),
                                        col=c('dimgrey', 'dimgrey', 'dimgrey', '#C83B38'),
                                        x = 'log2FoldChange',
                                        y = 'padj',
                                        pCutoff =  0.05,
                                        FCcutoff = 1,
                                        labSize = 7,
                                        drawConnectors = TRUE,
                                        widthConnectors = 1.3,
                                        typeConnectors= "open",
                                        lengthConnectors = unit(0.008, "npc"),
                                        title = "") + 
         theme(
           axis.title = element_text(size=25,face = "bold"),  # Bold x and y axis titles
           axis.text = element_text(size=16,face = "bold"),
           legend.text =  element_text(size=16),
           legend.title  =  element_text(size=16)
         ),
       device = "svg",
       path = "/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/06_ANALYSIS/EpiPotato_mRNA_HighStress/AnnabelleVsCamel",
       width = 10,
       height = 10,
       dpi = 300
       )





###Get sig genes up and downregulated for genotype comparison
dds_Genotype_res_sig_list <- list()
for (dds_name in names(dds_Genotype_res_list)) {
  
  dds_Genotype_res_sig_list[[dds_name]]  <- list(
 Upregulated = dds_Genotype_res_list[[dds_name]] %>% dplyr::filter(padj <= 0.05) %>% dplyr::filter(log2FoldChange >= 1),
 Downregulated = dds_Genotype_res_list[[dds_name]] %>% dplyr::filter(padj <= 0.05) %>% dplyr::filter(log2FoldChange <= -1))
}

##save the control AAvsCA comparison DEGs
write.table(
rbind(dds_Genotype_res_sig_list$control$Upregulated,
      dds_Genotype_res_sig_list$control$Downregulated),
"/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/06_ANALYSIS/EpiPotato_mRNA_HighStress/AnnabelleVsCamel/AAvsCA_Control_DEGs.txt",
row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

AAvsCA_Control_DEGs <- read.delim("/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/06_ANALYSIS/EpiPotato_mRNA_HighStress/AnnabelleVsCamel/AAvsCA_Control_DEGs.txt")
AAvsCA_Control_DEGs %>% dplyr::filter(padj <= 0.05) %>% dplyr::filter(log2FoldChange >= 1) %>% nrow()
AAvsCA_Control_DEGs %>% dplyr::filter(padj <= 0.05) %>% dplyr::filter(log2FoldChange <= -1) %>% nrow()
##save the Heat AAvsCA comparison DEGs
write.table(
  rbind(dds_Genotype_res_sig_list$Heat$Upregulated,
        dds_Genotype_res_sig_list$Heat$Downregulated),
  "/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/06_ANALYSIS/EpiPotato_mRNA_HighStress/AnnabelleVsCamel/AAvsCA_Heat_DEGs.txt",
  row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

ggvenn::ggvenn(
  list(
    control = dds_Genotype_res_sig_list$control$Upregulated$locusName,
    Heat = dds_Genotype_res_sig_list$Heat$Upregulated$locusName
  ),
  text_size = 8
)

ggvenn::ggvenn(
  list(
    control = dds_Genotype_res_sig_list$control$Downregulated$locusName,
    Heat = dds_Genotype_res_sig_list$Heat$Downregulated$locusName
  ),
  text_size = 8
)

GO_AAvsCA<- lapply(dds_Genotype_res_sig_list, GO_KO_enrichment_with_enricher_CP,
       term2gene_GO = term2gene_GO_DB,
       term2name_GO = term2name_GO_DB)



TOP10_GO_AAvsCA_Control<- CombineUpDownRegulation_SelectTopEnrichmentTermsByCounts_Plot(
  TopCounts = 10,
  GO_AAvsCA$control$Upregulated@result,
  GO_AAvsCA$control$Downregulated@result
)



ggsave(
  "GO_AAvsCA_Control.svg",
  plot = ggplot(TOP10_GO_AAvsCA_Control,
                aes(x=GeneRatioNumeric,y=Description, fill = p.adjust)) + 
    geom_bar(stat = "identity") + 
    #scale_fill_manual(values=GO_labelColors)+
    scale_fill_continuous(low="red", high="blue")+
    geom_text(aes(label = Count), vjust = 0.5, hjust= 0.001,colour = "black", size =4.5, fontface=2) +
    theme_bw()  + 
    facet_grid(~Regulation) +
    xlab("Gene ratio") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size=12,face="bold"),
          axis.text.y= element_text(size=15,face="bold"),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size=15,face = "bold"),
          legend.title = element_text(size=15,face = "bold"),
          legend.text  = element_text(size=15,face = "bold"),
          strip.text = element_text(size=15,face = "bold")),
  device = "svg",
  path = "/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/06_ANALYSIS/EpiPotato_mRNA_HighStress/AnnabelleVsCamel",
  width = 14,
  height = 10,
  dpi = 300
)

###COMPARE TOTAL DEGS with TOTAL DMR DEGs

###Read the DMR DEGs
library(readxl)
### file path
file_path <- "/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/06_ANALYSIS/EpiPotato_mRNA_HighStress/AnnabelleVsCamel/DMR_DEGs/DMR_DEGsListByRegion_2ndExp_AAvsCA_Control.xlsx"

# Read all sheets into a list of data frames
DMR_DEGsListByRegion_2ndExp_AAvsCA_Control <- lapply(excel_sheets(file_path), read_excel, path = file_path)

DMR_DEGsListByRegion_2ndExp_AAvsCA_Control
names(DMR_DEGsListByRegion_2ndExp_AAvsCA_Control) <- excel_sheets(file_path)

### Consider all Genes that are methylated, as just differentially methylated or not within three regions... without separating context....
library(ggvenn)

##plot venn diagramm to visualize how methylated the DEGs are...
ggvenn(
list(
  Upstream = unique(DMR_DEGsListByRegion_2ndExp_AAvsCA_Control$Upstream$locusName),
  Downstream = unique(DMR_DEGsListByRegion_2ndExp_AAvsCA_Control$Downstream$locusName),
  GeneBody = unique(DMR_DEGsListByRegion_2ndExp_AAvsCA_Control$GeneBody$locusName)
))

####Now consider all DEGs that are methylated regardless to be counted only once... 
## if DEG is differentially methylated in genebody and upstream.. (count as one)
##first combine all genes together...

Total_DMRDEGs<- unique(c(DMR_DEGsListByRegion_2ndExp_AAvsCA_Control$Upstream$locusName,
                  DMR_DEGsListByRegion_2ndExp_AAvsCA_Control$Downstream$locusName,
                  DMR_DEGsListByRegion_2ndExp_AAvsCA_Control$GeneBody$locusName))


###Read total DEGs identified
AAvsCA_Control_DEGs <- read.delim("/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/06_ANALYSIS/EpiPotato_mRNA_HighStress/AnnabelleVsCamel/AAvsCA_Control_DEGs.txt")

##REMOVE THE SCAFFOLDS
AAvsCA_Control_DEGs<- AAvsCA_Control_DEGs %>% dplyr::filter(!str_detect(locusName,"Soltu.DM.S"))
##Determine how many DEGs are not differentially methylated
ggvenn(list(
 Total_DMRDEGs = Total_DMRDEGs,
 Total_DEGs = AAvsCA_Control_DEGs$locusName
))

##GET non DMR DEGs
Non_DMR_DEGs <- AAvsCA_Control_DEGs[!AAvsCA_Control_DEGs$locusName %in% Total_DMRDEGs,]

ggvenn::ggvenn(
  list(
    DMR_DEGs = Total_DMRDEGs,
    Non_DMR_DEGs = Non_DMR_DEGs$locusName
  ),text_size = 8
)


library(clusterProfiler)
###Up- and downregulated genes together

EnrichmentDMR_DEGs_Non_DMR_DEGs<- GO_KO_enrichment_with_enricher_CP(
  list(
    DMR_DEGs = Total_DMRDEGs,
    Non_DMR_DEGs = Non_DMR_DEGs$locusName
  ),
  term2gene_GO = term2gene_GO_DB,
  term2name_GO = term2name_GO_DB
)



ggplot(data = data.frame(
  Set=c("Total_DEGs","DMR_DEGs","Non_DMR_DEGs"),
  Size=c(length(AAvsCA_Control_DEGs$locusName), length(Total_DMRDEGs), length(Non_DMR_DEGs$locusName))), 
  aes(x = Set, y = Size, fill = Set)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = c("skyblue3", "green4","purple4")) +
  theme_minimal(base_size = 14) +
  labs(x = "", y = "DEG Counts") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size=12,face="bold"),
        axis.text.y= element_text(size=15,face="bold"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=15,face = "bold"),
        legend.position = "none",
        strip.text = element_text(size=15,face = "bold"))

library(eulerr)



CombineUpDownRegulation_SelectTopEnrichmentTermsByCounts_Plot <- function(TopCounts=20, 
                                                                          DMR_DEGs_table,
                                                                          Non_DMR_DEGs_table
                                                                          ){
  CombineRegulation_DF <- rbind(DMR_DEGs_table %>% 
                                  mutate(Regulation="DMR_DEGs") %>% 
                                  select(Description,GeneRatio,Count,p.adjust,Regulation,geneID) %>% 
                                  arrange(desc(Count), .by_group =TRUE) %>% dplyr::filter(p.adjust < 0.05) %>%
                                  slice(1:TopCounts),
                                Non_DMR_DEGs_table %>% 
                                  mutate(Regulation="Non_DMR_DEGs") %>% 
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

TOP20_GO_Enrichment<- CombineUpDownRegulation_SelectTopEnrichmentTermsByCounts_Plot(
  TopCounts = 20,
  DMR_DEGs_table = EnrichmentDMR_DEGs_Non_DMR_DEGs$DMR_DEGs@result,
  Non_DMR_DEGs_table = EnrichmentDMR_DEGs_Non_DMR_DEGs$Non_DMR_DEGs@result
)
TOP10_GO_Enrichment<- CombineUpDownRegulation_SelectTopEnrichmentTermsByCounts_Plot(
  TopCounts = 10,
  DMR_DEGs_table = EnrichmentDMR_DEGs_Non_DMR_DEGs$DMR_DEGs@result,
  Non_DMR_DEGs_table = EnrichmentDMR_DEGs_Non_DMR_DEGs$Non_DMR_DEGs@result
)

ggplot(TOP20_GO_Enrichment,
       aes(x=Count,y=Description, fill = p.adjust)) + 
  geom_bar(stat = "identity") + 
  #scale_fill_manual(values=GO_labelColors)+
  scale_fill_continuous(low="red", high="blue")+
  geom_text(aes(label = Count), vjust = 0.5, hjust= 0.001,colour = "black", size =4.5, fontface=2) +
  theme_bw()  + 
  facet_grid(~Regulation) +
  xlab("Gene ratio") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size=12,face="bold"),
        axis.text.y= element_text(size=15,face="bold"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=15,face = "bold"),
        legend.title = element_text(size=15,face = "bold"),
        legend.text  = element_text(size=15,face = "bold"),
        strip.text = element_text(size=15,face = "bold"))


ggplot(TOP10_GO_Enrichment,
       aes(x=Count,y=Description, fill = p.adjust)) + 
  geom_bar(stat = "identity") + 
  #scale_fill_manual(values=GO_labelColors)+
  scale_fill_continuous(low="red", high="blue")+
  geom_text(aes(label = Count), vjust = 0.5, hjust= 0.001,colour = "black", size =4.5, fontface=2) +
  theme_bw()  + 
  facet_grid(~Regulation) +
  xlab("Gene count") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size=12,face="bold"),
        axis.text.y= element_text(size=15,face="bold"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=15,face = "bold"),
        legend.title = element_text(size=15,face = "bold"),
        legend.text  = element_text(size=15,face = "bold"),
        strip.text = element_text(size=15,face = "bold"))


###Up- and downregulated genes separate
Total_DMRDEGs_DF<- data.frame(locusName = Total_DMRDEGs)
Total_DMRDEGs_DF_L2FC <- inner_join(Total_DMRDEGs_DF,
                                    AAvsCA_Control_DEGs,by="locusName")


DMR_DEGs_and_non_DMR_DEGs_list<- list(
            DMR_DEGS_upregulated = Total_DMRDEGs_DF_L2FC %>% dplyr::filter(log2FoldChange >= 1),
            DMR_DEGS_downregulated = Total_DMRDEGs_DF_L2FC %>% dplyr::filter(log2FoldChange <= -1),
            
            Non_DMR_DEGs_upregulated = Non_DMR_DEGs %>% dplyr::filter(log2FoldChange >= 1),
            Non_DMR_DEGs_downregulated = Non_DMR_DEGs %>% dplyr::filter(log2FoldChange <= -1) )

DMR_DEGs_and_non_DMR_DEGs_Vector <- list()
for (comparison in names(DMR_DEGs_and_non_DMR_DEGs_list)) {
  
  DMR_DEGs_and_non_DMR_DEGs_Vector[[comparison]]<- as.character(DMR_DEGs_and_non_DMR_DEGs_list[[comparison]]$locusName)
  
}


ggvenn(DMR_DEGs_and_non_DMR_DEGs_Vector,
       text_size = 8)

EnrichmentDMR_DEGs_Non_DMR_DEGs_UpDown<- GO_KO_enrichment_with_enricher_CP(
  DMR_DEGs_and_non_DMR_DEGs_Vector,
  term2gene_GO = term2gene_GO_DB,
  term2name_GO = term2name_GO_DB
)


CombineUpDownRegulation_SelectTopEnrichmentTermsByCounts_Plot <- function(TopCounts=20, 
                                                                          DMR_DEGs_table_Up,
                                                                          DMR_DEGs_table_Down,
                                                                          Non_DMR_DEGs_table_Up,
                                                                          Non_DMR_DEGs_table_Dowm
){
  CombineRegulation_DF <- rbind(DMR_DEGs_table_Up %>% 
                                  mutate(Regulation="DMR_DEGs_Upreugulated") %>% 
                                  select(Description,GeneRatio,Count,p.adjust,Regulation,geneID) %>% 
                                  arrange(desc(Count), .by_group =TRUE) %>% dplyr::filter(p.adjust < 0.05) %>%
                                  slice(1:TopCounts),
                                DMR_DEGs_table_Down %>% 
                                  mutate(Regulation="DMR_DEGs_Downregulated") %>% 
                                  select(Description,GeneRatio,Count,p.adjust,Regulation,geneID) %>% 
                                  arrange(desc(Count), .by_group =TRUE) %>% dplyr::filter(p.adjust < 0.05) %>%
                                  slice(1:TopCounts),
                                Non_DMR_DEGs_table_Up %>% 
                                  mutate(Regulation="Non_DMR_DEGs_Upregulated") %>% 
                                  select(Description,GeneRatio,Count,p.adjust,Regulation,geneID) %>% 
                                  arrange(desc(Count), .by_group =TRUE) %>% dplyr::filter(p.adjust < 0.05) %>%
                                  slice(1:TopCounts),
                                Non_DMR_DEGs_table_Dowm %>% 
                                  mutate(Regulation="Non_DMR_DEGs_Downregulated") %>% 
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

TOP10_GO_Enrichment_UpDown <- CombineUpDownRegulation_SelectTopEnrichmentTermsByCounts_Plot(
  TopCounts=10,
  DMR_DEGs_table_Up = EnrichmentDMR_DEGs_Non_DMR_DEGs_UpDown$DMR_DEGS_upregulated@result,
  DMR_DEGs_table_Down = EnrichmentDMR_DEGs_Non_DMR_DEGs_UpDown$DMR_DEGS_downregulated@result,
  Non_DMR_DEGs_table_Up = EnrichmentDMR_DEGs_Non_DMR_DEGs_UpDown$Non_DMR_DEGs_upregulated@result,
  Non_DMR_DEGs_table_Dowm = EnrichmentDMR_DEGs_Non_DMR_DEGs_UpDown$Non_DMR_DEGs_downregulated@result)

ggplot(TOP10_GO_Enrichment_UpDown,
       aes(x=Count,y=Description, fill = p.adjust)) + 
  geom_bar(stat = "identity") + 
  #scale_fill_manual(values=GO_labelColors)+
  scale_fill_continuous(low="red", high="blue")+
  geom_text(aes(label = Count), vjust = 0.5, hjust= 0.001,colour = "black", size =4.5, fontface=2) +
  theme_bw()  + 
  facet_grid(~Regulation) +
  xlab("Gene count") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size=12,face="bold"),
        axis.text.y= element_text(size=15,face="bold"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=15,face = "bold"),
        legend.title = element_text(size=15,face = "bold"),
        legend.text  = element_text(size=15,face = "bold"),
        strip.text = element_text(size=15,face = "bold"))


CompareClusterEnrichmentDMR_DEGs_Non_DMR_DEGs_UpDown<- compareCluster(
  DMR_DEGs_and_non_DMR_DEGs_Vector, 
  fun = "enricher",
  pvalueCutoff = 0.05,
  pAdjustMethod = "fdr",
  TERM2GENE = term2gene_GO_DB,
  TERM2NAME = term2name_GO_DB)

col_fun <- colorRamp2(c(-2, 0, 2), c("#2A78C6","white","#C83B38"))

dotplot(CompareClusterEnrichmentDMR_DEGs_Non_DMR_DEGs_UpDown,
       showCategory = 10,
       by = "Count",
       size = "Count") +
  scale_fill_gradient2(
    low = "#C83B38", mid = "white", high = "#2A78C6",
    midpoint = 0.025) +
  theme_minimal() + 
  theme(axis.text.x = element_text(hjust = 1, size = 15, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        legend.text = element_text(size =20, face ="bold"),
        legend.title = element_text(size =20, face ="bold"))

##Get the genes of ORA analysis

EnrichmentDMR_DEGs_Non_DMR_DEGs_UpDown
###Get the functional annotation of the GO Terms
EnrichmentDMR_DEGs_Non_DMR_DEGs_UpDown_functannot <- list()
for (Intersections in names(EnrichmentDMR_DEGs_Non_DMR_DEGs_UpDown)) {
  
  #Separate table of both GO and KO Enrichment results and store to list
  GO_table <- EnrichmentDMR_DEGs_Non_DMR_DEGs_UpDown[[Intersections]]@result
  
  #rename the core_enrichment genes..
  colnames(GO_table)[colnames(GO_table) == "geneID"] <- "locusName"
  
  
  #Create the enrichment term table (summary)
  GO_table_summary<- GO_table %>% 
    dplyr::select(locusName,ID,Description,pvalue,p.adjust) %>% dplyr::filter(p.adjust < 0.05) %>%
    separate_rows(locusName, sep = "/") %>% left_join(S.tuberosum.v6.1_latest %>% 
                                                        dplyr::select(locusName,v6.1.Description,aktualisierte.annotation..bitte.ergänzen.,
                                                                      Mapman.description), by ="locusName")
  
  
  ##save everything to lsit
  EnrichmentDMR_DEGs_Non_DMR_DEGs_UpDown_functannot[[Intersections]] <- GO_table_summary
  
}
#####Ensure that each gene is unique...... and combine all GO terms and respective enrichment scores,pvalues etc.. of each gene together
EnrichmentDMR_DEGs_Non_DMR_DEGs_UpDown_functannot_unique <- list()
for (Intersections in names(EnrichmentDMR_DEGs_Non_DMR_DEGs_UpDown_functannot)) {
  
  #get the comparison name with both Pos and Neg table
  ORA_comparison <- EnrichmentDMR_DEGs_Non_DMR_DEGs_UpDown_functannot[[Intersections]]
  
  ORA_comparison_unique <- ORA_comparison %>%
    group_by(locusName, v6.1.Description,aktualisierte.annotation..bitte.ergänzen.,Mapman.description) %>%
    dplyr::rename(GO_description = Description) %>%
    summarize(
      ID = paste(ID, collapse = "/"),
      GO_description = paste(GO_description, collapse = "/"),
      pvalue = paste(pvalue, collapse = "/"),
      p.adjust = paste(p.adjust, collapse = "/")) %>% 
    dplyr::select(locusName,ID,GO_description,pvalue,p.adjust,v6.1.Description,aktualisierte.annotation..bitte.ergänzen.,Mapman.description)
  
  #save to the master list
  EnrichmentDMR_DEGs_Non_DMR_DEGs_UpDown_functannot_unique[[Intersections]] <- ORA_comparison_unique
  
}



View(EnrichmentDMR_DEGs_Non_DMR_DEGs_UpDown_functannot_unique$DMR_DEGS_upregulated)

#save all enrihcment DFs

library(openxlsx)

write.xlsx(EnrichmentDMR_DEGs_Non_DMR_DEGs_UpDown_functannot_unique, 
           file = "/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/06_ANALYSIS/EpiPotato_mRNA_HighStress/AnnabelleVsCamel/DMR_DEGs/EnrichmentDMR_DEGs.xlsx")



##### Determine significance of TOTAL genes that are DMRs and DEGs.....

#32819 that were annotated (Has genomic coordinates and is within chromosomes)

###Total number of genes in v6.1 
nrow(UniqueCountsDuplicatedHighStress)

library(openxlsx)

# Load workbook
wb <- loadWorkbook("/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/06_ANALYSIS/EpiPotato_mRNA_HighStress/AnnabelleVsCamel/AAvsCA_C_2nd.xlsx")

DMR_Genes_Total<- list(
  Upstream = read.xlsx(wb, sheet = 1),
  GeneBody = read.xlsx(wb, sheet = 2),
  Downstream = read.xlsx(wb, sheet = 3))

Total_DMR_Genes <- unique(c(DMR_Genes_Total$Upstream$locusName,
                         DMR_Genes_Total$Downstream$locusName,
                         DMR_Genes_Total$GeneBody$locusName))


###number of DMR GENES that are not DEGs (NOTE::DIFFERNTIALLY METHYLATED but differentially expressed or not)
Total_DMR_Genes
Total_DMRDEGs

ggvenn(list(
  DMR_Genes= Total_DMR_Genes,
 DMR_DEGs = Total_DMRDEGs
),text_size = 8)

####FIND out total non_DMRs GENES and non DEGs..

##FIRST DETERMIN NUMBER OF NON_DMR GENEs
#14450 Difference between total genes and DMR genes
TOTALGENE_ID<- UniqueCountsDuplicatedHighStress %>% 
  dplyr::rename(locusName=Geneid) %>% dplyr::filter(str_detect(Chr,"chr")) %>%
  dplyr::select(locusName) 

TOTALGENE_ID$locusName <- gsub(".v6.1","",TOTALGENE_ID$locusName)

### number of non DMR genes.. exclude scaffold genes
NON_DMR_GENES <- TOTALGENE_ID[!TOTALGENE_ID$locusName %in% Total_DMR_Genes,]

###Now determine of these non DMR genes, how many are differentially expressed.. mapped to total degs
nrow(AAvsCA_Control_DEGs)

NON_DMR_DEGs <-  NON_DMR_GENES[NON_DMR_GENES %in% AAvsCA_Control_DEGs$locusName]
NON_DMR_NON_DEGs <-  NON_DMR_GENES[!NON_DMR_GENES %in% AAvsCA_Control_DEGs$locusName]


##contingency table of total genes, dmr genes, dmr degs etc
ContingencyTableDMR_DEGs <- read.delim("/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/06_ANALYSIS/EpiPotato_mRNA_HighStress/AnnabelleVsCamel/ContingencyTableDMR_DEGs", comment.char="#")
ContingencyTableDMR_DEGs <-  ContingencyTableDMR_DEGs%>% column_to_rownames("X")

fisherTestDMR_DEGs<- fisher.test(as.table(as.matrix(ContingencyTableDMR_DEGs)))

# create dataframe from contingency table
x <- c()
for (row in rownames(ContingencyTableDMR_DEGs)) {
  for (col in colnames(ContingencyTableDMR_DEGs)) {
    x <- rbind(x, matrix(rep(c(row, col), ContingencyTableDMR_DEGs[row, col]), ncol = 2, byrow = TRUE))
  }
}
df <- as.data.frame(x)
colnames(df) <- c("ExpressionDifference", "MethylationDifference")
df

library(ggstatsplot)
ggbarstats(
  df, MethylationDifference, ExpressionDifference,
  results.subtitle = FALSE,
  subtitle = paste0(
    "Fisher's exact test", ", p-value = ",
    ifelse(fisherTestDMR_DEGs$p.value < 0.001, "< 0.001", round(test$p.value, 3))
  )
) +   theme(axis.text.x = element_text(angle = 60, hjust = 1, size=12,face="bold"),
            axis.text.y= element_text(size=15,face="bold"),
            axis.title.y = element_blank(),
            axis.title.x = element_text(size=15,face = "bold"),
            legend.title = element_text(size=15,face = "bold"),
            legend.text  = element_text(size=15,face = "bold"),
            strip.text = element_text(size=15,face = "bold"),
            title = element_text(size=15,face = "bold"))


flipped_table <- matrix(c(1200, 13250, 1764, 16605), nrow = 2, byrow = TRUE)


colnames(flipped_table) <- c("NON_DMR", "DMR")
rownames(flipped_table) <- c("DEGs", "NON_DEGs")

fisher.test(flipped_table)

library(gplots)
balloonplot(as.table(as.matrix(ContingencyTableDMR_DEGs)))

mosaicplot(ContingencyTableDMR_DEGs,
           main = "Mosaic plot",
           color = T
)


### pvalue was  <0.05, odds are the DEGs are 1.17 times more likely to be DMRs than non-DEGs
chisq.test(as.table(as.matrix(ContingencyTableDMR_DEGs)))$expected
fisher.test(as.table(as.matrix(ContingencyTableDMR_DEGs)))



###Determine if there is bias in Hyper and hypomethylation regarding up and downregulation
##Consider all DMR DEGs...
###First switch Annabelle as the comparison....
DMR_DEGsListByRegion_2ndExp_AAvsCA_Control_edit <- list()
for (DMR_type in names(DMR_DEGsListByRegion_2ndExp_AAvsCA_Control)) {
  
  ####MAKE ANNABELLE AS THE COMPARISON#####
  DMR_DEGsListByRegion_2ndExp_AAvsCA_Control_edit[[DMR_type]] <-  DMR_DEGsListByRegion_2ndExp_AAvsCA_Control[[DMR_type]] %>% 
     dplyr::mutate(Methylation_Difference = case_when(str_detect(Methylation_Difference,"gain")~"hypomethylation",
                                                      str_detect(Methylation_Difference,"loss")~"hypermethylation"),
                   
                   Expression_Direction = case_when(Expression_L2FC > 0 ~ "Positive",
                                                    Expression_L2FC < 0 ~ "Negative"))
  
}

View(DMR_DEGsListByRegion_2ndExp_AAvsCA_Control_edit$Upstream)

##ANNABELLE IS THE COMPARISON
##HYPERMETHYLATION means annablle is hypermethylation
##positive L2FC means ANnabelle is higher expressed
TotalHypoUpregulated<- rbind(DMR_DEGsListByRegion_2ndExp_AAvsCA_Control_edit$Upstream %>% dplyr::filter(Methylation_Difference =="hypomethylation" | Expression_L2FC > 0) ,
                             DMR_DEGsListByRegion_2ndExp_AAvsCA_Control_edit$Downstream %>% dplyr::filter(Methylation_Difference =="hypomethylation" | Expression_L2FC > 0) ,
                             DMR_DEGsListByRegion_2ndExp_AAvsCA_Control_edit$GeneBody %>% dplyr::filter(Methylation_Difference =="hypomethylation" | Expression_L2FC > 0))

CombineALLDMRTYPesNdetermineMethylationDirectionNexpression<- function(DMR_DEG_LIST, SelectMEthylationDIrection, ExpressionDirection){
  
  ## Combine all regions
  MethylationRegulationDirection <- rbind(
    DMR_DEG_LIST$Upstream %>%
      dplyr::filter(Methylation_Difference == SelectMEthylationDIrection & Expression_Direction == ExpressionDirection),
    
    DMR_DEG_LIST$Downstream %>%
      dplyr::filter(Methylation_Difference == SelectMEthylationDIrection & Expression_Direction == ExpressionDirection),
    
    DMR_DEG_LIST$GeneBody %>%
      dplyr::filter(Methylation_Difference == SelectMEthylationDIrection & Expression_Direction == ExpressionDirection)
  )
  
  
  ###return only neccessary columns
  MethylationRegulationDirection_final <- MethylationRegulationDirection %>% dplyr::select(locusName, Methylation_Difference, Expression_Direction)
  
  
  return(MethylationRegulationDirection_final)
  
}

Total_DMR_DEG_Regulation_Types<- list(
TotalHypoUpregulated = CombineALLDMRTYPesNdetermineMethylationDirectionNexpression(DMR_DEGsListByRegion_2ndExp_AAvsCA_Control_edit,SelectMEthylationDIrection = "hypomethylation", ExpressionDirection = "Positive"),
TotalHypoDownregulated = CombineALLDMRTYPesNdetermineMethylationDirectionNexpression(DMR_DEGsListByRegion_2ndExp_AAvsCA_Control_edit,SelectMEthylationDIrection = "hypomethylation", ExpressionDirection = "Negative"),
TotalHyperUpregulated = CombineALLDMRTYPesNdetermineMethylationDirectionNexpression(DMR_DEGsListByRegion_2ndExp_AAvsCA_Control_edit,SelectMEthylationDIrection = "hypermethylation", ExpressionDirection = "Positive"),
TotalHyperDownregulated = CombineALLDMRTYPesNdetermineMethylationDirectionNexpression(DMR_DEGsListByRegion_2ndExp_AAvsCA_Control_edit,SelectMEthylationDIrection = "hypermethylation", ExpressionDirection = "Negative")
)

##Count total.. genes
Total_DMR_DEG_Regulation_Types_unique <- list()
for (DMR_DEG_TYPE in names(Total_DMR_DEG_Regulation_Types)) {
  
  ##Ensure the gene row is unique
  Total_DMR_DEG_Regulation_Types_unique[[DMR_DEG_TYPE]] <- unique(Total_DMR_DEG_Regulation_Types[[DMR_DEG_TYPE]])
  
  ##count genes before unique function
  print(paste0("TOTAL BEFORE unique for ", DMR_DEG_TYPE,": ", nrow(Total_DMR_DEG_Regulation_Types[[DMR_DEG_TYPE]])))
  
  #count genes after unique function
  print(paste0("TOTAL AFTER unique for ", DMR_DEG_TYPE,": ", nrow(unique(Total_DMR_DEG_Regulation_Types[[DMR_DEG_TYPE]]))))
  
  #count only the locusName vector after the unique function
  print(paste0("TOTAL AFTER unique for just gene name ", DMR_DEG_TYPE,": ", length(unique(Total_DMR_DEG_Regulation_Types[[DMR_DEG_TYPE]]$locusName))))
  

}

nrow(TotalHyperUpregulated[TotalHyperUpregulated$locusName %in% TotalHypoUpregulated$locusName,])

##Look for overlaps
#04G004920 is an example of hyper and hypo methylated gene but positively upregulated 
ggvenn(list(
  HypoUpregulated= Total_DMR_DEG_Regulation_Types$TotalHypoUpregulated$locusName,
  HypoDownregulated = Total_DMR_DEG_Regulation_Types$TotalHypoDownregulated$locusName,
  HyperUpregulated= Total_DMR_DEG_Regulation_Types$TotalHyperUpregulated$locusName,
  HyperDownregulated= Total_DMR_DEG_Regulation_Types$TotalHyperDownregulated$locusName
),
       text_size = 10)


# Convert the list to lupset df
library(UpSetR)
geneSets <- list(
  HypoUpregulated= Total_DMR_DEG_Regulation_Types$TotalHypoUpregulated$locusName,
  HypoDownregulated = Total_DMR_DEG_Regulation_Types$TotalHypoDownregulated$locusName,
  HyperUpregulated= Total_DMR_DEG_Regulation_Types$TotalHyperUpregulated$locusName,
  HyperDownregulated= Total_DMR_DEG_Regulation_Types$TotalHyperDownregulated$locusName
)
GeneID <- sort(unique(unlist(geneSets)))

mat <- t(+sapply(geneSets, "%in%", x = GeneID))  ## matrix output
colnames(mat) <- GeneID

UpsetDF_DMRDEG_HYPERHYPO<- t(data.frame(mat))
UpsetDF_DMRDEG_HYPERHYPO <- UpsetDF_DMRDEG_HYPERHYPO %>% as.data.frame() %>% rownames_to_column("Gene")

library(micro.gen.extra)
HyperHypo_Downregulated<- data.frame(locusName = get_intersect_members(UpsetDF_DMRDEG_HYPERHYPO, c("HyperDownregulated","HypoDownregulated"), exclusive = TRUE))
HyperHypo_Upregulated<- data.frame(locusName = get_intersect_members(UpsetDF_DMRDEG_HYPERHYPO, c("HyperUpregulated","HypoUpregulated"), exclusive = TRUE))

HyperHypo_Downregulated_annot <- HyperHypo_Downregulated %>% inner_join(S.tuberosum.v6.1_latest)
HyperHypo_Upregulated_annot <- HyperHypo_Upregulated %>% inner_join(S.tuberosum.v6.1_latest)

write.table(HyperHypo_Downregulated_annot,"/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/06_ANALYSIS/EpiPotato_mRNA_HighStress/AnnabelleVsCamel/DMR_DEGs/HyperHypo_Downregulated_annot.txt",col.names = TRUE,row.names = FALSE)
write.table(HyperHypo_Upregulated_annot,"/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/06_ANALYSIS/EpiPotato_mRNA_HighStress/AnnabelleVsCamel/DMR_DEGs/HyperHypo_Upregulated_annot.txt",col.names = TRUE,row.names = FALSE)

UpSetR::upset(as.data.frame(UpsetDF_DMRDEG_HYPERHYPO),
              text.scale = c(2, 2, 2, 2, 2, 2),
              sets = c("HypoUpregulated",
                       "HypoDownregulated",
                       "HyperUpregulated",
                       "HyperDownregulated"),
              keep.order = TRUE)


#contingency table
ContingencyTableDMR_DEGSHyperHypoUnique <- ContingencyTableDMR_DEGSHyperHypoUnique %>% column_to_rownames("X")

ContingencyTableDMR_DEGSHyperHypo_NON_unique<- ContingencyTableDMR_DEGSHyperHypo_NON_unique %>% column_to_rownames("X")

chisq.test(ContingencyTableDMR_DEGSHyperHypoUnique)
fisher.test(ContingencyTableDMR_DEGSHyperHypoUnique)

chisq.test(ContingencyTableDMR_DEGSHyperHypo_NON_unique)
fisher.test(ContingencyTableDMR_DEGSHyperHypo_NON_unique)

nrow(TotalHypoUpregulated)
nrow(unique(TotalHypoUpregulated))
length(unique(TotalHypoUpregulated$locusName))

###upload the v3.4 annotation
Soltu_DM_und_PGSC_Nummern <- read_excel("/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/00_META_DATA/Soltu.DM und PGSC Nummern.xlsx")

##Get the soltu to v3.4 DMP annotation
Soltu_DM_v3.4_annotaion <- Soltu_DM_und_PGSC_Nummern %>% dplyr::select(`Feature ID`,`v3.4 DMP`) %>%
  dplyr::rename(locusName=`Feature ID`)

Total_DMRDEGs_L2FC<- data.frame(locusName=Total_DMRDEGs)

##Get the log2foldchange...
Total_DEG_types_L2FC<- list(
 Total_DMRDEGs_L2FC = inner_join(Total_DMRDEGs_L2FC, AAvsCA_Control_DEGs %>% dplyr::select(locusName,log2FoldChange, v6.1.Description,Mapman.description), by="locusName"),
 Non_DMR_DEGs_L2FC = Non_DMR_DEGs %>% dplyr::select(locusName, log2FoldChange, v6.1.Description, Mapman.description) )


##Annotate to v3.4 genome
Total_DEG_types_L2FC_v3.4 <- list()
for (DEG_types in names(Total_DEG_types_L2FC)) {

  ##annotate for both DF to PGSC annotation  
  Total_DEG_types_L2FC_v3.4[[DEG_types]] <- Total_DEG_types_L2FC[[DEG_types]] %>% 
                                    inner_join(Soltu_DM_v3.4_annotaion, by="locusName") %>% 
                                    dplyr::select(locusName,`v3.4 DMP`, log2FoldChange, v6.1.Description, Mapman.description)
  
}

##Check for missing v3.4 annotation
Total_DEG_types_L2FC_v3.4$Total_DMRDEGs_L2FC[!complete.cases(Total_DEG_types_L2FC_v3.4$Total_DMRDEGs_L2FC),]
Total_DEG_types_L2FC_v3.4$Non_DMR_DEGs_L2FC[!complete.cases(Total_DEG_types_L2FC_v3.4$Non_DMR_DEGs_L2FC),]

##omit genes with missing v3.4, all except one was hypothetical
Total_DMRDEGs_L2FC_v3.4<- Total_DEG_types_L2FC_v3.4$Total_DMRDEGs_L2FC %>% dplyr::select(`v3.4 DMP`, log2FoldChange) %>% na.omit()
Non_DMR_DEGs_L2FC_v3.4<- Total_DEG_types_L2FC_v3.4$Non_DMR_DEGs_L2FC %>% dplyr::select(`v3.4 DMP`, log2FoldChange) %>% na.omit()

Total_DMRDEGs_L2FC_v3.4[duplicated(Total_DMRDEGs_L2FC_v3.4),]

duplicated(Total_DMRDEGs_L2FC_v3.4$`v3.4 DMP`)
Total_DMRDEGs_L2FC_v3.4[duplicated(Total_DMRDEGs_L2FC_v3.4[[1]]), ]


Total_DMRDEGs_L2FC_v3.4_duplicates <- Total_DMRDEGs_L2FC_v3.4 %>%
  group_by(`v3.4 DMP`) %>%
  filter(n() > 1) %>%
  ungroup() %>% dplyr::arrange(`v3.4 DMP`)

##Get the v3.4 ids that conflict each other for DMR_DEGs
conflicting_ids_DMR_DEGs <- Total_DMRDEGs_L2FC_v3.4 %>%
  group_by(`v3.4 DMP`) %>%
  filter(n() > 1) %>%
  summarise(
    has_positive = any(log2FoldChange > 0),
    has_negative = any(log2FoldChange < 0),
    .groups = 'drop'
  ) %>%
  filter(has_positive & has_negative) %>%
  pull(`v3.4 DMP`)  # Extract just the IDs

# Step 2: Get all rows with those conflicting IDs
conflicting_duplicate_rows_DMR_DEGs <- Total_DMRDEGs_L2FC_v3.4 %>%
  filter(`v3.4 DMP` %in% conflicting_ids_DMR_DEGs)

##first remove duplicated rows with conflicting l2fc
Total_DMRDEGs_L2FC_v3.4_final<- Total_DMRDEGs_L2FC_v3.4[!Total_DMRDEGs_L2FC_v3.4$`v3.4 DMP` %in% conflicting_duplicate_rows_DMR_DEGs$`v3.4 DMP`,]

##Second, remove duplicated rows with log2fold change of similar magnitude in same direction
Total_DMRDEGs_L2FC_v3.4_final<- Total_DMRDEGs_L2FC_v3.4_final %>% 
#dplyr::arrange(log2FoldChange) %>% ###ignore for the moment 
  distinct(`v3.4 DMP`, .keep_all = TRUE)




##Get the v3.4 ids that conflict each other for DMR_DEGs
conflicting_ids_Non_DMR_DEGs <- Non_DMR_DEGs_L2FC_v3.4 %>%
  group_by(`v3.4 DMP`) %>%
  filter(n() > 1) %>%
  summarise(
    has_positive = any(log2FoldChange > 0),
    has_negative = any(log2FoldChange < 0),
    .groups = 'drop'
  ) %>%
  filter(has_positive & has_negative) %>%
  pull(`v3.4 DMP`)  # Extract just the IDs

# Step 2: Get all rows with those conflicting IDs
conflicting_duplicate_rows_Non_DMR_DEGs <- Non_DMR_DEGs_L2FC_v3.4 %>%
  filter(`v3.4 DMP` %in% conflicting_ids_DMR_DEGs)

##first remove duplicated rows with conflicting l2fc
Non_DMRDEGs_L2FC_v3.4_final<- Non_DMR_DEGs_L2FC_v3.4[!Non_DMR_DEGs_L2FC_v3.4$`v3.4 DMP` %in% conflicting_duplicate_rows_Non_DMR_DEGs$`v3.4 DMP`,]

##Second, remove duplicated rows with log2fold change of similar magnitude in same direction
Non_DMRDEGs_L2FC_v3.4_final<- Non_DMRDEGs_L2FC_v3.4_final %>% 
  #dplyr::arrange(log2FoldChange) %>% ###ignore for the moment 
  distinct(`v3.4 DMP`, .keep_all = TRUE)




write.table(Total_DMRDEGs_L2FC_v3.4_final,"/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/06_ANALYSIS/EpiPotato_mRNA_HighStress/AnnabelleVsCamel/DMR_DEGs/Total_DMRDEGs_L2FC_v3.4_final.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(Non_DMRDEGs_L2FC_v3.4_final,"/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/06_ANALYSIS/EpiPotato_mRNA_HighStress/AnnabelleVsCamel/DMR_DEGs/Non_DMRDEGs_L2FC_v3.4_final.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")



write.table(Total_DMRDEGs_L2FC_v3.4, "/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/06_ANALYSIS/EpiPotato_mRNA_HighStress/AnnabelleVsCamel/DMR_DEGs/Total_DMRDEGs_L2FC_v3.4.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(Non_DMR_DEGs_L2FC_v3.4, "/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/06_ANALYSIS/EpiPotato_mRNA_HighStress/AnnabelleVsCamel/DMR_DEGs/Non_DMR_DEGs_L2FC_v3.4.txt",col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")






####FIMO identification of promoter boxes for upregulated genes and downregulated genes####
## 1332 upregulated...
## 1640 downregulated...
##load fimo
fimoUP <- read.delim("/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/06_ANALYSIS/EpiPotato_mRNA_HighStress/AnnabelleVsCamel/FIMO/up/fimoUP.tsv")
fimoDOWN <- read.delim("/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/06_ANALYSIS/EpiPotato_mRNA_HighStress/AnnabelleVsCamel/FIMO/down/fimoDOWN.tsv")

unique(fimoUP$motif_alt_id)
unique(fimoDOWN$motif_alt_id)

fimoUP_edit<- fimoUP %>% dplyr::mutate(MotifGroup=case_when(str_detect(motif_alt_id,"BPC") ~"BPC",
                                                            str_detect(motif_alt_id,"DREB") ~"DREB",
                                                            str_detect(motif_alt_id,"ERF") ~"ERF",
                                                            str_detect(motif_alt_id,"RAMOSA") ~"RAMOSA",
                                                            str_detect(motif_alt_id,"BZIP|bZIP") ~"BZIP",
                                                            str_detect(motif_alt_id,"DOF") ~"DOF",
                                                            str_detect(motif_alt_id,"CDF") ~"CDF",
                                                            str_detect(motif_alt_id,"LRL") ~"LRL",
                                                            str_detect(motif_alt_id,"AGL") ~"AGL",
                                                            str_detect(motif_alt_id,"NAC") ~"NAC",
                                                            str_detect(motif_alt_id,"MYB") ~"MYB",
                                                            str_detect(motif_alt_id,"SVP") ~"SVP",
                                                            str_detect(motif_alt_id,"bHLH") ~"bHLH",
                                                            str_detect(motif_alt_id,"RAP") ~"RAP",
                                                            str_detect(motif_alt_id,"PI") ~"PI",
                                                            str_detect(motif_alt_id,"O2") ~"O2",
                                                            str_detect(motif_alt_id,"DPBF") ~"DPBF",
                                                            str_detect(motif_alt_id,"BAD") ~"BAD",
                                                            str_detect(motif_alt_id,"abi4") ~"abi4",
                                                            str_detect(motif_alt_id,"RAX") ~"RAX",
                                                            str_detect(motif_alt_id,"CAMTA5") ~"CAMTA5",
                                                            str_detect(motif_alt_id,"HY5") ~"HY5",
                                                            str_detect(motif_alt_id,"AP1") ~"AP1",
                                                            str_detect(motif_alt_id,"FLC") ~"FLC",
                                                            str_detect(motif_alt_id,"TCP") ~"TCP",
                                                            str_detect(motif_alt_id,"EREB") ~"EREB",
                                                            str_detect(motif_alt_id,"TRP") ~"TRP",
                                                            str_detect(motif_alt_id,"LEP") ~"LEP",
                                                            str_detect(motif_alt_id,"LOB") ~"LOB",
                                                            str_detect(motif_alt_id,"AS2") ~"AS2",
                                                            str_detect(motif_alt_id,"ABF1") ~"ABF1",
                                                            str_detect(motif_alt_id,"LOB") ~"LOB",
                                                            str_detect(motif_alt_id,"SOC1") ~"SOC1",
                                                            str_detect(motif_alt_id,"AT4G12670") ~"Homeodomain-like superfamily protein",
                                                            str_detect(motif_alt_id,"HSF") ~"HSF",
                                                            str_detect(motif_alt_id,"CRF") ~"CRF")) %>% 
                                                            na.omit() %>%
                                                            dplyr::select(sequence_name,MotifGroup) %>% unique()

fimoDOWN_edit<- fimoDOWN %>% dplyr::mutate(MotifGroup=case_when(str_detect(motif_alt_id,"BPC") ~"BPC",
                                                            str_detect(motif_alt_id,"DREB") ~"DREB",
                                                            str_detect(motif_alt_id,"ERF") ~"ERF",
                                                            str_detect(motif_alt_id,"RAMOSA") ~"RAMOSA",
                                                            str_detect(motif_alt_id,"BZIP|bZIP") ~"BZIP",
                                                            str_detect(motif_alt_id,"DOF") ~"DOF",
                                                            str_detect(motif_alt_id,"CDF") ~"CDF",
                                                            str_detect(motif_alt_id,"LRL") ~"LRL",
                                                            str_detect(motif_alt_id,"AGL") ~"AGL",
                                                            str_detect(motif_alt_id,"NAC") ~"NAC",
                                                            str_detect(motif_alt_id,"MYB") ~"MYB",
                                                            str_detect(motif_alt_id,"SVP") ~"SVP",
                                                            str_detect(motif_alt_id,"bHLH") ~"bHLH",
                                                            str_detect(motif_alt_id,"RAP") ~"RAP",
                                                            str_detect(motif_alt_id,"PI") ~"PI",
                                                            str_detect(motif_alt_id,"O2") ~"O2",
                                                            str_detect(motif_alt_id,"DPBF") ~"DPBF",
                                                            str_detect(motif_alt_id,"BAD") ~"BAD",
                                                            str_detect(motif_alt_id,"abi4") ~"abi4",
                                                            str_detect(motif_alt_id,"RAX") ~"RAX",
                                                            str_detect(motif_alt_id,"CAMTA5") ~"CAMTA5",
                                                            str_detect(motif_alt_id,"HY5") ~"HY5",
                                                            str_detect(motif_alt_id,"AP1") ~"AP1",
                                                            str_detect(motif_alt_id,"FLC") ~"FLC",
                                                            str_detect(motif_alt_id,"TCP") ~"TCP",
                                                            str_detect(motif_alt_id,"EREB") ~"EREB",
                                                            str_detect(motif_alt_id,"TRP") ~"TRP",
                                                            str_detect(motif_alt_id,"LEP") ~"LEP",
                                                            str_detect(motif_alt_id,"LOB") ~"LOB",
                                                            str_detect(motif_alt_id,"AS2") ~"AS2",
                                                            str_detect(motif_alt_id,"ABF1") ~"ABF1",
                                                            str_detect(motif_alt_id,"LOB") ~"LOB",
                                                            str_detect(motif_alt_id,"SOC1") ~"SOC1",
                                                            str_detect(motif_alt_id,"AT4G12670") ~"Homeodomain-like superfamily protein",
                                                            str_detect(motif_alt_id,"HSF") ~"HSF",
                                                            str_detect(motif_alt_id,"CRF") ~"CRF")) %>% na.omit() %>%
                                                            dplyr::select(sequence_name,MotifGroup) %>% unique()

unique(fimoUP_edit$MotifGroup)
is.na(fimoUP_edit)

test <- fimoUP_edit[is.na(fimoUP_edit$MotifGroup),]

table(fimoUP_edit$MotifGroup)
table(fimoDOWN_edit$MotifGroup)

bind_rows(
  fimoUP_edit %>% mutate(Direction = "UP"),
  fimoDOWN_edit %>% mutate(Direction = "DOWN")
) %>% dplyr::filter(!str_detect(MotifGroup,"BPC")) %>%
  ggplot(aes(x = MotifGroup, fill = Direction)) +
  geom_bar(position = position_dodge(width = 0.9), width = 0.7) +
  labs(y = "Number of DEGs containing described TFS binding site", x = "MotifGroup")+
  theme_classic() +
  theme(axis.text = element_text(size = 14, hjust = 1, face = "bold"),
        axis.text.x = element_text(angle = 45),
        legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        axis.title  = element_text(size = 16, face = "bold"))

bind_rows(
  fimoUP_edit %>% mutate(Direction = "UP"),
  fimoDOWN_edit %>% mutate(Direction = "DOWN")
) %>% #dplyr::filter(!str_detect(MotifGroup,"BPC")) %>%
  ggplot(aes(x = MotifGroup, fill = Direction)) +
  geom_bar(position = position_dodge(width = 0.9), width = 0.7) +
  labs(y = "Number of DEGs containing described TFS binding site", x = "MotifGroup")+
  theme_classic() +
  theme(axis.text = element_text(size = 14, hjust = 1, face = "bold"),
        axis.text.x = element_text(angle = 45),
        legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        axis.title  = element_text(size = 16, face = "bold"))

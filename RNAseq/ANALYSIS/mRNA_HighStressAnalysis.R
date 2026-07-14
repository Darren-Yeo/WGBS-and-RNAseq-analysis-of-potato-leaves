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

#edit the annotation table to get gene ID to AT symbol
GeneID2GeneLabel<- S.tuberosum.v6.1_latest %>% dplyr::select(locusName,v6.1.Description,aktualisierte.annotation..bitte.ergänzen.,Mapman.description)


#GeneID2ATSymbol <- GeneID2ATSymbol %>%
#  mutate(
#    AT_ID = str_extract(Mapman.description, "AT[1-5]G\\d{5}"),
#    Symbols = str_extract(Mapman.description, "(?<=Symbols: )(.*?)(?= \\|)")  # Capture everything after "Symbols: " until " |"
#  )
#

# Merge v6.1 description and actual description into a new column
GeneID2GeneLabel <- GeneID2GeneLabel %>%
  mutate(
    Gene_label = paste(v6.1.Description, 
                       aktualisierte.annotation..bitte.ergänzen.,
                       locusName, sep = " | ")
  )

GeneID2GeneLabel$Gene_label <- gsub("Soltu.DM.","",GeneID2GeneLabel$Gene_label)

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


###Analysis with interacting terms.. conditions as the main comparison
dds <-  DESeqDataSetFromMatrix(countData=UniqueCountsDuplicatedHighStress_counts, 
          colData=MetaData, 
          design= ~ Conditions + Genotype + Conditions:Genotype)

#Run DEseq..
dds_full <- DESeq(dds)

#Get normalized counts
##get vst counts
dds_norm <- counts(dds_full, normalized =TRUE)



##get vst counts
dds_vst<- vst(dds_full)
dds_vst <- assay(dds_vst)

#Scale dataset, calculate Z-scores..
dds_vst_Scaled<- dds_vst %>% t() %>% scale() %>% t() %>%
  as.data.frame() %>%
  rownames_to_column("locusName")
###save vst and norm counts
write.table(dds_vst, "/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/06_ANALYSIS/EpiPotato_mRNA_HighStress/dds_vst_highHeat.txt",col.names = TRUE, row.names = TRUE)
write.table(dds_vst_Scaled, "/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/06_ANALYSIS/EpiPotato_mRNA_HighStress/dds_vst_Scaled_highHeat.txt",col.names = TRUE, row.names = TRUE)
write.table(dds_norm, "/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/06_ANALYSIS/EpiPotato_mRNA_HighStress/dds_norm_highHeat.txt",col.names = TRUE, row.names = TRUE)
write.table(MetaData, "/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/06_ANALYSIS/EpiPotato_mRNA_HighStress/MetaData.txt",col.names = TRUE, row.names = TRUE)

PCA <- PCA(t(dds_vst), scale.unit = TRUE, graph = FALSE)

fviz_pca_ind(PCA, 
             addEllipses = FALSE, label = "var",
             col.var = "black", repel = TRUE) + 
  geom_point(aes(color=MetaData$Conditions), size = 3) +
  geom_text(aes(label = MetaData$Genotype), 
            size = 4, vjust = 1, hjust = 1, color = "black",
            fontface= "bold")+
  labs(color="Conditions")+
  #scale_shape_manual(values=c(19,17))+
  scale_colour_manual(values=c("Heat"="#CC0000",
                               "Control"="#66CCFF"))+
  labs(title = "PCA projection of transcriptome") +
  #labs(x = "PC1:18.9%",
  #     y = "PC2:15.1%",
  #     subtitle="PCA transcriptome")+
  theme_bw()+
  theme(legend.position="bottom",
        plot.title = element_text(size = 20,face = "bold"),
        legend.title = element_text(size = 16,face = "bold"),
        legend.text = element_text(size = 16,face = "bold"),
        axis.title = element_text(size = 16,face="bold"),
        axis.text = element_text(size = 16,face = "bold"))

fviz_pca_ind(PCA, 
             addEllipses = FALSE, label = "var",
             col.var = "black", repel = TRUE) + 
  geom_point(aes(shape = MetaData$Conditions), size = 3) +
  geom_text(aes(label = MetaData$Genotype), 
            size = 4, vjust = 1, hjust = 1, color = "black",
            fontface= "bold") +
  labs(shape = "Conditions") +
  scale_shape_manual(values = c("Heat" = 17, "Control" = 19)) +  # 17 = triangle, 19 = circle
  labs(title = "PCA projection of transcriptome") +
  theme_bw() +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"))

#met1
#Soltu.DM.11G013230
#Soltu.DM.04G014710

#CMT3
Soltu.DM.12G000130

##DRM1
Soltu.DM.10G030090

dds_Zscore_selected_long<- ConvertWide2LongGGplot2_selectedGenes(c(
                                        "Soltu.DM.05G026370",
                                        "Soltu.DM.02G030260",
                                        "Soltu.DM.02G030280",
                                        "Soltu.DM.02G030300",
                                        "Soltu.DM.05G024030",
                                        "Soltu.DM.03G011110",
                                        "Soltu.DM.06G029500",
                                        "Soltu.DM.05G005140",
                                        "Soltu.DM.09G004240",
                                        "Soltu.DM.10G024770",
                                        "Soltu.DM.11G005260",
                                        "Soltu.DM.04G021680",
                                        "Soltu.DM.04G017430"),as.data.frame(t(scale(t(dds_vst)))), "Z_score")

dds_Zscore_selected_long_label <- inner_join(dds_Zscore_selected_long,
                                           S.tuberosum.v6.1_latest %>% dplyr::select(locusName,aktualisierte.annotation..bitte.ergänzen.))

##plot the expression of key genes
genotype_order <-c("AA_C","AA_H","CA_C","CA_H") 
pdf("/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/06_ANALYSIS/ZscoreExpressionSelectedGenes.pdf",
    width = 16, height = 12)

for (genes in c("Soltu.DM.05G026370",
                "Soltu.DM.02G030260",
                "Soltu.DM.02G030280",
                "Soltu.DM.02G030300",
                "Soltu.DM.05G024030",
                "Soltu.DM.03G011110",
                "Soltu.DM.06G029500",
                "Soltu.DM.05G005140")) {
  
  print(paste0("Plotting expression for gene ",genes))
  
  plot <- ggplot() +
    geom_boxplot(data = dds_Zscore_selected_long_label %>% filter(str_detect(locusName,genes)), 
                 aes(x = factor(Genotype, levels = genotype_order), y = Z_score, fill = Conditions),
                 position = position_dodge(1), width = 0.4) +
    facet_wrap(~aktualisierte.annotation..bitte.ergänzen., nrow = 3) +
    scale_fill_manual(values=c("lightblue","#B5251C")) + 
    ylab("Norm Counts") +
    xlab("Genotypes (Sensitive to tolerant)") +
    theme_bw() +
    theme(axis.text = element_text(size = 10,face="bold"),
          strip.text = element_text(size = 10,face="bold"),
          axis.title = element_text(size = 14,face="bold"),
          legend.title = element_text(size = 14,face="bold"),
          legend.text = element_text(size = 10,face="bold"))
  
  print(plot)
  
}

dev.off()

SP6A_plot <-  dds_Zscore_selected_long_label %>% filter(str_detect(locusName,"Soltu.DM.05G026370"))
ROS1_plot <-  dds_Zscore_selected_long_label %>% filter(str_detect(locusName,"Soltu.DM.10G024770"))

plot <- ggplot() +
  geom_boxplot(data = dds_Zscore_selected_long_label %>% filter(str_detect(locusName,"Soltu.DM.05G026370")), 
               aes(x = factor(Genotype, levels = genotype_order), y = Z_score, fill = Genotype),
               position = position_dodge(1), width = 0.4) +
  facet_wrap(~aktualisierte.annotation..bitte.ergänzen., 
             labeller = labeller(aktualisierte.annotation..bitte.ergänzen. = c("StPEBP9, SP6A" = "SP6A - Soltu.DM.05G026370"))) +
  scale_fill_manual(values = c(
    "AA_C" = "lightblue3",
    "AA_H" = "blue",
    "CA_C" = "pink3",
    "CA_H" = "red4"
  )) +
  ylab("Expression (Z scores)") +
  xlab("Genotypes") +
  scale_x_discrete(labels = c(
    "AA_C" = "AA",
    "AA_H" = "AA37C",
    "CA_C" = "CA",
    "CA_H" = "CA37C"
  )) +
  theme_bw() +
  theme(title = element_text(face = "bold",size = 25),
        legend.title = element_text(face = "bold",size = 18),
        legend.text = element_text(face = "bold",size = 16),
        legend.position = "none",
        axis.title = element_text(face = "bold",size = 18),
        axis.text = element_text(face = "bold",size=16),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,size=16),
        strip.text =element_text(face = "bold",size=16))

(plot2 <- ggplot() +
  geom_boxplot(data = dds_Zscore_selected_long_label %>% 
                 filter(str_detect(locusName,"Soltu.DM.10G024770|Soltu.DM.09G004240|Soltu.DM.11G005260")), 
               aes(x = factor(Genotype, levels = genotype_order), y = Z_score, fill = Genotype),
               position = position_dodge(1), width = 0.4) +
  facet_wrap(~locusName) +
  scale_fill_manual(values = c(
    "AA_C" = "lightblue3",
    "AA_H" = "blue",
    "CA_C" = "pink3",
    "CA_H" = "red4"
  )) +
  ylab("Expression (Z scores)") +
  xlab("Genotypes") +
  scale_x_discrete(labels = c(
    "AA_C" = "AA",
    "AA_H" = "AA37C",
    "CA_C" = "CA",
    "CA_H" = "CA37C"
  )) +
  theme_bw() +
  theme(title = element_text(face = "bold",size = 25),
        legend.title = element_text(face = "bold",size = 18),
        legend.text = element_text(face = "bold",size = 16),
        legend.position = "none",
        axis.title = element_text(face = "bold",size = 18),
        axis.text = element_text(face = "bold",size=16),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,size=16),
        strip.text =element_text(face = "bold",size=16)))

ggsave("SP6A_plot_colors.svg",
       plot = plot,
       device = "svg",
       path = "/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/06_ANALYSIS/EpiPotato_mRNA_HighStress/",
       width = 8,
       height = 8,
       dpi = 300)

### methylase and demethylase
dds_Zscore_selected_long_metDemet<- ConvertWide2LongGGplot2_selectedGenes(c(
  "Soltu.DM.04G014710",
  "Soltu.DM.11G013230",
  "Soltu.DM.01G001630",
  "Soltu.DM.08G003560",
  "Soltu.DM.12G000130",
  "Soltu.DM.08G015570",
  "Soltu.DM.08G015580",
  "Soltu.DM.02G006560",
  "Soltu.DM.04G000230",
  "Soltu.DM.10G030090",
  "Soltu.DM.11G005260",
  "Soltu.DM.09G004240",
  "Soltu.DM.10G024770",
  "Soltu.DM.03G037240",
  "Soltu.DM.04G021680",
  "Soltu.DM.06G028750",
  "Soltu.DM.09G025490",
  "Soltu.DM.03G027780"),as.data.frame(t(scale(t(dds_vst)))), "Z_score")

dds_Zscore_selected_long_metDemet_label <- inner_join(dds_Zscore_selected_long_metDemet,
                                             S.tuberosum.v6.1_latest %>% dplyr::select(locusName,v6.1.Description))

ggplot() +
  geom_boxplot(data = dds_Zscore_selected_long_metDemet_label %>% 
                 filter(str_detect(locusName,"Soltu.DM.11G005260|Soltu.DM.09G004240|Soltu.DM.10G024770|Soltu.DM.03G037240|Soltu.DM.04G021680|Soltu.DM.06G028750|Soltu.DM.09G025490|Soltu.DM.03G027780")), 
               aes(x = factor(Genotype, levels = genotype_order), y = Z_score, fill = Genotype),
               position = position_dodge(1), width = 0.4) +
  facet_wrap(~locusName) +
  scale_fill_manual(values = c(
    "AA_C" = "lightblue3",
    "AA_H" = "blue",
    "CA_C" = "pink3",
    "CA_H" = "red4"
  )) +
  ylab("Expression (Z scores)") +
  xlab("Genotypes") +
  scale_x_discrete(labels = c(
    "AA_C" = "AA",
    "AA_H" = "AA37C",
    "CA_C" = "CA",
    "CA_H" = "CA37C"
  )) +
  theme_bw() +
  theme(title = element_text(face = "bold",size = 25),
        legend.title = element_text(face = "bold",size = 18),
        legend.text = element_text(face = "bold",size = 16),
        legend.position = "none",
        axis.title = element_text(face = "bold",size = 18),
        axis.text = element_text(face = "bold",size=16),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,size=16),
        strip.text =element_text(face = "bold",size=16))


###plot overlay of constans 1,2,3 and sp6a


ggsave("SP6A_CONSTANS.svg",
         plot = ggplot(dds_Zscore_selected_long_label %>% 
                         filter(str_detect(locusName,"Soltu.DM.05G026370|Soltu.DM.02G030260|Soltu.DM.02G030280|Soltu.DM.02G030300")), 
                       aes(x = factor(Genotype, levels = genotype_order), y = Z_score, 
                           fill = aktualisierte.annotation..bitte.ergänzen.)) +
         scale_fill_manual(values=c("red4","orange3","yellow3","lightblue3")) +
         geom_boxplot(width = 0.6, position = position_dodge()) +  # Adjust width and transparency
         ylab("Expression (Z scores)") +
         xlab("Sample names") +
         labs(fill = "Gene") +
         theme_bw() +
         theme(axis.text = element_text(size = 13,face="bold"),
               strip.text = element_text(size = 10,face="bold"),
               axis.title = element_text(size = 18,face="bold"),
               axis.title.x = element_blank(),
               axis.text.x = element_blank(),
               legend.title = element_text(size = 18,face="bold"),
               legend.text = element_text(size = 13,face="bold")),
         device = "svg",
         path = "/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/06_ANALYSIS/EpiPotato_mRNA_HighStress/",
         width = 10,
         height = 8,
         dpi = 300)

###run reseq, reduced model
# Get the reduced model without the interaction term
dds_red <- DESeq(dds, test = "LRT", reduced = ~ Conditions + Genotype)
# Get the results for the reduced model
res_inter <- results(dds_red)
row.names(res_inter) <- gsub(".v6.1","",row.names(res_inter))
res_inter_labelled<- res_inter %>% 
  as.data.frame() %>% 
  rownames_to_column("locusName") %>% 
  inner_join(S.tuberosum.v6.1_latest)

##get only significant genes
res_inter_labelled_sig <- res_inter_labelled %>% dplyr::filter(padj <=0.05)

#get L2FC for the significant genes
res_inter_labelled_sig_L2FC <- res_inter_labelled_sig$log2FoldChange

#label names to respective L2FC 
names(res_inter_labelled_sig_L2FC) <-res_inter_labelled_sig$locusName

#omit the NA values
res_inter_labelled_sig_L2FC <- na.omit(res_inter_labelled_sig_L2FC)

#Sort from highest to lowest
res_inter_labelled_sig_L2FC <- sort(res_inter_labelled_sig_L2FC,decreasing = TRUE)

##do GSEA

##Compare between genotypes

###Analysis without interacting terms.. conditions as the main comparison


###conditions comparison
##only four significant results
GSEA_res_inter<- RunGSEA(res_inter_labelled_sig_L2FC,
                          term2gene = term2gene_GO_DB,
                          term2name = term2name_GO_DB)

GSEA_res_inter_df <- GSEA_res_inter@result

##ORA
#separate upregulated and downregulated genes
res_inter_labelled_sig_UP<- res_inter_labelled_sig_L2FC[res_inter_labelled_sig_L2FC > 1] 
res_inter_labelled_sig_DOWN<- res_inter_labelled_sig_L2FC[res_inter_labelled_sig_L2FC < 1] 

ORA_res_inter_list<- GO_KO_enrichment_with_enricher_CP(
  list(
    upregulated = names(res_inter_labelled_sig_UP),
    downregulated = names(res_inter_labelled_sig_DOWN)),
  term2gene_GO=term2gene_GO_DB,
  term2name_GO=term2name_GO_DB
)

ORA_res_inter_list$upregulated@result

#Conditions comparison
dds_list<- list(
dds_AA_HvsC= DESeqDataSetFromMatrix(countData=UniqueCountsDuplicatedHighStress_counts %>% as.data.frame() %>%
                         dplyr::select(contains("15"),contains("4"),contains("5")) %>% as.matrix(), 
                       colData=MetaData %>% dplyr::filter(str_detect(Genotype,"AA")), 
                       design= ~ Conditions),

dds_CA_HvsC= DESeqDataSetFromMatrix(countData=UniqueCountsDuplicatedHighStress_counts %>% as.data.frame() %>%
                         dplyr::select(contains("11"),contains("2"),contains("3")) %>% as.matrix(), 
                       colData=MetaData %>% dplyr::filter(str_detect(Genotype,"CA")), 
                       design= ~ Conditions),

#genoytpe comparison
dds_C_AAvsCA= DESeqDataSetFromMatrix(countData=UniqueCountsDuplicatedHighStress_counts %>% as.data.frame() %>%
                         dplyr::select(contains("K")) %>% as.matrix(), 
                       colData=MetaData %>% dplyr::filter(str_detect(Conditions,"Control")), 
                       design= ~ Genotype),

dds_H_AAvsCA= DESeqDataSetFromMatrix(countData=UniqueCountsDuplicatedHighStress_counts %>% as.data.frame() %>%
                         dplyr::select(contains("H")) %>% as.matrix(), 
                       colData=MetaData %>% dplyr::filter(str_detect(Conditions,"Heat")), 
                       design= ~ Genotype)

)

##Get results from the comparisons
dds_list_res <- list()
for (dds_name in names(dds_list)) {
  
  #run the DEseq function
  dds <- DESeq(dds_list[[dds_name]])
  
  design_formula <- as.character(dds@design)
  
  #get the results
  if ("Conditions" %in% all.vars(dds@design)){
     
    res <- results(dds, contrast = c("Conditions","Heat","Control"), pAdjustMethod = "fdr")
    
  } else if ("Genotype" %in% all.vars(dds@design)) {
    
    res <- results(dds, contrast = c("Genotype","AA","CA"), pAdjustMethod = "fdr")
    
  } else {
    next
  }
  
  rownames(res) <- gsub(".v6.1","",rownames(res))
  
  ##label the genes
  dds_list_res[[dds_name]] <- res %>% 
    as.data.frame() %>% 
    rownames_to_column("locusName") %>%
    inner_join(GeneID2GeneLabel)
  
}

dds_list_res_sig <- list()
for (dds_name in names(dds_list_res)) {
  
  dds_list_res_sig[[dds_name]] <- dds_list_res[[dds_name]] %>% dplyr::filter(padj <= 0.05)
  
}

dds_DEGs_Counts<- data.frame(
  Comparison = c("AA_norm", "CA_norm" , "AA_high", "CA_high"),
  Upregulated = c(2887,1967,4644,4368),
  Downregulated = c(1792,1237,3742,3605)
)

dds_DEGs_Counts_long <- dds_DEGs_Counts %>%
  pivot_longer(cols = c(Upregulated, Downregulated), 
               names_to = "Regulation", 
               values_to = "Count") %>% dplyr::mutate(
                 Batch=case_when(str_detect(Comparison,"norm") ~ paste0("30", "\u00B0", "C","/","28","\u00B0", "C"),
                                 str_detect(Comparison,"high") ~ paste0("37", "\u00B0", "C","/","35","\u00B0", "C")),
                 Genotype=case_when(str_detect(Comparison,"AA")~ "AA",
                                    str_detect(Comparison,"CA")~ "CA"))

ggplot(
  dds_DEGs_Counts_long %>% dplyr::filter(str_detect(Batch,"30")),aes(x = Genotype, y = Count,
                           fill = Regulation)) + 
  geom_bar(stat = "identity", position = "dodge") +  # Side-by-side bars
  facet_wrap(~ Batch) + 
  geom_text(aes(label = Count), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, 
            size = 5, fontface = "bold") + 
  scale_y_continuous(limits = c(0, 5000))+
  scale_fill_manual(values = c("Upregulated" = "darkgreen", 
                               "Downregulated" = "purple3")) +
  ylab("DEG counts") +
  theme_bw() +
  theme(axis.text = element_text(size = 20,face="bold"),
        strip.text = element_text(size = 20,face="bold"),
        axis.title = element_blank(),
        legend.title = element_text(size = 20,face="bold"),
        legend.text = element_text(size = 18,face="bold"))



dds_list_res_sig_L2FC <- list()
for (dds_name in names(dds_list_res_sig)) {
  
  ##
  L2FC_genes <- dds_list_res_sig[[dds_name]]$log2FoldChange
  
  names(L2FC_genes) <- dds_list_res_sig[[dds_name]]$locusName
  
  L2FC_genes <- na.omit(L2FC_genes)
  
  L2FC_genes <- sort(L2FC_genes,decreasing = TRUE)
  
  dds_list_res_sig_L2FC[[dds_name]] <- L2FC_genes
  
}


RunGSEA<- function(L2FC_decreasing_GeneNames, term2gene, term2name,
                   exponent = 1,
                   minGSSize = 10,
                   maxGSSize = 500,
                   eps = 0,
                   pvalueCutoff = 0.05, 
                   pAdjustMethod = "fdr",
                   gson = NULL,
                   seed = FALSE,
                   by = "fgsea"){
  
  GSEA_results<- GSEA(
    L2FC_decreasing_GeneNames,
    exponent = exponent,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    eps = eps, #For some pathways, in reality P-values are less than 1e-10. You can set the `eps` argument to zero for better estimation.
    pvalueCutoff = pvalueCutoff,
    pAdjustMethod = pAdjustMethod,
    gson = gson,
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    verbose = TRUE,
    seed = seed,
    by = by,
    #There were 1 pathways for which P-values were not calculated properly due to unbalanced (positive and negative) 
    #gene-level statistic values. For such pathways pval, padj, NES, log2err are set to NA.
    #You can try to increase the value of the argument nPermSimple (for example set it nPermSimple = 10000)
    nPermSimple = 10000)
  
  return(GSEA_results)
}


dds_list_res_GSEA <- lapply(dds_list_res_sig_L2FC, RunGSEA,
       term2gene = term2gene_GO_DB,
       term2name = term2name_GO_DB)

##tried to make labels for the gene ID, not successful, mapping was all wrong
#id2symbol <- setNames(dds_list_res_sig$dds_AA_HvsC$Gene_label, 
#                      dds_list_res_sig$dds_AA_HvsC$locusName)
#
#dds_list_res_GSEA$dds_AA_HvsC@result$core_enrichment <- sapply(strsplit(dds_list_res_GSEA$dds_AA_HvsC@result$core_enrichment, "/"), function(genes) {
#  paste(na.omit(id2symbol[genes]), collapse = "/")
#})


test<- dds_list_res_GSEA$dds_AA_HvsC@result
#plot the L2FC with the enrichment AA HvsC
gseaplot(dds_list_res_GSEA$dds_AA_HvsC,
         "GO:0009408",
         title = "response to heat")

gseaplot2(dds_list_res_GSEA$dds_AA_HvsC,
         geneSetID = c("GO:0009408","GO:0006457","GO:0042542","GO:0009768","GO:0009834","GO:0009755"),
         base_size = 14,
         pvalue_table = FALSE,
         rel_heights = c(1, 0.15, 0.2),
         subplots = 1:3)


gseaplot(dds_list_res_GSEA$dds_AA_HvsC,
         "GO:0010345",
         title = "suberin biosynthetic process")

View(dds_list_res_GSEA$dds_AA_HvsC@result)

###treeplot of enriched terms
treeplot(pairwise_termsim(dds_list_res_GSEA$dds_AA_HvsC),
         method = "JC")
treeplot(pairwise_termsim(dds_list_res_GSEA$dds_CA_HvsC),
         method = "JC")


emapplot(pairwise_termsim(dds_list_res_GSEA$dds_AA_HvsC),
         showCategory = 20)

emapplot(pairwise_termsim(dds_list_res_GSEA$dds_CA_HvsC),
         showCategory = 20)

ridgeplot(dds_list_res_GSEA$dds_AA_HvsC,
          showCategory = 20,
          label_format = 100)

ridgeplot(dds_list_res_GSEA$dds_CA_HvsC,
          showCategory = 20,
          label_format = 100)

heatplot(dds_list_res_GSEA$dds_AA_HvsC, 
         foldChange=dds_list_res_sig_L2FC$dds_AA_HvsC, showCategory=20) +
  scale_fill_gradientn(colors = c("blue", "white", "red"))

#plot the L2FC with the enrichment CA HvsC
heatplot(dds_list_res_GSEA$dds_CA_HvsC, 
         foldChange=dds_list_res_sig_L2FC$dds_AA_HvsC, showCategory=15) +
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  theme(axis.text.x = element_blank())


##GET the core enrichment genes


library(openxlsx)
dds_list_res_GSEA_functannot_unique$dds_AA_HvsC
dds_list_res_GSEA_functannot_unique$dds_CA_HvsC


dds_list_res_GSEA_functannot_unique_HvsC_L2FC_list <- list(
#Merge respective L2FC to respective comparisons for the GSEA table (ANNABELLE)
Annabelle_HvsC = inner_join(dds_list_res_GSEA_functannot_unique$dds_AA_HvsC %>% dplyr::select(-pvalue,-p.adjust), dds_list_res_sig$dds_AA_HvsC %>% dplyr::select(locusName,log2FoldChange),by="locusName"),
#Merge respective L2FC to respective comparisons for the GSEA table (CAMEL)
Camel_HvsC = inner_join(dds_list_res_GSEA_functannot_unique$dds_CA_HvsC %>% dplyr::select(-pvalue,-p.adjust), dds_list_res_sig$dds_CA_HvsC %>% dplyr::select(locusName,log2FoldChange),"locusName")

)

#rearrange table
for (tables in names(dds_list_res_GSEA_functannot_unique_HvsC_L2FC_list)) {

  #rearrange table
  dds_list_res_GSEA_functannot_unique_HvsC_L2FC_list[[tables]] <- dds_list_res_GSEA_functannot_unique_HvsC_L2FC_list[[tables]] %>% 
    dplyr::select(locusName,
                  v6.1.Description,
                  aktualisierte.annotation..bitte.ergänzen.,
                  log2FoldChange,
                  NES,
                  ID,
                  Description,
                  Mapman.description) 
}


# Function to save two data frames to an Excel file with separate sheets
save_dfs_to_excel <- function(df1, df2, df3, df4, df5, df6, SheetName1="Sheet1",
                              SheetName2="Sheet2",
                              SheetName3="Sheet3",
                              SheetName4="Sheet4",
                              SheetName5="Sheet5",
                              SheetName6="Sheet6",
                              file_name = "output.xlsx") {
  # Create a new workbook
  wb <- createWorkbook()
  
  # Add the first data frame to the first sheet
  addWorksheet(wb, SheetName1)
  writeData(wb, sheet = SheetName1, df1)
  
  # Check if the second data frame is not NULL before adding the second sheet
  if (!is.null(df2)) {
    addWorksheet(wb, SheetName2)
    writeData(wb, sheet = SheetName2, df2)
  }
  
  if (!is.null(df3)) {
    addWorksheet(wb, SheetName3)
    writeData(wb, sheet = SheetName3, df3)
  }
  
  if (!is.null(df4)) {
    addWorksheet(wb, SheetName4)
    writeData(wb, sheet = SheetName4, df4)
  }
  
  if (!is.null(df5)) {
    addWorksheet(wb, SheetName5)
    writeData(wb, sheet = SheetName5, df5)
  }
  
  if (!is.null(df6)) {
    addWorksheet(wb, SheetName6)
    writeData(wb, sheet = SheetName6, df6)
  }
  
  # Save the workbook to the specified file
  saveWorkbook(wb, file_name, overwrite = TRUE)
  
  # Print confirmation message
  cat("Data frames have been saved to", file_name, "\n")
}

###save the dfs into Excel
##save_dfs_to_excel(dds_list_res_GSEA_functannot_unique_HvsC_L2FC_list$Annabelle_HvsC,
##                  dds_list_res_GSEA_functannot_unique_HvsC_L2FC_list$Camel_HvsC,
##                  SheetName1 = "Annabelle_HvsC",
##                  SheetName2 = "Camel_HvsC",
##                  file_name = "/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/06_ANALYSIS/test.xlsx")
##
##dds_list_res_GSEA_functannot_unique$dds_AA_HvsC %>% dplyr::select(locusName,v6.1.Description,aktualisierte.annotation..bitte.ergänzen.,NES,Mapman.description)
##dds_list_res_GSEA_functannot_unique$dds_CA_HvsC %>% dplyr::select(locusName,NES,Mapman.description)
##
###join both CA and AA together
##Intersection_GSEA_AA_CA_HvsC <- inner_join(dds_list_res_GSEA_functannot_unique$dds_AA_HvsC %>% ungroup() %>% dplyr::select(locusName,v6.1.Description,Description,aktualisierte.annotation..bitte.ergänzen.,NES,Mapman.description),
##                   dds_list_res_GSEA_functannot_unique$dds_CA_HvsC %>% ungroup() %>% dplyr::select(locusName,Description,NES),
##                   by="locusName")
##
###sort the NES by - and =, by character sorting
##Intersection_GSEA_AA_CA_HvsC<- rbind(Intersection_GSEA_AA_CA_HvsC %>% dplyr::filter(str_detect(NES.x,"\\-")),
##                                    Intersection_GSEA_AA_CA_HvsC %>% dplyr::filter(!str_detect(NES.x,"\\-")))
##
##Intersection_GSEA_AA_CA_HvsC <- Intersection_GSEA_AA_CA_HvsC %>% dplyr::select(locusName,Description.x,NES.x,Description.y,NES.y,v6.1.Description,aktualisierte.annotation..bitte.ergänzen.,Mapman.description)
##
###save to excel and send
##save_dfs_to_excel(df1 = Intersection_GSEA_AA_CA_HvsC %>% dplyr::select(locusName,Description.x,NES.x,NES.y,v6.1.Description,aktualisierte.annotation..bitte.ergänzen.,Mapman.description),
##                  df2 = NULL,
##                  SheetName1 = "GSEA_Annabelle_Camel_overlap",
##                  SheetName2 = NULL,
##                  file_name = "/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/06_ANALYSIS/Intersection_GSEA_AA_CA_HvsC.xlsx")


# get heatstress terms of genes from both CA and Annabelle
GSEA_HvsC_AA_CA_overlaps_heatStressTermsGenes <- Intersection_GSEA_AA_CA_HvsC %>% 
  dplyr::filter(
    str_detect(Description.x, "hydrogen peroxide|protein complex|protein folding|heat") |
      str_detect(Description.y, "hydrogen peroxide|protein complex|protein folding|heat")
  ) %>% 
  mutate(AT_Locus = str_extract(Mapman.description, "AT\\dG\\d{5}"))

save_dfs_to_excel(df1 = GSEA_HvsC_AA_CA_overlaps_heatStressTermsGenes %>% 
                    dplyr::select(locusName,AT_Locus,Description.x,NES.x,NES.y,v6.1.Description,aktualisierte.annotation..bitte.ergänzen.,Mapman.description),
                  df2 = NULL,
                  SheetName1 = "GSEA_AA_CA_HS_overlap",
                  SheetName2 = NULL,
                  file_name = "/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/06_ANALYSIS/GSEA_HvsC_AA_CA_overlaps_heatStressTermsGenes.xlsx")
devtools::install_github("dkainer/RWRtoolkit")


##GO enrichment with ORA




#upregulated genes
dds_list_res_sig_up <- list()
for (upregulated_genes in names(dds_list_res_sig_L2FC)) {
  
  dds_list_res_sig_up[[upregulated_genes]] <- names(dds_list_res_sig_L2FC[[upregulated_genes]][dds_list_res_sig_L2FC[[upregulated_genes]] >1 ])
  
}

#downregulated genes
dds_list_res_sig_down <- list()
for (downregulated_genes in names(dds_list_res_sig_L2FC)) {
  
  dds_list_res_sig_down[[downregulated_genes]] <- names(dds_list_res_sig_L2FC[[downregulated_genes]][dds_list_res_sig_L2FC[[downregulated_genes]] < -1 ])
  
}


##Running the ORA
ORA_dds_list_res_sig_up <- GO_KO_enrichment_with_enricher_CP(dds_list_res_sig_up, term2gene_GO = term2gene_GO_DB,term2name_GO = term2name_GO_DB)
ORA_dds_list_res_sig_down<- GO_KO_enrichment_with_enricher_CP(dds_list_res_sig_down, term2gene_GO = term2gene_GO_DB,term2name_GO = term2name_GO_DB)

test <- dds_list_res_sig_ORA_up$dds_AA_HvsC@result

dotplot(ORA_dds_list_res_sig_up$dds_AA_HvsC)
dotplot(ORA_dds_list_res_sig_up$dds_CA_HvsC)






############Compare between normal and high heat-stress############ Different experiment, may have batch effects.. not sure

UniqueCountsDuplicated_counts2 <- read.csv("/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/06_ANALYSES/RawCountsInput/UniqueCountsDuplicated_counts2.txt", 
                                           sep="")
## Metadata... 
MetaData_tolerance_Comnbined <- read.csv("/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/06_ANALYSES/RawCountsInput/MetaData_tolerance_Comnbined.txt", sep="")

##select the meta data for only the parents normal heat stress
MetaData_tolerance_Comnbined_parents <- MetaData_tolerance_Comnbined %>% 
  dplyr::filter(str_detect(Genotype,"AA|CA")) %>% 
  dplyr::select(sampleID,Conditions,Genotype)

## 
UniqueCountsDuplicated_counts2_normal_parents <- UniqueCountsDuplicated_counts2 %>% 
                                                  as.data.frame() %>% column_to_rownames("locusName") %>%
                                                  dplyr::select(contains("AA"),contains("CA"))




UniqueCountsDuplicatedHighStress_counts


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
                            str_detect(SampleID,"15H") ~ "AA"))

colnames(MetaData)[1] <- "sampleID"

##Combine the meta data together....
MetaDataOverallNormalHighheat <- rbind(
  MetaData_tolerance_Comnbined_parents,
  MetaData
)

MetaDataOverallNormalHighheat <-MetaDataOverallNormalHighheat %>% 
  dplyr::mutate( 
    Batch = case_when(str_detect(sampleID,"_C") ~ paste0("Batch_1 (","30", "\u00B0", "C","/","28","\u00B0", "C", ")"),
                      str_detect(sampleID,"_H") ~ paste0("Batch_1 (","30", "\u00B0", "C","/","28","\u00B0", "C", ")"),
                      str_detect(sampleID,"RNA") ~ paste0("Batch_2 (","37", "\u00B0", "C","/","35","\u00B0", "C", ")"))) 


UniqueCountsDuplicatedCombined <- cbind(UniqueCountsDuplicated_counts2_normal_parents,UniqueCountsDuplicatedHighStress_counts)

###Create DEseq2 object
dds_batch <-  DESeqDataSetFromMatrix(countData=UniqueCountsDuplicatedCombined, 
                                     colData=MetaDataOverallNormalHighheat, 
                                     design= ~ Batch + Conditions + Genotype + Conditions:Genotype)

#run Deseq model
dds_batch <- DESeq(dds_batch)
dds_batch_vst <- vst(dds_batch, blind = FALSE)
dds_batch_vst <- assay(dds_batch_vst)

###DO a PCA comparison to determine if there is batch effects... yes there was
PCA_batch <- PCA(t(dds_batch_vst), scale.unit = TRUE, graph = FALSE)

fviz_pca_ind(PCA_batch, 
             addEllipses = FALSE, label = "var",
             col.var = "black", repel = TRUE) + 
  geom_point(aes(color=MetaDataOverallNormalHighheat$Conditions, shape = MetaDataOverallNormalHighheat$Batch), size = 3) +
  geom_text(aes(label = MetaDataOverallNormalHighheat$Genotype), 
            size = 4, vjust = 1, hjust = 1, color = "black",fontface="bold")+
  labs(color="Conditions",
       shape= "Difference in heat stress temperatures")+
  scale_shape_manual(values=c(19,17))+
  scale_colour_manual(values=c("Heat"="#CC0000",
                               "Control"="#66CCFF"))+
  labs(title  ="PCA including heat-stress von two different experiments")+
  theme_bw()+
  theme(legend.position="bottom",
        plot.title = element_text(size = 20,face = "bold"),
        legend.title = element_text(size = 16,face = "bold"),
        legend.text = element_text(size = 16,face = "bold"),
        axis.title = element_text(size = 16,face="bold"),
        axis.text = element_text(size = 16,face = "bold"))


# Remove batch effect using removeBatchEffect
dds_batch_vst_corrected <- limma::removeBatchEffect(dds_batch_vst, batch = MetaDataOverallNormalHighheat$Batch, 
                                    design = model.matrix(~ MetaDataOverallNormalHighheat$Conditions + MetaDataOverallNormalHighheat$Genotype))

PCA_batch_corrected <- PCA(t(dds_batch_vst_corrected), scale.unit = TRUE, graph = FALSE)

fviz_pca_ind(PCA_batch_corrected, 
             addEllipses = FALSE, label = "var",
             col.var = "black", repel = TRUE) + 
  geom_point(aes(color=MetaDataOverallNormalHighheat$Conditions, shape = MetaDataOverallNormalHighheat$Batch), size = 3) +
  geom_text(aes(label = MetaDataOverallNormalHighheat$Genotype), 
            size = 4, vjust = 1, hjust = 1, color = "black",fontface="bold")+
  labs(color="Conditions",
       shape= "Batch and difference in heat stress temperatures")+
  scale_shape_manual(values=c(19,17))+
  scale_colour_manual(values=c("Heat"="#CC0000",
                               "Control"="#66CCFF"))+
  labs(title  ="PCA Including heat-stress von two different experiments (Batch corrected)")+
  theme_bw()+
  theme(legend.position="bottom",
        plot.title = element_text(size = 20,face = "bold"),
        legend.title = element_text(size = 16,face = "bold"),
        legend.text = element_text(size = 16,face = "bold"),
        axis.title = element_text(size = 16,face="bold"),
        axis.text = element_text(size = 16,face = "bold"))


##Do DEG analysis between parents


UniqueCountsDuplicated_counts2_normal_parents

MetaData_tolerance_Comnbined_parents

##run deseq for the parents, normal stress
dds_normalstress<- DESeqDataSetFromMatrix(countData = UniqueCountsDuplicated_counts2_normal_parents, 
                       colData = MetaData_tolerance_Comnbined_parents, 
                       design= ~ Conditions + Genotype + Conditions:Genotype)

dds_normalstress_full <- DESeq(dds_normalstress)

#get norm counts
dds_normalstress_norm <- counts(dds_normalstress_full,normalized=TRUE)
##get vst counts
dds_normalstress_vst<- vst(dds_normalstress_full)
dds_normalstress_vst <- assay(dds_normalstress_vst)

##filter by norm counts for the vst counts


PCA_normalstress <- PCA(t(dds_normalstress_vst), scale.unit = TRUE, graph = FALSE)
test <- MetaData_tolerance_Comnbined_parents 

fviz_pca_ind(PCA_normalstress, 
             addEllipses = FALSE, label = "var",
             col.var = "black", repel = TRUE) + 
  geom_point(aes(color=MetaData_tolerance_Comnbined_parents$Conditions), size = 3) +
  geom_text(aes(label = MetaData_tolerance_Comnbined_parents$sampleID), 
            size = 4, vjust = 1, hjust = 1, color = "black",
            fontface= "bold")+
  labs(color="Conditions")+
  #scale_shape_manual(values=c(19,17))+
  scale_colour_manual(values=c("Heat"="grey",
                               "Control"="grey"))+
  labs(title = "PCA projection of transcriptome") +
  #labs(x = "PC1:18.9%",
  #     y = "PC2:15.1%",
  #     subtitle="PCA transcriptome")+
  theme_bw()+
  theme(legend.position="bottom",
        plot.title = element_text(size = 20,face = "bold"),
        legend.title = element_text(size = 16,face = "bold"),
        legend.text = element_text(size = 16,face = "bold"),
        axis.title = element_text(size = 16,face="bold"),
        axis.text = element_text(size = 16,face = "bold"))

fviz_pca_ind(PCA_normalstress, 
             addEllipses = FALSE, label = "var",
             col.var = "black", repel = TRUE) + 
  geom_point(aes(shape = MetaData_tolerance_Comnbined_parents$Conditions), size = 3) +
  geom_text(aes(label = MetaData_tolerance_Comnbined_parents$sampleID), 
            size = 4, vjust = 1, hjust = 1, color = "black",
            fontface= "bold") +
  labs(shape = "Conditions") +
  scale_shape_manual(values = c("Heat" = 17, "Control" = 19)) +  # 17 = triangle, 19 = circle
  labs(title = "PCA projection of transcriptome") +
  theme_bw() +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"))


###Plot the COnstans1 and SP6A overlay but only with respect with annabelle and camel
ConvertWide2LongGGplot2_selectedGenes<- function(selected_gene_vector,Read_Counts_table, ZScore_or_NormCounts){
  
  ##Remove the extras
  row.names(Read_Counts_table) <- gsub(".v6.1","",row.names(Read_Counts_table))
  
  #Get counts of the selected genes
  Read_Counts_table_selectedGenes<- Read_Counts_table[row.names(Read_Counts_table) %in% selected_gene_vector,]
  
  ###Convert the table to long, only specific to my dataset
  Read_Counts_table_selectedGenes_long <- Read_Counts_table_selectedGenes %>% 
    rownames_to_column("locusName") %>% 
    pivot_longer(
      cols = c(contains("_C"),contains("_H")),
      names_to = "Sample",
      values_to = ZScore_or_NormCounts
    )  %>% 
    mutate(Conditions=case_when(str_detect(Sample,"C1|C2|C3|C4") ~ "Control",
                                str_detect(Sample,"H1|H2|H3|H4") ~ "Heat"),
           Genotype=case_when(str_detect(Sample,"130") ~ "130",
                              str_detect(Sample,"162") ~ "162",
                              str_detect(Sample,"239") ~ "239",
                              str_detect(Sample,"249") ~ "249",
                              str_detect(Sample,"327") ~ "327",
                              str_detect(Sample,"360") ~ "360",
                              str_detect(Sample,"365") ~ "365",
                              str_detect(Sample,"476") ~ "476",
                              str_detect(Sample,"70") ~ "70",
                              str_detect(Sample,"88") ~ "88",
                              str_detect(Sample,"AA") ~ "AA",
                              str_detect(Sample,"CA") ~ "CA"),
           Tolerance=case_when(str_detect(Sample,"AA") ~ "Tolerant",
                               str_detect(Sample,"CA") ~ "Less tolerant",
                               str_detect(Sample,"239") ~ "Tolerant",
                               str_detect(Sample,"249") ~ "Tolerant",
                               str_detect(Sample,"327") ~ "Tolerant",
                               str_detect(Sample,"360") ~ "Tolerant",
                               str_detect(Sample,"365") ~ "Tolerant",
                               str_detect(Sample,"70") ~ "Less tolerant",
                               str_detect(Sample,"88") ~ "Less tolerant",
                               str_detect(Sample,"130") ~ "Less tolerant",
                               str_detect(Sample,"162") ~ "Less tolerant",
                               str_detect(Sample,"476") ~ "Less tolerant"))
  
  
  Read_Counts_table_selectedGenes_long$Sample <- gsub("X","",Read_Counts_table_selectedGenes_long$Sample)
  
  Read_Counts_table_selectedGenes_long$Sample <- substring(Read_Counts_table_selectedGenes_long$Sample,
                                                           1, 
                                                           nchar(Read_Counts_table_selectedGenes_long$Sample) - 1)
  
  
  return(Read_Counts_table_selectedGenes_long)
}


rownames(dds_normalstress_vst) <- gsub(".v6.1","",rownames(dds_normalstress_vst))

CONSTANS_SP6A_long<- ConvertWide2LongGGplot2_selectedGenes(c("Soltu.DM.02G030260","Soltu.DM.02G030280","Soltu.DM.02G030300",
                                        "Soltu.DM.05G026370"),
                                      as.data.frame(t(scale(t(dds_normalstress_vst)))),
                                      "Z_score")

CONSTANS_SP6A_long <- CONSTANS_SP6A_long %>% inner_join(S.tuberosum.v6.1_latest,by="locusName")

genotype_order <- c("AA_C","AA_H","CA_C","CA_H")

ggplot(CONSTANS_SP6A_long %>% dplyr::filter(str_detect(locusName,"Soltu.DM.05G026370")), 
       aes(x = factor(Sample, levels = genotype_order), y = Z_score, 
           fill = Genotype)) +
  scale_fill_manual(values=c("lightblue","#B5251C")) +
  geom_boxplot(width = 0.6, position = position_dodge()) +  # Adjust width and transparency
  ylab("Expression (Z scores)") +
  xlab("Sample names") +
  labs(fill = "Gene") +
  theme_bw() +
  theme(axis.text = element_text(size = 13,face="bold"),
        strip.text = element_text(size = 10,face="bold"),
        axis.title = element_text(size = 18,face="bold"),
        axis.title.x = element_blank(),
        legend.title = element_text(size = 18,face="bold"),
        legend.text = element_text(size = 13,face="bold"))


###Create the DEseqmatrix again
##Conditions comparison

##RUNNN EVERYTHING AGAIN THIS TIME IN LIST.. to aovid mistakes again
dds_list<- list(
  AA_HvsC_Normal= DESeqDataSetFromMatrix(countData=UniqueCountsDuplicated_counts2_normal_parents %>% as.data.frame() %>%
                                        dplyr::select(contains("AA")) %>% as.matrix(), 
                                      colData=MetaData_tolerance_Comnbined_parents %>% 
                                        dplyr::filter(str_detect(Genotype,"AA")), 
                                      design= ~ Conditions),
  
  CA_HvsC_Normal= DESeqDataSetFromMatrix(countData=UniqueCountsDuplicated_counts2_normal_parents %>% as.data.frame() %>%
                                        dplyr::select(contains("CA")) %>% as.matrix(), 
                                      colData=MetaData_tolerance_Comnbined_parents %>% 
                                        dplyr::filter(str_detect(Genotype,"CA")), 
                                      design= ~ Conditions),
  
  AA_HvsC_High= DESeqDataSetFromMatrix(countData=UniqueCountsDuplicatedHighStress_counts %>% as.data.frame() %>%
                                        dplyr::select(contains("15"),contains("4"),contains("5")) %>% as.matrix(), 
                                      colData=MetaData %>% dplyr::filter(str_detect(Genotype,"AA")), 
                                      design= ~ Conditions),
  
  CA_HvsC_High= DESeqDataSetFromMatrix(countData=UniqueCountsDuplicatedHighStress_counts %>% as.data.frame() %>%
                                        dplyr::select(contains("11"),contains("2"),contains("3")) %>% as.matrix(), 
                                      colData=MetaData %>% dplyr::filter(str_detect(Genotype,"CA")), 
                                      design= ~ Conditions)

  
)

##Get results from the comparisons
dds_list_res <- list()
for (dds_name in names(dds_list)) {
  
  #run the DEseq function
  dds <- DESeq(dds_list[[dds_name]])
  
  #get the results
    
    res <- results(dds, contrast = c("Conditions","Heat","Control"), pAdjustMethod = "fdr")
  
  rownames(res) <- gsub(".v6.1","",rownames(res))
  
  ##label the genes
  dds_list_res[[dds_name]] <- res %>% 
    as.data.frame() %>% 
    rownames_to_column("locusName") %>%
    inner_join(GeneID2GeneLabel)
  
}
## Plot Volcano and label key genes like SP6A and COL1

test <- dds_list_res$AA_HvsC_Normal %>% dplyr::filter(padj <10^-110)
test <- dds_list_res$CA_HvsC_Normal %>% dplyr::filter(padj <10^-150)

annabelle_Normal_Labels <- dds_list_res$AA_HvsC_Normal %>% dplyr::mutate(
  Gene_label_vol=case_when(str_detect(locusName,"04G033650") ~ "cold regulated gene",
                           str_detect(locusName,"06G022960") ~ "subfamily IIIB acid phosphatase",
                           str_detect(locusName,"09G025880") ~ "Serine protease inhibitor",
                           str_detect(locusName,"05G026370") ~ "StPEBP9, SP6A") 
)

Camel_Normal_Labels <- dds_list_res$CA_HvsC_Normal %>% dplyr::mutate(
  Gene_label_vol=case_when(str_detect(locusName,"01G031960") ~ "sterol-4alpha-methyl oxidase 1-1",
                           str_detect(locusName,"01G048600") ~ "maternal effect embryo arrest",
                           str_detect(locusName,"03G021500") ~ "PRR",
                           str_detect(locusName,"04G026390") ~ "hypothetical",
                           str_detect(locusName,"06G003760") ~ "sterol 4-alpha-methyl-oxidase 2-2",
                           str_detect(locusName,"07G022840") ~ "SOUL heme",
                           str_detect(locusName,"10G026720") ~ "hypothetical",
                           str_detect(locusName,"05G026370") ~ "StPEBP9, SP6A")
)

###save to folder
write.table(dds_list_res$AA_HvsC_Normal,"/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/06_ANALYSIS/ExpressionDEGsFullTable/AA_HvsC_Normal.txt", col.names = TRUE, row.names = TRUE)
write.table(dds_list_res$CA_HvsC_Normal,"/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/06_ANALYSIS/ExpressionDEGsFullTable/CA_HvsC_Normal.txt", col.names = TRUE, row.names = TRUE)
write.table(dds_list_res$AA_HvsC_High,"/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/06_ANALYSIS/ExpressionDEGsFullTable/AA_HvsC_High.txt", col.names = TRUE, row.names = TRUE)
write.table(dds_list_res$CA_HvsC_High,"/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/06_ANALYSIS/ExpressionDEGsFullTable/CA_HvsC_High.txt", col.names = TRUE, row.names = TRUE)


CA_HvsC_High<- EnhancedVolcano::EnhancedVolcano(dds_list_res$CA_HvsC_High,
                 lab = dds_list_res$AA_HvsC_High$aktualisierte.annotation..bitte.ergänzen.,
                 selectLab = c("StPEBP9, SP6A"),
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
      legend.text =  element_text(size=16,),
      legend.title  =  element_text(size=16)
    )

AA_HvsC_High<- EnhancedVolcano::EnhancedVolcano(dds_list_res$AA_HvsC_High,
                                                lab = dds_list_res$AA_HvsC_High$aktualisierte.annotation..bitte.ergänzen.,
                                                selectLab = c("StPEBP9, SP6A"),
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
  )

?ggsave
ggsave(
  "AA_HvsC_High.svg",
  plot = AA_HvsC_High,
  device = "svg",
  path = "/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/06_ANALYSIS/EpiPotato_mRNA_HighStress/",
  width = 10,
  height = 10,
  dpi = 300
)



                                 


dds_list_res_sig <- list()
for (dds_name in names(dds_list_res)) {
  
  dds_list_res_sig[[dds_name]] <- dds_list_res[[dds_name]] %>% dplyr::filter(padj <= 0.05)
  
}


dds_list_res_sig$


##GET THE L2FC
dds_list_res_sig_L2FC <- list()
for (dds_name in names(dds_list_res_sig)) {
  
  ##
  L2FC_genes <- dds_list_res_sig[[dds_name]]$log2FoldChange
  
  names(L2FC_genes) <- dds_list_res_sig[[dds_name]]$locusName
  
  L2FC_genes <- na.omit(L2FC_genes)
  
  L2FC_genes <- sort(L2FC_genes,decreasing = TRUE)
  
  dds_list_res_sig_L2FC[[dds_name]] <- L2FC_genes
  
}


Convert2UpsetDF <- function(gene_list){
  # Get unique gene names
  
  
  for (comparisons in names(gene_list)) {
    
    gene_list[[comparisons]] <- gene_list[[comparisons]]$locusName
  }
  
  
  all_genes <- unique(unlist(gene_list))
  
  # Create an empty dataframe
  upset_df <- data.frame(Gene = all_genes)
  
  # Populate the dataframe with binary values
  for (comparison in names(gene_list)) {
    upset_df[[comparison]] <- ifelse(upset_df$Gene %in% gene_list[[comparison]], 1, 0)
  }
  
  return(upset_df)
  
}

##Separate all > 1 L2FC
dds_list_res_sig_Up <- list(
  AA_HvsC_Normal = dds_list_res_sig$AA_HvsC_Normal %>% dplyr::filter(log2FoldChange > 1),
  AA_HvsC_High = dds_list_res_sig$AA_HvsC_High %>% dplyr::filter(log2FoldChange > 1),
  CA_HvsC_Normal = dds_list_res_sig$CA_HvsC_Normal %>% dplyr::filter(log2FoldChange > 1),
  CA_HvsC_High = dds_list_res_sig$CA_HvsC_High %>% dplyr::filter(log2FoldChange > 1))

##Separate all < -1 L2FC
dds_list_res_sig_Down <- list(
  AA_HvsC_Normal = dds_list_res_sig$AA_HvsC_Normal %>% dplyr::filter(log2FoldChange < -1),
  AA_HvsC_High = dds_list_res_sig$AA_HvsC_High %>% dplyr::filter(log2FoldChange < -1),
  CA_HvsC_Normal = dds_list_res_sig$CA_HvsC_Normal %>% dplyr::filter(log2FoldChange < -1),
  CA_HvsC_High = dds_list_res_sig$CA_HvsC_High %>% dplyr::filter(log2FoldChange < -1))

###Venn diagramm of all DEGs between comparisons Annabelle and Camel (High heat)

DEGs_dds_highHeat_list_res_Sig <- list(
  AA_HvsC_High = dds_list_res_sig$AA_HvsC_High %>% dplyr::filter(log2FoldChange > 1 | log2FoldChange < -1),
    CA_HvsC_High = dds_list_res_sig$CA_HvsC_High %>% dplyr::filter(log2FoldChange > 1 | log2FoldChange < -1)
)

DEGs_dds_highHeat_list_res_Sig$AA_HvsC_High %>% dplyr::filter(log2FoldChange > 1 ) %>% nrow()
DEGs_dds_highHeat_list_res_Sig$CA_HvsC_High %>% dplyr::filter(log2FoldChange < -1 ) %>% nrow()
ggsave(
  "VennDiagrammHvC_HighHeat.svg",
  plot = ggvenn::ggvenn(
    list(
      Annabelle = DEGs_dds_highHeat_list_res_Sig$AA_HvsC_High$locusName,
      Camel = DEGs_dds_highHeat_list_res_Sig$CA_HvsC_High$locusName
    ),
    fill_color = c("orange4", "green4"),
    fill_alpha = 0.3,
    set_name_size = 11,
    text_size = 13
  ) ,
  device = "svg",
  path = "/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/06_ANALYSIS/EpiPotato_mRNA_HighStress/",
  width = 10,
  height = 10,
  dpi = 300
)

###Get the overlaps between Annabelle and Camel heat-stress reaction transcripts
UpregulatedCommon_AA_CA_HighHeat<- inner_join(dds_list_res_sig_Up$AA_HvsC_High %>% dplyr::select(locusName),
                                              dds_list_res_sig_Up$CA_HvsC_High %>% dplyr::select(locusName),
                                              by="locusName")


DownregulatedCommon_AA_CA_HighHeat<- inner_join(dds_list_res_sig_Down$AA_HvsC_High %>% dplyr::select(locusName),
                                                dds_list_res_sig_Down$CA_HvsC_High %>% dplyr::select(locusName),
                                                 by="locusName")
 
###GET the anti overlaps.... (unique to Annabelle)

UpregulatedANTI_AA_HighHeat <- anti_join(DEGs_dds_highHeat_list_res_Sig$AA_HvsC_High %>% dplyr::select(locusName,log2FoldChange),
                                          DEGs_dds_highHeat_list_res_Sig$CA_HvsC_High %>% dplyr::select(locusName,log2FoldChange),
                                        by = "locusName") %>% dplyr::filter(log2FoldChange >= 1)

DownregulatedANTI_AA_HighHeat <- anti_join(DEGs_dds_highHeat_list_res_Sig$AA_HvsC_High %>% dplyr::select(locusName,log2FoldChange),
                                           DEGs_dds_highHeat_list_res_Sig$CA_HvsC_High %>% dplyr::select(locusName,log2FoldChange),
                                          by = "locusName") %>% dplyr::filter(log2FoldChange <= -1)
##valdiate to venn diagramm
nrow(rbind(UpregulatedANTI_AA_HighHeat,
           DownregulatedANTI_AA_HighHeat))

###GET the anti overlaps.... (unique to Camel)
UpregulatedANTI_CA_HighHeat <- anti_join(DEGs_dds_highHeat_list_res_Sig$CA_HvsC_High %>% dplyr::select(locusName,log2FoldChange),
                                         DEGs_dds_highHeat_list_res_Sig$AA_HvsC_High %>% dplyr::select(locusName,log2FoldChange),
                                         by = "locusName") %>% dplyr::filter(log2FoldChange >= 1)

DownregulatedANTI_CA_HighHeat <- anti_join(DEGs_dds_highHeat_list_res_Sig$CA_HvsC_High %>% dplyr::select(locusName,log2FoldChange),
                                           DEGs_dds_highHeat_list_res_Sig$AA_HvsC_High %>% dplyr::select(locusName,log2FoldChange),
                                           by = "locusName") %>% dplyr::filter(log2FoldChange <= -1)

##valdiate to venn diagramm
nrow(rbind(UpregulatedANTI_CA_HighHeat,
     DownregulatedANTI_CA_HighHeat))
 
###Comparecluster of all sets... up and down regulated

CombinedEnrichmentSetsHighHeat_HvC_AA_CA<- compareCluster(
  list(
    Upregulated_AA_CA_= UpregulatedCommon_AA_CA_HighHeat$locusName,
    Downregulated_AA_CA = DownregulatedCommon_AA_CA_HighHeat$locusName,
    Upregulated_AA = UpregulatedANTI_AA_HighHeat$locusName,
    Downregulated_AA = DownregulatedANTI_AA_HighHeat$locusName,
    Upregulated_CA = UpregulatedANTI_CA_HighHeat$locusName,
    Downregulated_CA = DownregulatedANTI_CA_HighHeat$locusName
  ), 
  fun = "enricher",
  pvalueCutoff = 0.05,
  pAdjustMethod = "fdr",
  TERM2GENE = term2gene_GO_DB,
  TERM2NAME = term2name_GO_DB
  
)

col_fun <- colorRamp2(c(-2, 0, 2), c("#2A78C6","white","#C83B38"))
CombinedEnrichmentSetsHighHeat_HvC_AA_CA@compareClusterResult
##PLOT terms, top 10 by counts


ggsave(
  "GO_CommonUnique_AA_CA_HvC.svg",
  plot = dotplot(CombinedEnrichmentSetsHighHeat_HvC_AA_CA,
                 showCategory = 10,
                 by = "Count",
                 size = "Count") +
    scale_fill_gradient2(
      low = "#C83B38", mid = "white", high = "#2A78C6",
      midpoint = 0.025) +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20, face = "bold"),
          axis.text.y = element_text(size = 20, face = "bold"),
          legend.text = element_text(size =20, face ="bold"),
          legend.title = element_text(size =20, face ="bold")),
  device = "svg",
  path = "/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/06_ANALYSIS/EpiPotato_mRNA_HighStress/",
  width = 17,
  height = 22,
  dpi = 300
)
###Common and unique DEGs Annabelle and Camel under heat-stress GO enrichment
UncommonCommonDEGs_AA_CA_HvC<- GO_KO_enrichment_with_enricher_CP(
  list(
    Upregulated_AA_CA_= UpregulatedCommon_AA_CA_HighHeat$locusName,
    Downregulated_AA_CA = DownregulatedCommon_AA_CA_HighHeat$locusName,
    Upregulated_AA = UpregulatedANTI_AA_HighHeat$locusName,
    Downregulated_AA = DownregulatedANTI_AA_HighHeat$locusName,
    Upregulated_CA = UpregulatedANTI_CA_HighHeat$locusName,
    Downregulated_CA = DownregulatedANTI_CA_HighHeat$locusName
  ),
  term2gene_GO = term2gene_GO_DB,
  term2name_GO = term2name_GO_DB
)

#Functional annotaiton 
UncommonCommonDEGs_AA_CA_HvC_functannot <- list()
for (Intersections in names(UncommonCommonDEGs_AA_CA_HvC)) {
  
  #Separate table of both GO and KO Enrichment results and store to list
  GO_table <- UncommonCommonDEGs_AA_CA_HvC[[Intersections]]@result
  
  #rename the core_enrichment genes..
  colnames(GO_table)[colnames(GO_table) == "geneID"] <- "locusName"
  
  
  #Create the enrichment term table (summary)
  GO_table_summary<- GO_table %>% 
    dplyr::select(locusName,ID,Description,pvalue,p.adjust) %>% dplyr::filter(p.adjust < 0.05) %>%
    separate_rows(locusName, sep = "/") %>% left_join(S.tuberosum.v6.1_latest %>% 
                                                        dplyr::select(locusName,v6.1.Description,aktualisierte.annotation..bitte.ergänzen.,
                                                                      Mapman.description), by ="locusName")
  
  
  ##save everything to lsit
  UncommonCommonDEGs_AA_CA_HvC_functannot[[Intersections]] <- GO_table_summary
  
}


#####Ensure that each gene is unique...... and combine all GO terms and respective enrichment scores,pvalues etc.. of each gene together
UncommonCommonDEGs_AA_CA_HvC_functannot_unique <- list()
for (Intersections in names(UncommonCommonDEGs_AA_CA_HvC_functannot)) {
  
  #get the comparison name with both Pos and Neg table
  ORA_comparison <- UncommonCommonDEGs_AA_CA_HvC_functannot[[Intersections]]
  
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
  UncommonCommonDEGs_AA_CA_HvC_functannot_unique[[Intersections]] <- ORA_comparison_unique
  
}

##save file excel for johnanna
write.xlsx(UncommonCommonDEGs_AA_CA_HvC_functannot_unique, 
    file = "/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/06_ANALYSIS/EpiPotato_mRNA_HighStress/UncommonCommonDEGs_AA_CA_HvC_functannot_unique.xlsx")



###Common DEGs Annabelle and Camel under heat-stress
#####GO enrichment up and Downregulated genes high heat-stress. ###
CommonDEGs_AA_CA_HvC<- GO_KO_enrichment_with_enricher_CP(list(
  Upregulated = UpregulatedCommon_AA_CA_HighHeat$locusName,
 Downregulated= DownregulatedCommon_AA_CA_HighHeat$locusName
), 
term2gene_GO = term2gene_GO_DB,
term2name_GO = term2name_GO_DB)

CommonDEGs_AA_CA_HvC_TOP20<- CombineUpDownRegulation_SelectTopEnrichmentTermsByCounts_Plot(TopCounts = 10,
                                                                                           CommonDEGs_AA_CA_HvC$Upregulated@result,
                                                                                           CommonDEGs_AA_CA_HvC$Downregulated@result)



ggsave(
  "GO_Common_AA_CA_HvC.svg",
  plot = ggplot(CommonDEGs_AA_CA_HvC_TOP20,
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
  path = "/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/06_ANALYSIS/EpiPotato_mRNA_HighStress/",
  width = 14,
  height = 10,
  dpi = 300
)

GSEA_Combined <- compareCluster(
  dds_list_res_sig_L2FC, 
  fun = "GSEA",
  nPermSimple = 10000,
  pvalueCutoff = 0.05,
  pAdjustMethod = "fdr",
  TERM2GENE = term2gene_GO_DB,
  TERM2NAME = term2name_GO_DB
  
)


for (comparison in names(dds_list_res_sig_Up)) {
  
  print(paste(comparison,":",nrow(dds_list_res_sig_Up[[comparison]])))
  
}

for (comparison in names(dds_list_res_sig_Down)) {
  
  print(paste(comparison,":",nrow(dds_list_res_sig_Down[[comparison]])))
  
}

##convert the list to and upset DF
dds_list_res_sig_upset<- list(
upregulated_genes = Convert2UpsetDF(dds_list_res_sig_Up),
downregulated_genes = Convert2UpsetDF(dds_list_res_sig_Down))

#plot the upset plots (DEGs L2FC >1 and < -1)
UpSetR::upset(dds_list_res_sig_upset$upregulated_genes,
              text.scale = c(2, 2, 2, 2, 2, 2),
              sets = c("AA_HvsC_Normal",
                       "AA_HvsC_High",
                       "CA_HvsC_Normal",
                       "CA_HvsC_High"),
              keep.order = TRUE)
UpSetR::upset(dds_list_res_sig_upset$downregulated_genes,
              text.scale = c(2, 2, 2, 2, 2, 2),
              sets = c("AA_HvsC_Normal",
                       "AA_HvsC_High",
                       "CA_HvsC_Normal",
                       "CA_HvsC_High"),
              keep.order = TRUE)


####GET the selected intersections of the selected comparisons
DEGs_NormalHighHeatStress<- list(
  upregulated = list(
   Normal_heat_stress = micro.gen.extra::get_intersect_members(dds_list_res_sig_upset$upregulated_genes,
                                           c("AA_HvsC_Normal","CA_HvsC_Normal")),
   High_heat_stress = micro.gen.extra::get_intersect_members(dds_list_res_sig_upset$upregulated_genes,
                                           c("AA_HvsC_High","CA_HvsC_High"))
    
  ),
  dowregulated = list(
    Normal_heat_stress = micro.gen.extra::get_intersect_members(dds_list_res_sig_upset$downregulated_genes,
                                           c("AA_HvsC_Normal","CA_HvsC_Normal")),
    High_heat_stress = micro.gen.extra::get_intersect_members(dds_list_res_sig_upset$downregulated_genes,
                                           c("AA_HvsC_High","CA_HvsC_High"))
  )
)

##annotate them...
DEGs_NormalHighHeatStress_annot <- list()
for (Regulation in names(DEGs_NormalHighHeatStress)) {
  
  RegulationDirection<- DEGs_NormalHighHeatStress[[Regulation]]

  RegulationDirection_annot <- list()
  for (Heat_stress in names(RegulationDirection)) {
    
    ##convert each gene vector to dataframe
    GeneLocus <- as.data.frame(RegulationDirection[[Heat_stress]])
    
    #rename the column as locusName
    colnames(GeneLocus) <- "locusName"
    
    GeneLocus_annot<-  GeneLocus %>% 
      inner_join(
      S.tuberosum.v6.1_latest, by="locusName"
    )
    
    RegulationDirection_annot[[Heat_stress]] <- GeneLocus_annot
    
  }
  DEGs_NormalHighHeatStress_annot[[Regulation]] <- RegulationDirection_annot
}


## ORA enrichment, of the selected intersections (Common DEGs identified in both AA and CA under heatstress for normal HS and High HS)
ORA_NormalHigh_HS <- GO_KO_enrichment_with_enricher_CP(list(
      Normal_HS_up = DEGs_NormalHighHeatStress_annot$upregulated$Normal_heat_stress$locusName,
      High_HS_up = DEGs_NormalHighHeatStress_annot$upregulated$High_heat_stress$locusName,
      Normal_HS_down = DEGs_NormalHighHeatStress_annot$dowregulated$Normal_heat_stress$locusName,
      High_HS_down = DEGs_NormalHighHeatStress_annot$dowregulated$High_heat_stress$locusName
  ), term2gene_GO = term2gene_GO_DB,term2name_GO = term2name_GO_DB
  )

#### Summarize enrichment table, get only the top 20 enriched terms
ORA_NormalHS <- CombineUpDownRegulation_SelectTopEnrichmentTermsByCounts_Plot(
      TopCounts=20,
      ORA_NormalHigh_HS$Normal_HS_up@result,
      ORA_NormalHigh_HS$Normal_HS_down@result)

#### Summarize enrichment table, get only the top 20 enriched terms
ORA_HighHS <- CombineUpDownRegulation_SelectTopEnrichmentTermsByCounts_Plot(
  TopCounts=20,
  ORA_NormalHigh_HS$High_HS_up@result,
  ORA_NormalHigh_HS$High_HS_down@result)

ggplot(ORA_HighHS,
       aes(x=GeneRatioNumeric,y=Description, fill = p.adjust)) + 
  geom_bar(stat = "identity") + 
  #scale_fill_manual(values=GO_labelColors)+
    scale_fill_continuous(low="red", high="blue")+
  geom_text(aes(label = Count), vjust = 0.5, hjust= 0.001,colour = "black", size =4.5, fontface=2) +
  theme_bw()  + facet_grid(~Regulation) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size=12,face="bold"),
        axis.text.y= element_text(size=15,face="bold"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_text(size=15,face = "bold"),
        legend.text  = element_text(size=15,face = "bold"),
        strip.text = element_text(size=15,face = "bold"))

### Do the venn comparison between between Annabelle and Camel in response to heatsetress (H vs C)
ggvenn(
  list(
    Annabelle = dds_list_res_sig_Up$AA_HvsC_Normal$locusName,
    Camel = dds_list_res_sig_Up$CA_HvsC_Normal$locusName),
  fill_color = c("lightblue", "#B5251C"),
  set_name_size = 9,
  text_size = 8,
  fill_alpha = 0.8,
  stroke_alpha = 1,
) + theme(
  text = element_text(face="bold")
)

ggvenn(
  list(
    Annabelle = dds_list_res_sig_Down$AA_HvsC_Normal$locusName,
    Camel = dds_list_res_sig_Down$CA_HvsC_Normal$locusName),
  fill_color = c("lightblue", "#B5251C"),
  set_name_size = 9,
  text_size = 8,
  fill_alpha = 0.8,
  stroke_alpha = 1,
) + theme(
  text = element_text(face="bold")
)

ggvenn(
  list(
    Annabelle = dds_list_res_sig_Up$AA_HvsC_Normal$locusName,
    Camel = dds_list_res_sig_Down$CA_HvsC_Normal$locusName),
  fill_color = c("lightblue", "#B5251C"),
  set_name_size = 9,
  text_size = 8,
  fill_alpha = 0.8,
  stroke_alpha = 1,
) + theme(
  text = element_text(face="bold")
)

ggvenn(
  list(
    Annabelle = dds_list_res_sig_Down$AA_HvsC_Normal$locusName,
    Camel = dds_list_res_sig_Up$CA_HvsC_Normal$locusName),
  fill_color = c("lightblue", "#B5251C"),
  set_name_size = 9,
  text_size = 8,
  fill_alpha = 0.8,
  stroke_alpha = 1,
) + theme(
  text = element_text(face="bold")
)



ggvenn(
  list(
    AA_upregulated = dds_list_res_sig_Up$AA_HvsC_Normal$locusName,
    CA_upregulated = dds_list_res_sig_Up$CA_HvsC_Normal$locusName,
    
    AA_downregulated = dds_list_res_sig_Down$AA_HvsC_Normal$locusName,
    CA_downregulated = dds_list_res_sig_Down$CA_HvsC_Normal$locusName),
  fill_color = c("lightblue", "#B5251C"),
  set_name_size = 9,
  text_size = 8,
  fill_alpha = 0.8,
  stroke_alpha = 1,
) + theme(
  text = element_text(face="bold")
)


UpSetR::upset(
  fromList(
    list(
      AA_upregulated = dds_list_res_sig_Up$AA_HvsC_Normal$locusName,
      CA_upregulated = dds_list_res_sig_Up$CA_HvsC_Normal$locusName,
      
      AA_downregulated = dds_list_res_sig_Down$AA_HvsC_Normal$locusName,
      CA_downregulated = dds_list_res_sig_Down$CA_HvsC_Normal$locusName
    )
  ),text.scale = c(3, 3, 3, 1.5, 3, 5),
  sets = c("AA_upregulated",
           "AA_downregulated",
           "CA_upregulated",
           "CA_downregulated"),
  keep.order = TRUE
)

# Convert to a wide binary format
 NormalStress_DEG_DF<- stack( list(
  AA_upregulated = dds_list_res_sig_Up$AA_HvsC_Normal$locusName,
  CA_upregulated = dds_list_res_sig_Up$CA_HvsC_Normal$locusName,
  
  AA_downregulated = dds_list_res_sig_Down$AA_HvsC_Normal$locusName,
  CA_downregulated = dds_list_res_sig_Down$CA_HvsC_Normal$locusName
)) %>%
  rename(Gene = values, Category = ind) %>%
   mutate(Present = 1) %>%
   tidyr::pivot_wider(names_from = Category, values_from = Present, values_fill = 0)



UpSetR::upset(upset_df %>% column_to_rownames("Gene"))
####Get the list
NormalStress_DEG_Intersect_list<- list(
    AA_CA_Upregulated = micro.gen.extra::get_intersect_members(NormalStress_DEG_DF,groups = c("AA_upregulated","CA_upregulated")),
    AA_CA_Downregulated = micro.gen.extra::get_intersect_members(NormalStress_DEG_DF,groups = c("AA_downregulated","CA_downregulated")),
    AA_upregulated = micro.gen.extra::get_intersect_members(NormalStress_DEG_DF,groups = c("AA_upregulated")),
    CA_upregulated = micro.gen.extra::get_intersect_members(NormalStress_DEG_DF,groups = c("CA_upregulated")),
    
    AA_downregulated = micro.gen.extra::get_intersect_members(NormalStress_DEG_DF,groups = c("AA_downregulated")),
    CA_downregulated = micro.gen.extra::get_intersect_members(NormalStress_DEG_DF,groups = c("CA_downregulated")))

##run enrichment for each of the intersect....
NormalStress_DEG_Intersect_GO<- GO_KO_enrichment_with_enricher_CP(NormalStress_DEG_Intersect_list,
                                  term2gene_GO = term2gene_GO_DB,term2name_GO = term2name_GO_DB)

##Plot, dot plot with compare clust cluster
NormalStress_DEG_Intersect_GO_compareCluster <- compareCluster(
  NormalStress_DEG_Intersect_list, 
  fun = "enricher",
  pvalueCutoff = 0.05,
  pAdjustMethod = "fdr",
  TERM2GENE = term2gene_GO_DB,
  TERM2NAME = term2name_GO_DB
  
)

dotplot(NormalStress_DEG_Intersect_GO_compareCluster,
        color= "p.adjust",
        showCategory = 5) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        legend.text = element_text(size =14, face ="bold"),
        legend.title = element_text(size =12, face ="bold"))


###Get the functional annotation of the GO Terms
NormalStress_DEG_Intersect_GO_FuncAnnot <- list()
for (Intersections in names(NormalStress_DEG_Intersect_GO)) {
  
  #Separate table of both GO and KO Enrichment results and store to list
  GO_table <- NormalStress_DEG_Intersect_GO[[Intersections]]@result
  
  #rename the core_enrichment genes..
  colnames(GO_table)[colnames(GO_table) == "geneID"] <- "locusName"
  
  
  #Create the enrichment term table (summary)
  GO_table_summary<- GO_table %>% 
    dplyr::select(locusName,ID,Description,pvalue,p.adjust) %>% dplyr::filter(p.adjust < 0.05) %>%
    separate_rows(locusName, sep = "/") %>% left_join(S.tuberosum.v6.1_latest %>% 
                                                        dplyr::select(locusName,v6.1.Description,aktualisierte.annotation..bitte.ergänzen.,
                                                                      Mapman.description), by ="locusName")
  
  
  ##save everything to lsit
  NormalStress_DEG_Intersect_GO_FuncAnnot[[Intersections]] <- GO_table_summary
  
}
#####Ensure that each gene is unique...... and combine all GO terms and respective enrichment scores,pvalues etc.. of each gene together
NormalStress_DEG_Intersect_GO_FuncAnnot_unique <- list()
for (Intersections in names(NormalStress_DEG_Intersect_GO_FuncAnnot)) {
  
  #get the comparison name with both Pos and Neg table
  ORA_comparison <- NormalStress_DEG_Intersect_GO_FuncAnnot[[Intersections]]
  
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
  NormalStress_DEG_Intersect_GO_FuncAnnot_unique[[Intersections]] <- ORA_comparison_unique
  
}

library(openxlsx)
save_dfs_to_excel(df1 = NormalStress_DEG_Intersect_GO_FuncAnnot_unique$AA_CA_Upregulated, SheetName1 = "AA_CA_Upregulated",
                  df2 = NormalStress_DEG_Intersect_GO_FuncAnnot_unique$AA_CA_Downregulated, SheetName2 = "AA_CA_Downregulated",
                  df3 = NormalStress_DEG_Intersect_GO_FuncAnnot_unique$AA_upregulated, SheetName3  = "AA_upregulated",
                  df4 = NormalStress_DEG_Intersect_GO_FuncAnnot_unique$CA_upregulated, SheetName4 = "CA_upregulated",
                  df5 = NormalStress_DEG_Intersect_GO_FuncAnnot_unique$AA_downregulated,SheetName5 = "AA_downregulated",
                  df6 = NormalStress_DEG_Intersect_GO_FuncAnnot_unique$CA_downregulated, SheetName6 = "CA_downregulated",
                  file_name = "/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/06_ANALYSIS/ORA_IntersectionsNormalHeatStress.xlsx")


length(dds_list_res_sig_L2FC$AA_HvsC_Normal > 1 )
length(dds_list_res_sig_L2FC$AA_HvsC_Normal < -1 )

##Run GSEA
GSEA_Combined <- compareCluster(
    dds_list_res_sig_L2FC, 
  fun = "GSEA",
  nPermSimple = 10000,
  pvalueCutoff = 0.05,
  pAdjustMethod = "fdr",
  TERM2GENE = term2gene_GO_DB,
  TERM2NAME = term2name_GO_DB
  
)

#Count total core genes in each identified term
GSEA_Combined@compareClusterResult$Gene_count <- sapply(strsplit(as.character(GSEA_Combined@compareClusterResult$core_enrichment), "/"), length)

##PLOT terms, top 10 by counts
dotplot(GSEA_Combined,
        color= "NES",
        showCategory = 10,
        by = "Gene_count",
        size = "Gene_count") +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", 
    midpoint = 0, name = "NES") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        legend.text = element_text(size =14, face ="bold"),
        legend.title = element_text(size =12, face ="bold"))


# Extract the terms for each comparison (assuming 'compareClusterResult' contains the results for all comparisons)
 unique(GSEA_Combined@compareClusterResult$ID[which(GSEA_Combined@compareClusterResult$Cluster == "AA_HvsC_High")])
 unique(GSEA_Combined@compareClusterResult$ID[which(GSEA_Combined@compareClusterResult$Cluster == "CA_HvsC_High")])
 unique(GSEA_Combined@compareClusterResult$ID[which(GSEA_Combined@compareClusterResult$Cluster == "AA_HvsC_Normal")])
 unique(GSEA_Combined@compareClusterResult$ID[which(GSEA_Combined@compareClusterResult$Cluster == "CA_HvsC_Normal")])

# Combine the terms from all comparisons into a list
all_comparisons <- list(AA_HvsC_High = unique(GSEA_Combined@compareClusterResult$ID[which(GSEA_Combined@compareClusterResult$Cluster == "AA_HvsC_High")]), 
                        CA_HvsC_High =  unique(GSEA_Combined@compareClusterResult$ID[which(GSEA_Combined@compareClusterResult$Cluster == "CA_HvsC_High")]), 
                        AA_HvsC_Normal =  unique(GSEA_Combined@compareClusterResult$ID[which(GSEA_Combined@compareClusterResult$Cluster == "AA_HvsC_Normal")]), 
                        CA_HvsC_Normal =  unique(GSEA_Combined@compareClusterResult$ID[which(GSEA_Combined@compareClusterResult$Cluster == "CA_HvsC_Normal")]))

ggvenn::ggvenn(
  all_comparisons,
  set_name_size = 5,
  text_size = 8,
  fill_alpha = 0.6,
  stroke_alpha = 1,
) + theme(
  text = element_text(face="bold")
)

?ggvenn::geom_venn

ggplot() + 
  geom_venn(all_comparisons,
            show_percentage = TRUE) +
  theme_void()


common_terms_all <- Reduce(intersect, all_comparisons)
GSEA_Combined_common <- GSEA_Combined

GSEA_Combined_common@compareClusterResult<- GSEA_Combined_common@compareClusterResult[GSEA_Combined_common@compareClusterResult$ID %in% common_terms_all,]


##PLOT terms, only non unqiue terms
dotplot(GSEA_Combined_common,
        color= "NES",
        showCategory = 10,
        by = "Gene_count",
        size = "Gene_count") +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", 
    midpoint = 0, name = "NES") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "bold"),
        axis.text.y = element_text(size = 17, face = "bold"),
        legend.text = element_text(size =17, face ="bold"),
        legend.title = element_text(size =12, face ="bold"))



dds_list_res_GSEA <- lapply(dds_list_res_sig_L2FC, RunGSEA,
                            term2gene = term2gene_GO_DB,
                            term2name = term2name_GO_DB)



##plot the ridge plot
enrichplot::ridgeplot(dds_list_res_GSEA$CA_HvsC_High,
                      showCategory = c("protein folding",
                                       "response to hydrogen peroxide",
                                       "response to heat",
                                       "protein complex oligomerization",
                                       "hormone-mediated signaling pathway",
                                       "pectin catabolic process"),
                      label_format = 100) +  theme(axis.text.x = element_text(size = 14, face = "bold"),
                                                   axis.text.y = element_text(size = 17, face = "bold"),
                                                   legend.text = element_text(size =17, face ="bold"),
                                                   legend.title = element_text(size =12, face ="bold"))

gseaplot2(dds_list_res_GSEA$CA_HvsC_High,
         geneSetID="GO:0009620")

dds_list_res_GSEA$dds_AA_HvsC
dds_list_res_GSEA$dds_CA_HvsC
####### get the genes in the GSEA for the normal stress
dds_list_res_GSEA_functannot <- list()
for (comparison in names(dds_list_res_GSEA)) {
  
  #Separate table of both GO and KO Enrichment results and store to list
  GO_table <- dds_list_res_GSEA[[comparison]]@result
  
  #rename the core_enrichment genes..
  colnames(GO_table)[colnames(GO_table) == "core_enrichment"] <- "locusName"
  
  
  #Create the enrichment term table (summary)
  GO_table_summary<- GO_table %>% 
    dplyr::select(locusName,ID,Description,NES,pvalue,p.adjust) %>% 
    separate_rows(locusName, sep = "/") %>% left_join(S.tuberosum.v6.1_latest %>% 
                                                        dplyr::select(locusName,v6.1.Description,aktualisierte.annotation..bitte.ergänzen.,
                                                                      Mapman.description), by ="locusName")
  
  
  ##save everything to lsit
  dds_list_res_GSEA_functannot[[comparison]] <- GO_table_summary
  
}



#Ensure that each gene is unique...... and combine all GO terms and respective enrichment scores,pvalues etc.. of each gene together
dds_list_res_GSEA_functannot_unique <- list()
for (Comparison in names(dds_list_res_GSEA_functannot)) {
  
  #get the comparison name with both Pos and Neg table
  GSEA_comparison <- dds_list_res_GSEA_functannot[[Comparison]]
  
  GSEA_comparison_unique <- GSEA_comparison %>%
    group_by(locusName, v6.1.Description,aktualisierte.annotation..bitte.ergänzen.,Mapman.description) %>%
    summarize(
      ID = paste(ID, collapse = "/"),
      Description = paste(Description, collapse = "/"),
      NES = paste(NES, collapse = "/"),
      pvalue = paste(pvalue, collapse = "/"),
      p.adjust = paste(p.adjust, collapse = "/"))
  
  #save to the master list
  dds_list_res_GSEA_functannot_unique[[Comparison]] <- GSEA_comparison_unique
  
}



grid.arrange(
  ggvenn::ggvenn(
    list(
      AA_HvsC_Normal = dds_list_res_GSEA_functannot_unique$AA_HvsC_Normal$locusName,
      CA_HvsC_Normal = dds_list_res_GSEA_functannot_unique$CA_HvsC_Normal$locusName
    ), set_name_size = 6,
    text_size = 6,
    fill_alpha = 0.6
  ),
  ggvenn::ggvenn(
    list(
      AA_HvsC_High = dds_list_res_GSEA_functannot_unique$AA_HvsC_High$locusName,
      CA_HvsC_High = dds_list_res_GSEA_functannot_unique$CA_HvsC_High$locusName
    ), set_name_size = 6,
    text_size = 6,
    fill_alpha = 0.6
  )
)


ggvenn::ggvenn(
  list(
    AA_HvsC_Normal = dds_list_res_GSEA_functannot_unique$AA_HvsC_Normal$locusName,
    CA_HvsC_Normal = dds_list_res_GSEA_functannot_unique$CA_HvsC_Normal$locusName,
    AA_HvsC_High = dds_list_res_GSEA_functannot_unique$AA_HvsC_High$locusName,
    CA_HvsC_High = dds_list_res_GSEA_functannot_unique$CA_HvsC_High$locusName
  ), set_name_size = 6,
  text_size = 6,
  fill_alpha = 0.6
)

#Annabelle
dds_list_res_GSEA$AA_HvsC_Normal@result$Gene_count <- sapply(strsplit(as.character(dds_list_res_GSEA$AA_HvsC_Normal@result$core_enrichment), "/"), length)
#Canem
dds_list_res_GSEA$CA_HvsC_Normal@result$Gene_count <- sapply(strsplit(as.character(dds_list_res_GSEA$CA_HvsC_Normal@result$core_enrichment), "/"), length)

dotplot(dds_list_res_GSEA$AA_HvsC_Normal,
        color= "NES",
        showCategory = 10)

dotplot(dds_list_res_GSEA$AA_HvsC_Normal,
        color= "NES",
        showCategory = 15) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", 
    midpoint = 0, name = "NES") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "bold"),
        axis.text.y = element_text(size = 17, face = "bold"),
        legend.text = element_text(size =17, face ="bold"),
        legend.title = element_text(size =12, face ="bold"))

SelectTopTermsGSEA<- function(GSEA_list,TopCounts=10) {
  
  dds_list_res_GSEA_TopTermsCounts <- list()
  for (GSEA_Comparison in names(GSEA_list)) {
    
    ##Select the GSEA table
    GSEA_table <- GSEA_list[[GSEA_Comparison]]
    
    ##Count the core enrichment Genes
    GSEA_table@result$Gene_count <- sapply(strsplit(as.character(GSEA_table@result$core_enrichment), "/"), length)
    
    ##Add a column named regulation based on the NES score
    
    GSEA_table@result <- GSEA_table@result %>% dplyr::mutate(
      Regulation=case_when(NES > 0 ~ "Upregulated",
                           NES < 0  ~ "Downregulated"))
    
    ##Separate the upregulated from downregulated terms
    
    ##get only the top terms upregulated
    Upregulated_top_terms <- GSEA_table@result %>% 
      dplyr::filter(str_detect(Regulation,"Upregulated")) %>% 
      dplyr::filter(p.adjust < 0.05) %>%
      arrange(desc(Gene_count), .by_group =TRUE)  %>%
      slice(1:TopCounts)
    
    #Get only the top terms downregulated
    Downregulated_top_terms <- GSEA_table@result %>% 
      dplyr::filter(str_detect(Regulation,"Downregulated")) %>% 
      dplyr::filter(p.adjust < 0.05) %>%
      arrange(desc(Gene_count), .by_group =TRUE)  %>%
      slice(1:TopCounts)
    
    rbind(Upregulated_top_terms,
          Downregulated_top_terms)
    
    
    dds_list_res_GSEA_TopTermsCounts[[GSEA_Comparison]] <- rbind(Upregulated_top_terms,
                                                       Downregulated_top_terms)
  }
  
  return(dds_list_res_GSEA_TopTermsCounts)
  
}

dds_list_res_GSEA_topTermsCounts <- SelectTopTermsGSEA(dds_list_res_GSEA)

ggplot(dds_list_res_GSEA_topTermsCounts$CA_HvsC_High,
       aes(x=Gene_count,y=Description, fill = p.adjust)) + 
  geom_bar(stat = "identity") + 
  #scale_fill_manual(values=GO_labelColors)+
  scale_fill_continuous(low="red", high="blue")+
 # geom_text(aes(label = Gene_count), vjust = 0.5, hjust= 0.001,colour = "black", size =4.5, fontface=2) +
  theme_bw()  + facet_grid(~Regulation) + 
  xlab("Gene counts") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size=12,face="bold"),
        axis.text.y= element_text(size=16,face="bold"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=17,face = "bold"),
        legend.title = element_text(size=17,face = "bold"),
        legend.text  = element_text(size=17,face = "bold"),
        strip.text = element_text(size=17,face = "bold"))

View(test$AA_HvsC_Normal)

test$AA_HvsC_Normal@result$Gene_count <- sapply(strsplit(as.character(test$AA_HvsC_Normal@result$core_enrichment), "/"), length)



library(UpSetR)

# Get unique gene names
AllGeneComparisons <- list(
  AA_HvsC_Normal = dds_list_res_GSEA_functannot_unique$AA_HvsC_Normal %>% ungroup() %>% dplyr:::select(locusName, NES),
  CA_HvsC_Normal = dds_list_res_GSEA_functannot_unique$CA_HvsC_Normal %>% ungroup() %>% dplyr::select(locusName, NES),
  AA_HvsC_High = dds_list_res_GSEA_functannot_unique$AA_HvsC_High %>% ungroup() %>% dplyr::select(locusName, NES),
  CA_HvsC_High = dds_list_res_GSEA_functannot_unique$CA_HvsC_High %>% ungroup() %>% dplyr::select(locusName, NES)
)



# Get unique gene names across all comparisons
all_genes <- unique(unlist(lapply(AllGeneComparisons, function(df) df$locusName)))

# Create a dataframe to store NES values
nes_df <- data.frame(Gene = all_genes)

# Populate the dataframe with NES values or NA
for (comparison in names(AllGeneComparisons)) {
  # Create a named vector of NES values
  nes_vector <- setNames(AllGeneComparisons[[comparison]]$NES, 
                         AllGeneComparisons[[comparison]]$locusName)
  
  # Match genes and fill NES values or NA
  nes_df[[comparison]] <- nes_vector[nes_df$Gene]
}

# Set Gene names as rownames and remove the Gene column
rownames(nes_df) <- nes_df$Gene
nes_df$Gene <- NULL

# View the resulting dataframe
head(nes_df)

# Convert NES values to binary presence/absence for UpSetR
AllGeneComparisons_df <- as.data.frame(ifelse(!is.na(nes_df), 1, 0))
AllGeneComparisons_df<- AllGeneComparisons_df %>% rownames_to_column("Gene")

UpSetR::upset(
  AllGeneComparisons_df,
  text.scale = c(2, 2, 2, 2, 2, 2)
)

###Get the selected comparisons
 GeneListofselectedGeneComparisons<- list(
   AA_HvsC_Normal = as.data.frame(micro.gen.extra::get_intersect_members(AllGeneComparisons_df, "AA_HvsC_Normal")),
   CA_HvsC_Normal = as.data.frame(micro.gen.extra::get_intersect_members(AllGeneComparisons_df, "CA_HvsC_Normal")),
   AA_HvsC_High = as.data.frame(micro.gen.extra::get_intersect_members(AllGeneComparisons_df, "AA_HvsC_High")),
   CA_HvsC_High = as.data.frame(micro.gen.extra::get_intersect_members(AllGeneComparisons_df, "CA_HvsC_High")),
 IntersectAll = as.data.frame(micro.gen.extra::get_intersect_members(AllGeneComparisons_df, c("AA_HvsC_Normal","CA_HvsC_Normal","AA_HvsC_High","CA_HvsC_High"))))

nes_df <- nes_df %>%rownames_to_column("locusName")


###change the column names to ""locusName
for (comparisons in names(GeneListofselectedGeneComparisons)) {
  
  colnames(GeneListofselectedGeneComparisons[[comparisons]]) <- "locusName"
  
}


test <- dds_list_res_GSEA_functannot_unique$AA_HvsC_Normal %>% ungroup() %>% dplyr::select(locusName,Description)

test2<- GeneListofselectedGeneComparisons$IntersectAll %>% 
  inner_join(test,by="locusName") %>%
  inner_join(nes_df,by="locusName") 


View(dds_list_res_GSEA_functannot_unique$CA_HvsC_Normal)

GeneListofselectedGeneComparisons$IntersectAll

##Intersect for comparisons
GSEA_intersect_Genes_All_comparisons<- as.data.frame(Reduce(intersect, 
       list(
         AA_HvsC_Normal = dds_list_res_GSEA_functannot_unique$AA_HvsC_Normal$locusName,
         CA_HvsC_Normal = dds_list_res_GSEA_functannot_unique$CA_HvsC_Normal$locusName,
         AA_HvsC_High = dds_list_res_GSEA_functannot_unique$AA_HvsC_High$locusName,
         CA_HvsC_High = dds_list_res_GSEA_functannot_unique$CA_HvsC_High$locusName
       )))


dfs <- list(
  dds_list_res_GSEA_functannot_unique$AA_HvsC_Normal %>% ungroup() %>% select(locusName, Description),
  dds_list_res_GSEA_functannot_unique$CA_HvsC_Normal %>% ungroup() %>% select(locusName, Description),
  dds_list_res_GSEA_functannot_unique$AA_HvsC_High %>% ungroup() %>% select(locusName, Description),
  dds_list_res_GSEA_functannot_unique$CA_HvsC_High %>% ungroup() %>% select(locusName, Description)
)

# Reduce: Apply inner_join() across all data frames based on 'locusName'
GSEA_intersect_Genes_All_comparisons <- Reduce(function(x, y) inner_join(x, y, by = "locusName"), dfs)


inner_join(dds_list_res_GSEA_functannot_unique$dds_AA_HvsC,
           dds_normalstress_res_GSEA_functannot_unique$dds_AA_HvsC)

?join_by


nrow(inner_join(dds_list_res_GSEA_functannot_unique$dds_AA_HvsC,
                dds_normalstress_res_GSEA_functannot_unique$dds_AA_HvsC))





##### account for the fact that samples are paired ######

MetaData_sample<- MetaData %>% mutate(
  sample= case_when(str_detect(SampleID,"15")~ "15",
                    str_detect(SampleID,"11")~ "11",
                    str_detect(SampleID,"2")~ "2",
                    str_detect(SampleID,"3")~ "3",
                    str_detect(SampleID,"4")~ "4",
                    str_detect(SampleID,"5")~ "5")
)


#Conditions comparison
dds_list_sample<- list(
  dds_AA_HvsC= DESeqDataSetFromMatrix(countData=UniqueCountsDuplicatedHighStress_counts %>% as.data.frame() %>%
                                        dplyr::select(contains("15"),contains("4"),contains("5")) %>% as.matrix(), 
                                      colData=MetaData_sample %>% dplyr::filter(str_detect(Genotype,"AA")), 
                                      design= ~ sample + Conditions),
  
  dds_CA_HvsC= DESeqDataSetFromMatrix(countData=UniqueCountsDuplicatedHighStress_counts %>% as.data.frame() %>%
                                        dplyr::select(contains("11"),contains("2"),contains("3")) %>% as.matrix(), 
                                      colData=MetaData_sample %>% dplyr::filter(str_detect(Genotype,"CA")), 
                                      design= ~ sample + Conditions)
  
  
)


dds_list_res_SAMPLE <- list()
for (dds_name in names(dds_list_sample)) {
  
  #run the DEseq function
  dds <- DESeq(dds_list_sample[[dds_name]])
  
  design_formula <- as.character(dds@design)
  
  #get the results
  res <- results(dds, contrast = c("Conditions","Heat","Control"), pAdjustMethod = "fdr")
    

  
  rownames(res) <- gsub(".v6.1","",rownames(res))
  
  ##label the genes
  dds_list_res_SAMPLE[[dds_name]] <- res %>% 
    as.data.frame() %>% 
    rownames_to_column("locusName") %>%
    inner_join(GeneID2GeneLabel)
  
}

dds_list_res_sig_SAMPLE <- list()
for (dds_name in names(dds_list_res_SAMPLE)) {
  
  dds_list_res_sig_SAMPLE[[dds_name]] <- dds_list_res_SAMPLE[[dds_name]] %>% dplyr::filter(padj <= 0.05)
  
}




dds_list_res_sig_SAMPLE$dds_AA_HvsC %>% dplyr::filter(log2FoldChange >= 1) %>% nrow()
dds_list_res_sig_SAMPLE$dds_CA_HvsC %>% dplyr::filter(log2FoldChange >= 1) %>% nrow()

dds_list_res_sig_SAMPLE$dds_AA_HvsC %>% dplyr::filter(log2FoldChange <= -1) %>% nrow()
dds_list_res_sig_SAMPLE$dds_CA_HvsC %>% dplyr::filter(log2FoldChange <= -1) %>% nrow()


ggvenn(list(
  Corrected_sample = dds_list_res_sig_SAMPLE$dds_AA_HvsC$locusName,
norm = dds_list_res_sig$AA_HvsC_High$locusName
))


ggvenn(list(
  Corrected_sample = dds_list_res_sig_SAMPLE$dds_CA_HvsC$locusName,
  norm = dds_list_res_sig$CA_HvsC_High$locusName
))

#### For heat vs control use these valuess... accounts for the fact they are paired samples
write.table(dds_list_res_sig_SAMPLE$dds_AA_HvsC,
            "/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/06_ANALYSIS/EpiPotato_mRNA_HighStress/DEGs_HighHeat_CorrectedPaired/dds_AA_HvsC.txt",
            col.names = TRUE, row.names = FALSE,sep = "\t",quote = FALSE)

write.table(dds_list_res_sig_SAMPLE$dds_CA_HvsC,
            "/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/06_ANALYSIS/EpiPotato_mRNA_HighStress/DEGs_HighHeat_CorrectedPaired/dds_CA_HvsC.txt",
            col.names = TRUE, row.names = FALSE,sep = "\t",quote = FALSE)

####
library(tidyverse)

ggsave(
  "VennDiagrammHvC_HighHeat.svg",
  plot = ggvenn::ggvenn(
    list(
      Annabelle = dds_list_res_sig_SAMPLE$dds_AA_HvsC$locusName,
      Camel = dds_list_res_sig_SAMPLE$dds_CA_HvsC$locusName
    ),
    fill_color = c("orange4", "green4"),
    fill_alpha = 0.3,
    set_name_size = 11,
    text_size = 13
  ) ,
  device = "svg",
  path = "/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/06_ANALYSIS/EpiPotato_mRNA_HighStress/DEGs_HighHeat_CorrectedPaired/",
  width = 10,
  height = 10,
  dpi = 300)


#### GET DEGs L2FC >1 and <-1
dds_list_res_sig_SAMPLE_DEGs <- list(
 AA_upregulated = dds_list_res_sig_SAMPLE$dds_AA_HvsC %>% dplyr::filter(log2FoldChange >= 1) ,
 CA_upregulated = dds_list_res_sig_SAMPLE$dds_CA_HvsC %>% dplyr::filter(log2FoldChange >= 1) ,
 AA_downregulated = dds_list_res_sig_SAMPLE$dds_AA_HvsC %>% dplyr::filter(log2FoldChange <= -1),
 CA_downregulated = dds_list_res_sig_SAMPLE$dds_CA_HvsC %>% dplyr::filter(log2FoldChange <= -1)
)

library(ggvenn)

ggsave(
  "VennDiagrammHvC_HighHeat_all_compare.svg",
  plot = ggvenn::ggvenn(
    list(
      AA_upregulated = dds_list_res_sig_SAMPLE_DEGs$AA_upregulated$locusName,
      CA_upregulated = dds_list_res_sig_SAMPLE_DEGs$CA_upregulated$locusName,
      AA_downregulated = dds_list_res_sig_SAMPLE_DEGs$AA_downregulated$locusName,
      CA_downregulated = dds_list_res_sig_SAMPLE_DEGs$CA_downregulated$locusName
    ),
    fill_color = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A"),
    fill_alpha = 0.3,
    set_name_size = ,
    text_size = 8
  ) + theme(
    text = element_text(face = "bold")  # make all text bold
  ),
  device = "svg",
  path = "/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/06_ANALYSIS/EpiPotato_mRNA_HighStress/DEGs_HighHeat_CorrectedPaired/",
  width = 10,
  height = 10,
  dpi = 300)

ggvenn(list(
 AA_upregulated = dds_list_res_sig_SAMPLE_DEGs$AA_upregulated$locusName,
 CA_upregulated = dds_list_res_sig_SAMPLE_DEGs$CA_upregulated$locusName,
 AA_downregulated = dds_list_res_sig_SAMPLE_DEGs$AA_downregulated$locusName,
 CA_downregulated = dds_list_res_sig_SAMPLE_DEGs$CA_downregulated$locusName
),
fill_color = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A"))

### Create table for summary
DEG_CountsSummary <- data.frame(
  Genotype = c("Annabelle", "Camel"),
  Upregulated = c(4328, 4087),
  Downregulated = c(-3576, -3455)
)

# Convert to long format for ggplot2
DEG_CountsSummary_long <- DEG_CountsSummary %>%
  pivot_longer(cols = c("Upregulated", "Downregulated"),
               names_to = "Regulation", values_to = "Count")

ggplot(DEG_CountsSummary_long, aes(x = Genotype, y = Count, fill = Regulation))+
  geom_bar(stat = "identity", position = "stack", width = 0.5) + 
  scale_fill_manual(values = c("Upregulated" = "grey50", "Downregulated" = "black")) +
  geom_hline(yintercept = 0, color = "black") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        legend.text = element_text(size =20, face ="bold"),
        legend.title = element_text(size =20, face ="bold"))
  

ggsave(
  "DEG_Counts.svg",
  plot = ggplot(DEG_CountsSummary_long, aes(x = Genotype, y = Count, fill = Regulation))+
    geom_bar(stat = "identity", position = "stack", width = 0.5) + 
    scale_fill_manual(values = c("Upregulated" = "grey50", "Downregulated" = "black")) +
    geom_hline(yintercept = 0, color = "black") +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20, face = "bold"),
          axis.text.y = element_text(size = 20, face = "bold"),
          axis.title.y = element_text(size = 20, face = "bold"),
          legend.text = element_text(size =20, face ="bold"),
          legend.title = element_text(size =20, face ="bold")),
  device = "svg",
  path = "/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/06_ANALYSIS/EpiPotato_mRNA_HighStress/DEGs_HighHeat_CorrectedPaired",
  width = 6,
  height = 8,
  dpi = 300
)



nrow(dds_list_res_sig_SAMPLE_DEGs$AA_upregulated)
nrow(dds_list_res_sig_SAMPLE_DEGs$CA_upregulated)
nrow(dds_list_res_sig_SAMPLE_DEGs$AA_downregulated)
nrow(dds_list_res_sig_SAMPLE_DEGs$CA_downregulated)


nrow(dds_list_res_sig_SAMPLE_DEGs$AA_upregulated)
nrow(dds_list_res_sig_SAMPLE_DEGs$CA_upregulated)
nrow(dds_list_res_sig_SAMPLE_DEGs$AA_downregulated)
nrow(dds_list_res_sig_SAMPLE_DEGs$CA_downregulated)

library(ggvenn)
library(UpSetR)

# Convert to a wide binary format
dds_list_res_sig_SAMPLE_DEGs_DF<- stack( list(
  AA_upregulated = dds_list_res_sig_SAMPLE_DEGs$AA_upregulated$locusName,
  CA_upregulated = dds_list_res_sig_SAMPLE_DEGs$CA_upregulated$locusName,
  AA_downregulated = dds_list_res_sig_SAMPLE_DEGs$AA_downregulated$locusName,
  CA_downregulated = dds_list_res_sig_SAMPLE_DEGs$CA_downregulated$locusName
)) %>%
  rename(Gene = values, Category = ind) %>%
  mutate(Present = 1) %>%
  tidyr::pivot_wider(names_from = Category, values_from = Present, values_fill = 0)


dds_list_res_sig_SAMPLE_DEGs_intersectList<- list(
  AA_CA_Upregulated = micro.gen.extra::get_intersect_members(dds_list_res_sig_SAMPLE_DEGs_DF,groups = c("AA_upregulated","CA_upregulated")),
  AA_CA_Downregulated = micro.gen.extra::get_intersect_members(dds_list_res_sig_SAMPLE_DEGs_DF,groups = c("AA_downregulated","CA_downregulated")),
  AA_upregulated = micro.gen.extra::get_intersect_members(dds_list_res_sig_SAMPLE_DEGs_DF,groups = c("AA_upregulated")),
  CA_upregulated = micro.gen.extra::get_intersect_members(dds_list_res_sig_SAMPLE_DEGs_DF,groups = c("CA_upregulated")),
  
  AA_downregulated = micro.gen.extra::get_intersect_members(dds_list_res_sig_SAMPLE_DEGs_DF,groups = c("AA_downregulated")),
  CA_downregulated = micro.gen.extra::get_intersect_members(dds_list_res_sig_SAMPLE_DEGs_DF,groups = c("CA_downregulated")))

length(dds_list_res_sig_SAMPLE_DEGs_intersectList$AA_CA_Upregulated)
length(dds_list_res_sig_SAMPLE_DEGs_intersectList$AA_CA_Downregulated)
length(dds_list_res_sig_SAMPLE_DEGs_intersectList$AA_upregulated)
length(dds_list_res_sig_SAMPLE_DEGs_intersectList$CA_upregulated)
length(dds_list_res_sig_SAMPLE_DEGs_intersectList$AA_downregulated)
length(dds_list_res_sig_SAMPLE_DEGs_intersectList$CA_downregulated)

##run enrichment for each of the intersect....
library(clusterProfiler)
CombinedEnrichmentSets_dds_list_res_sig_SAMPLE_DEGs_intersectList<- compareCluster(
  dds_list_res_sig_SAMPLE_DEGs_intersectList, 
  fun = "enricher",
  pvalueCutoff = 0.05,
  pAdjustMethod = "fdr",
  TERM2GENE = term2gene_GO_DB,
  TERM2NAME = term2name_GO_DB
  
)

col_fun <- colorRamp2(c(-2, 0, 2), c("#2A78C6","#C83B38"))
CombinedEnrichmentSetsHighHeat_HvC_AA_CA@compareClusterResult
##PLOT terms, top 10 by counts

library(clusterProfiler)
ggsave(
  "GO_CommonUnique_AA_CA_HvC_SAMPLE.svg",
  plot = dotplot(CombinedEnrichmentSets_dds_list_res_sig_SAMPLE_DEGs_intersectList,
                 showCategory = 10,
                 by = "Count",
                 size = "Count") +
    scale_fill_gradient2(
      low = "#C83B38", mid = "#E9DDDD", high = "#2A78C6",
      midpoint = 0.025,
      name = "FDR"   ) +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20, face = "bold"),
          axis.text.y = element_text(size = 20, face = "bold"),
          legend.text = element_text(size =20, face ="bold"),
          legend.title = element_text(size =20, face ="bold")),
  device = "svg",
  path = "/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/06_ANALYSIS/EpiPotato_mRNA_HighStress/DEGs_HighHeat_CorrectedPaired",
  width = 17,
  height = 22,
  dpi = 300
)

write.table(CombinedEnrichmentSets_dds_list_res_sig_SAMPLE_DEGs_intersectList@compareClusterResult,
            "/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/06_ANALYSIS/EpiPotato_mRNA_HighStress/DEGs_HighHeat_CorrectedPaired/EnrichmentResultsCluster.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

### get hsps, only GO terms with response to heat, so unbiased
term2gene_GO_DB_HS <- term2gene_GO_DB %>% dplyr::filter(GO=="GO:0009408")



inner_join(
  
)

dds_combined_samples<- inner_join(
dds_list_res_SAMPLE$dds_AA_HvsC %>% dplyr::select(locusName,log2FoldChange, padj, v6.1.Description) %>% dplyr::rename(AA_L2FC=log2FoldChange, AA_padj=padj),
dds_list_res_SAMPLE$dds_CA_HvsC %>% dplyr::select(locusName,log2FoldChange, padj) %>% dplyr::rename(CAL2FC=log2FoldChange, CA_padj=padj),
by="locusName")

dds_combined_samples_HSPs <- dds_combined_samples %>% inner_join(
  term2gene_GO_DB_HS,
  by="locusName"
)

dds_combined_samples_HSPs

## check GO terms of all DMR DEGs
test <- c("Soltu.DM.10G020640",
  "Soltu.DM.07G021500",
  "Soltu.DM.12G001080",
  "Soltu.DM.02G017060",
  "Soltu.DM.12G001070",
  "Soltu.DM.02G012310",
  "Soltu.DM.09G026950",
  "Soltu.DM.01G004990",
  "Soltu.DM.03G013950",
  "Soltu.DM.03G016990",
  "Soltu.DM.04G035950",
  "Soltu.DM.07G025250",
  "Soltu.DM.08G004830",
  "Soltu.DM.09G026950",
  "Soltu.DM.10G011440",
  "Soltu.DM.10G020640",
  "Soltu.DM.03G013870",
  "Soltu.DM.04G026840",
  "Soltu.DM.06G016280",
  "Soltu.DM.07G015360",
  "Soltu.DM.07G021500",
  "Soltu.DM.08G027180",
  "Soltu.DM.10G016350",
  "Soltu.DM.12G001080",
  "Soltu.DM.01G036040",
  "Soltu.DM.03G016470",
  "Soltu.DM.02G017060",
  "Soltu.DM.12G001070",
  "Soltu.DM.04G035230",
  "Soltu.DM.03G003110",
  "Soltu.DM.03G016980",
  "Soltu.DM.03G036230",
  "Soltu.DM.07G025240",
  "Soltu.DM.05G010580",
  "Soltu.DM.06G015140")

test2<- term2gene_GO_DB[term2gene_GO_DB$locusName %in% unique(test),]

test3<- test2 %>% unique() %>%  inner_join(term2name_GO_DB, by ="GO") %>% unique()


ComplexHeatmap::Heatmap(dds_combined_samples_HSPs %>% dplyr::select(AA_L2FC, CA_padj) %>% as.matrix() %>% na.omit())


View(CombinedEnrichmentSets_dds_list_res_sig_SAMPLE_DEGs_intersectList@compareClusterResult)

AAvsCA_Control_DEGs %>% dplyr::filter(log2FoldChange >1 ) %>% nrow()
AAvsCA_Control_DEGs %>% dplyr::filter(log2FoldChange < -1 ) %>% nrow()


#### summarise all DMR DEGs
library(readxl)
path <- "/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/06_ANALYSIS/DMR_DEGs_BT/"
Upstream_DMR_DEGs <- read_excel(paste0(path,"DMR_DEGs_regions_BT_AA_comparisonGeno.xlsx"), sheet = "Upstream")
GB_DMR_DEGs <- read_excel(paste0(path,"DMR_DEGs_regions_BT_AA_comparisonGeno.xlsx"), sheet = "Genebody")
Downstream_DMR_DEGs <- read_excel(paste0(path,"DMR_DEGs_regions_BT_AA_comparisonGeno.xlsx"), sheet = "Downstream")


Upstream_DMR_DEGs <- Upstream_DMR_DEGs %>% dplyr::mutate(Region="Upstream")
GB_DMR_DEGs <- GB_DMR_DEGs %>% dplyr::mutate(Region="Genebody")
Downstream_DMR_DEGs <- Downstream_DMR_DEGs %>% dplyr::mutate(Region="Downstream")

DMR_DEG_ALL<- rbind(Upstream_DMR_DEGs,GB_DMR_DEGs,Downstream_DMR_DEGs)


  
DMR_DEG_ALL_long <- DMR_DEG_ALL %>% dplyr::select(Region,Context,Meth_status,regulation) %>% 
    group_by(Region,Context,Meth_status,regulation) %>% 
    summarise(Count = n()) %>%
    ungroup()
  

fisherExactTest <- list()
for (context in c("CG","CHG","CHH")) {
  
  context_df <- subset(DMR_DEG_ALL_long, Context == context)
  
  fisherExactTest[[context]] <- context_df %>%
    group_by(Region) %>%
    summarise(
      test = list(
        fisher.test(xtabs(Count ~ Meth_status + regulation, data = cur_data()))
      )
    )
  
}


ggplot(DMR_DEG_ALL_long, 
       aes(x = Meth_status, y = Count, fill = regulation)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = Count),
            position = position_dodge(width = 0.8),
            vjust = -0.3, size = 6,
            fontface = "bold")+
  facet_grid(Region~Context) +
  scale_y_continuous(
    limits = c(0, 500),        # change upper limit as needed
    breaks = seq(0, 500, 100)  # Y-axis tick marks every 100
  )+
  scale_fill_manual(values = c(
    "Upregulated" = "black",
    "Downregulated" = "grey40"
  )) +
  theme_bw()+
  theme(title = element_text(face = "bold",size = 25),
        legend.title = element_text(face = "bold",size = 18),
        legend.text = element_text(face = "bold",size = 16),
        axis.title = element_text(face = "bold",size = 18),
        axis.text = element_text(face = "bold",size=16),
        strip.text =element_text(face = "bold",size=16))

#### count by the DEGs side unique

# Assuming your tibble is named `dmr_data`
dmr_data_non_unique <- DMR_DEG_ALL %>%
  group_by(locusName) %>%
  filter(n() > 1) %>%
  ungroup()

dmr_data_non_unique$locusName %>% unique() %>% length()

## DEGs differentially methylated in both directions
test <- dmr_data_non_unique %>%
  dplyr::group_by(locusName) %>%
  dplyr::filter(all(c(-1, 1) %in% Methylation_direction)) %>%
  dplyr::ungroup()

dmr_data_non_unique_oneDirection <- dmr_data_non_unique %>%
  dplyr::group_by(locusName) %>%
  dplyr::filter(!all(c(-1, 1) %in% Methylation_direction)) %>%
  dplyr::ungroup()

dmr_data_non_unique_oneDirection_sum<- dmr_data_non_unique_oneDirection %>% dplyr::select(locusName,Meth_status,regulation) %>% unique() %>%
  dplyr::group_by(Meth_status,regulation) %>% dplyr::summarise(count = n()) %>% ungroup()

dmr_data_non_unique_oneDirection_sum %>% dplyr::filter(Meth_status == "Hypo" & regulation =="Downregulated") %>% nrow()

dmr_data_unique_sum <- ALL_DMR_DEGs %>%
  group_by(locusName) %>%
  filter(n() == 1) %>% dplyr::group_by(Meth_status,regulation) %>% 
  dplyr::summarise(count = n()) %>% ungroup()

##### create summary table

deg_category <- tibble(
  Category = c("Total", "DEGs with Multi DMRs", "DEGs with Single DMR"),
  Count = c(1786, 1040, 746)
)

# Table 2: Directionality
directionality <- tibble(
  Category = c("Bi-directional", "Uni-directional"),
  Count = c(541, 499)
)

ggplot(directionality, 
       aes(x = Category, y = Count)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = Count))


## Table 3 Unique uni-directional DMR DEGs
dmr_data_non_unique_oneDirection_sum
ggplot(dmr_data_non_unique_oneDirection_sum, 
       aes(x = Meth_status, y = count, fill = regulation)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = count),
            position = position_dodge(width = 0.8),
            vjust = -0.3, size = 4,
            fontface = "bold")+ theme_bw()+
  theme(title = element_text(face = "bold",size = 25),
        legend.title = element_text(face = "bold",size = 18),
        legend.text = element_text(face = "bold",size = 16),
        axis.title = element_text(face = "bold",size = 18),
        axis.text = element_text(face = "bold",size=16),
        axis.text.x = element_text(size=16),
        strip.text =element_text(face = "bold",size=16))

######## Harverst data ####
Ernte_Hochhitze <- read_excel("/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/Ernte_Hochhitze.xlsx")

Ernte_Hochhitze<- Ernte_Hochhitze[,c(1:4,7)]
colnames(Ernte_Hochhitze) <- c("Sample","Sample_ID","Stress_length","Tuber_weight","Tuber_number")

Ernte_Hochhitze_summary<- Ernte_Hochhitze %>% group_by(Sample,Stress_length) %>%
  dplyr::summarise(
    mean_Tuber_weight = mean(Tuber_weight, na.rm = TRUE),
    sd_Tuber_weight   = sd(Tuber_weight, na.rm = TRUE),
    mean_Tuber_number = mean(Tuber_number, na.rm = TRUE),
    sd_Tuber_number   = sd(Tuber_number, na.rm = TRUE)
  )

library(tidyverse)

Ernte_Hochhitze_long<- Ernte_Hochhitze %>%
  pivot_longer(cols = c(Tuber_weight, Tuber_number),
               names_to = "Variable",
               values_to = "Value")

plot1<- ggplot(Ernte_Hochhitze_long, 
       aes(x = as.factor(Stress_length), y = Value, fill = Sample)) +
  geom_boxplot(width = 0.6, position = position_dodge())  +
  facet_wrap(.~Variable, nrow =2, scales = "free_y") +
  scale_fill_manual(values=c("lightblue3","pink3"))+
  theme_bw()+
  theme(title = element_text(face = "bold",size = 25),
        legend.title = element_text(face = "bold",size = 18),
        legend.text = element_text(face = "bold",size = 16),
        axis.title = element_text(face = "bold",size = 18),
        axis.text = element_text(face = "bold",size=16),
        axis.text.x = element_text(size=16),
        strip.text =element_text(face = "bold",size=16))

ggsave(plot = plot1)
ggsave(
  "AA_CA_harvest.svg",
  plot = plot1,
  device = "svg",
  path = "/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/mRNA_HighStress/",
  width = 8,
  height = 10,
  dpi = 300
)
ggplot(Ernte_Hochhitze_summary, aes(x = as.factor(Stress_length), 
                       y = mean_Tuber_weight, 
                       fill = Sample)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(aes(ymin = mean_Tuber_weight - sd_Tuber_weight, 
                    ymax = mean_Tuber_weight + sd_Tuber_weight),
                position = position_dodge(width = 0.7), width = 0.2) +
  scale_fill_manual(values = c("lightblue3", "pink3")) +
  theme_bw()+   
  theme(title = element_text(face = "bold",size = 25),
                      legend.title = element_text(face = "bold",size = 18),
                      legend.text = element_text(face = "bold",size = 16),
                      axis.title = element_text(face = "bold",size = 18),
                      axis.text = element_text(face = "bold",size=16),
                      axis.text.x = element_text(size=16),
                      strip.text =element_text(face = "bold",size=16))

ggplot(Ernte_Hochhitze_summary, aes(x = as.factor(Stress_length), 
                                    y = mean_Tuber_number, 
                                    fill = Sample)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(aes(ymin = mean_Tuber_number - sd_Tuber_number, 
                    ymax = mean_Tuber_number + sd_Tuber_number),
                position = position_dodge(width = 0.7), width = 0.2) +
  scale_fill_manual(values = c("lightblue3", "pink3")) +
  theme_bw()+   
  theme(title = element_text(face = "bold",size = 25),
        legend.title = element_text(face = "bold",size = 18),
        legend.text = element_text(face = "bold",size = 16),
        axis.title = element_text(face = "bold",size = 18),
        axis.text = element_text(face = "bold",size=16),
        axis.text.x = element_text(size=16),
        strip.text =element_text(face = "bold",size=16))


### look at all methylase and demethylase genes..
#StMET1-4A   Soltu.DM.04G014710
#StMET2-11A  Soltu.DM.11G013230
#StCMT3-1A   Soltu.DM.01G001630
#StCMT3-8A   Soltu.DM.08G003560
#StCMT3-12A  Soltu.DM.12G000130
#StDNMT2-8A  Soltu.DM.08G015570
#StDNMT2-8B  Soltu.DM.08G015580
#StDRM1-2A   Soltu.DM.02G006560
#StDRM2-4A   Soltu.DM.04G000230
#StDRM3-10A  Soltu.DM.10G030090
#
StMet_genes <- c("Soltu.DM.04G014710",
                  "Soltu.DM.11G013230",
                  "Soltu.DM.01G001630",
                  "Soltu.DM.08G003560",
                  "Soltu.DM.12G000130",
                  "Soltu.DM.08G015570",
                  "Soltu.DM.08G015580",
                  "Soltu.DM.02G006560",
                  "Soltu.DM.04G000230",
                  "Soltu.DM.10G030090")

#DNA demethylases
#StDeMet1 Soltu.DM.11G005260
#StDeMet2 Soltu.DM.09G004240
#StDeMet3 Soltu.DM.10G024770
#StDeMet4 Soltu.DM.03G037240
#StDeMet5 Soltu.DM.04G021680
#StDeMet6 Soltu.DM.06G028750
#StDeMet7 Soltu.DM.09G025490
#StDeMet8 Soltu.DM.03G027780
StDeMets_genes <- c("Soltu.DM.11G005260",
                    "Soltu.DM.09G004240",
                    "Soltu.DM.10G024770",
                    "Soltu.DM.03G037240",
                    "Soltu.DM.04G021680",
                    "Soltu.DM.06G028750",
                    "Soltu.DM.09G025490",
                    "Soltu.DM.03G027780")


StMet_genes_L2FC<- inner_join(dds_list_res_SAMPLE$dds_AA_HvsC %>% dplyr::filter(locusName %in% StMet_genes) %>% dplyr::select(locusName,log2FoldChange,Gene_label) %>% dplyr::rename(AA=log2FoldChange),
      dds_list_res_SAMPLE$dds_CA_HvsC %>% dplyr::filter(locusName %in% StMet_genes) %>% dplyr::select(locusName,log2FoldChange,Gene_label) %>% dplyr::rename(CA=log2FoldChange),
      by=c("locusName","Gene_label"))

StMet_genes_L2FC<- StMet_genes_L2FC %>% dplyr::mutate(
  Gene_label=case_when(str_detect(locusName,"Soltu.DM.04G014710") ~ "StMET1-4A, Soltu.DM.04G014710",
                       str_detect(locusName,"Soltu.DM.11G013230") ~ "StMET2-11A, Soltu.DM.11G013230",
                       str_detect(locusName,"Soltu.DM.01G001630") ~ "StCMT3-1A, Soltu.DM.01G001630",
                       str_detect(locusName,"Soltu.DM.08G003560") ~ "StCMT3-8A, Soltu.DM.08G003560",
                       str_detect(locusName,"Soltu.DM.12G000130") ~ "StCMT3-12A, Soltu.DM.12G000130",
                       str_detect(locusName,"Soltu.DM.08G015570") ~ "StDNMT2-8A, Soltu.DM.08G015570",
                       str_detect(locusName,"Soltu.DM.08G015580") ~ "StDNMT2-8B, Soltu.DM.08G015580",
                       str_detect(locusName,"Soltu.DM.02G006560") ~ "StDRM1-2A, Soltu.DM.02G006560",
                       str_detect(locusName,"Soltu.DM.04G000230") ~ "StDRM2-4A, Soltu.DM.04G000230",
                       str_detect(locusName,"Soltu.DM.10G030090") ~ "StDRM3-10A, Soltu.DM.10G030090")
)

library(ComplexHeatmap)
col_fun <- colorRamp2(c(-2, 0, 2), c("#2A78C6","white","#C83B38"))
##StMets
Heatmap(as.matrix(StMet_genes_L2FC %>% column_to_rownames("Gene_label") %>% dplyr::select(AA,CA)),
        column_names_gp = grid::gpar(fontsize = 15),
        row_names_gp = grid::gpar(fontsize = 15),
        col = col_fun,
        cluster_columns = F,
        heatmap_legend_param = list(title = expression(bold(log["2"]~FC)),
                                    at = c(-2, -1, 0, 1, 2),   # <- Add this line to define legend ticks
                                    labels = c("-2", "-1", "0", "1", "2"),
                                    title_gp = gpar(fontsize = 15, fontface = "bold"),
                                    labels_gp = gpar(fontsize = 15, fontface = "bold"))) 


StDeMets_genes_L2FC<- inner_join(dds_list_res_SAMPLE$dds_AA_HvsC %>% dplyr::filter(locusName %in% StDeMets_genes) %>% dplyr::select(locusName,log2FoldChange,Gene_label) %>% dplyr::rename(AA=log2FoldChange),
                              dds_list_res_SAMPLE$dds_CA_HvsC %>% dplyr::filter(locusName %in% StDeMets_genes) %>% dplyr::select(locusName,log2FoldChange,Gene_label) %>% dplyr::rename(CA=log2FoldChange),
                              by=c("locusName","Gene_label"))

StDeMets_genes_L2FC<- StDeMets_genes_L2FC %>% dplyr::mutate(
  Gene_label=case_when(str_detect(locusName,"Soltu.DM.11G005260") ~ "StDeMet1, Soltu.DM.11G005260",
                       str_detect(locusName,"Soltu.DM.09G004240") ~ "StDeMet2, Soltu.DM.09G004240",
                       str_detect(locusName,"Soltu.DM.10G024770") ~ "StDeMet3, Soltu.DM.10G024770",
                       str_detect(locusName,"Soltu.DM.03G037240") ~ "StDeMet4, Soltu.DM.03G037240",
                       str_detect(locusName,"Soltu.DM.04G021680") ~ "StDeMet5, Soltu.DM.04G021680",
                       str_detect(locusName,"Soltu.DM.06G028750") ~ "StDeMet6, Soltu.DM.06G028750",
                       str_detect(locusName,"Soltu.DM.09G025490") ~ "StDeMet7, Soltu.DM.09G025490",
                       str_detect(locusName,"Soltu.DM.03G027780") ~ "StDeMet8, Soltu.DM.03G027780")
)

##ST DeMet
Heatmap(as.matrix(StDeMets_genes_L2FC %>% column_to_rownames("Gene_label") %>% dplyr::select(AA,CA)),
        column_names_gp = grid::gpar(fontsize = 15),
        row_names_gp = grid::gpar(fontsize = 15),
        col = col_fun,
        cluster_columns = F,
        heatmap_legend_param = list(title = expression(bold(log["2"]~FC)),
                                    at = c(-2, -1, 0, 1, 2),   # <- Add this line to define legend ticks
                                    labels = c("-2", "-1", "0", "1", "2"),
                                    title_gp = gpar(fontsize = 15, fontface = "bold"),
                                    labels_gp = gpar(fontsize = 15, fontface = "bold")))

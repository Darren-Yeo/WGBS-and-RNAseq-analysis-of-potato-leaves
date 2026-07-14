library(tidyverse)

data <- read.table("/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/10_TE_script/Results_AA/AA_TEtranscripts_out.cntTable",header=T,row.names=1)
min_read <- 1
data <- data[apply(data,1,function(x){max(x)}) > min_read,]

filtered <- data[!grepl("^Soltu\\.DM\\.", rownames(data)), ]

filtered_1 <- filtered[!grepl("^\\d+:\\d+:\\d+$", rownames(filtered)),]
filtered_2 <- filtered_1[!grepl("^\\([ACGT]+\\)n", rownames(filtered_1)),]


rownames(data)

clean_TE_names <- data[grepl("^\\d+:\\d+:\\d+$", rownames(data)),]
filtered <- data[grepl("^\\([ACGT]+\\)n", rownames(data)),]

filtered <- filtered[!grepl("^Soltu.DM.", filtered)]

data <- data[rownames(data) %in% filtered,]





### high heat stres
Transcripts<- list(
AA_TEtranscripts_out_gene_TE_analysis =read.delim("/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/10_TE_script/Results_AA/AA_TEtranscripts_out_gene_TE_analysis.txt"),
CA_TEtranscripts_out_gene_TE_analysis =read.delim("/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/10_TE_script/Results_CA/CA_TEtranscripts_out_gene_TE_analysis.txt"))

## combine the results 

AA_results<- AA_TEtranscripts_out_gene_TE_analysis %>% rownames_to_column("locusName") %>% select(locusName, log2FoldChange, padj) %>% rename(AA_L2FC = log2FoldChange, AA_Padj = padj)
CA_results<- CA_TEtranscripts_out_gene_TE_analysis %>% rownames_to_column("locusName") %>% select(locusName, log2FoldChange, padj) %>% rename(CA_L2FC = log2FoldChange, CA_Padj = padj)

ResultsCombined<- inner_join(AA_results,CA_results,by="locusName")




TranscriptsTE_only <- list()
for (TE in names(Transcripts)) {
  
  data <- Transcripts[[TE]]
  
  clean_TE_names <- rownames(data)[!grepl("^\\d+:\\d+:\\d+$", rownames(data))]
  filtered <- clean_TE_names[!grepl("^\\([ACGT]+\\)n", clean_TE_names)]
  
  filtered <- filtered[!grepl("^Soltu.DM.", filtered)]
  
  data <- data[rownames(data) %in% filtered,]
  
  TranscriptsTE_only[[TE]] <- data
  
}


AA_Up<- TranscriptsTE_only$AA_TEtranscripts_out_gene_TE_analysis %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 1)
AA_Down<- TranscriptsTE_only$AA_TEtranscripts_out_gene_TE_analysis %>% filter(padj <= 0.05) %>% filter(log2FoldChange <= -1)


CA_Up<- TranscriptsTE_only$CA_TEtranscripts_out_gene_TE_analysis %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 1)
CA_Down<- TranscriptsTE_only$CA_TEtranscripts_out_gene_TE_analysis %>% filter(padj <= 0.05) %>% filter(log2FoldChange <= -1)

AA_sig <- TranscriptsTE_only$AA_TEtranscripts_out_gene_TE_analysis %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 1 | log2FoldChange <= -1)
CA_sig <- TranscriptsTE_only$CA_TEtranscripts_out_gene_TE_analysis %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 1 | log2FoldChange <= -1)

CA_sig <- CA_sig %>% rownames_to_column("TE")
AA_sig <- AA_sig %>% rownames_to_column("TE")

ggvenn::ggvenn(
  list(
  Annabelle  = AA_sig$TE,
  Camel  = CA_sig$TE
  ),
  text_size = 8,
  set_name_size = 8
)

###get the l2fc

SigTEs<- unique(c(AA_sig$TE,
                CA_sig$TE))

SigTEs_sort<- data.frame(TEs = sort(SigTEs))
SigTEs_sort %>% filter(str_detect(TEs,"Gyp|GYP")) %>% nrow()
SigTEs_sort %>% filter(str_detect(TEs,"Copia|COP|SHACOP")) %>% nrow()
SigTEs_sort %>% filter(str_detect(TEs,"rnd-")) %>% nrow()
SigTEs_sort %>% filter(str_detect(TEs,"LINE|L1|RTE")) %>% nrow()

SigTEs_sort %>% filter(str_detect(TEs,"CACTA|DNA|EnSPM|hAT|HARB|Vandal|Mariner|MuDR")) %>% nrow()


writeLines(SigTEs, "/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/10_TE_script/SigTEs.txt")

ResultsCombined_sig <- ResultsCombined[ResultsCombined$locusName %in% SigTEs,]
rownames(ResultsCombined_sig) <- NULL
ResultsCombined_sig_l2fc <- ResultsCombined_sig %>% select(locusName,AA_L2FC,CA_L2FC) %>% column_to_rownames("locusName")

ComplexHeatmap::Heatmap(ResultsCombined_sig_l2fc,
                        show_row_names = FALSE)


#### Get the norm counts
AA_TEtranscripts_out <- read.delim("/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/10_TE_script/Results_AA/AA_TEtranscripts_out.cntTable")
CA_TEtranscripts_out <- read.delim("/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/10_TE_script/Results_CA/CA_TEtranscripts_out.cntTable")

transcripts_out <- inner_join(AA_TEtranscripts_out,CA_TEtranscripts_out, by="gene.TE")

colnames(transcripts_out) <- gsub("_Run2_Aligned.sortedByCoord.out.bam.T","",colnames(transcripts_out))
colnames(transcripts_out) <- gsub("_Run2_Aligned.sortedByCoord.out.bam.C","",colnames(transcripts_out))
colnames(transcripts_out) <- gsub("X.media.rna.BAIJI.Darren_Doktorand.Darren_Epipotato.mRNA_HighStress.03_STAR.","",colnames(transcripts_out))

transcripts_out <- transcripts_out %>% column_to_rownames("gene.TE")
MetaData<- data.frame(SampleID=colnames(transcripts_out))

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



library(DESeq2)
dds<- DESeqDataSetFromMatrix(
  countData = transcripts_out,
  colData = MetaData,
  design = ~Conditions +Genotype
)

dds <- DESeq(dds)

dds_vst <- assay(vst(dds))
dds_norm <- counts(dds, normalized=TRUE)

###look at overall expression of TEs
dds_vst_out <- dds_vst[!grepl("Soltu\\.DM\\.", rownames(dds_vst)), ]
#dds_vst_out <- dds_vst_out[!grepl("^\\d+:\\d+:\\d+$", rownames(dds_vst_out)),]
#dds_vst_out <- dds_vst_out[!grepl("^\\([ACGT]+\\)n", rownames(dds_vst_out)),]
colnames(dds_vst_out) <- c("AA37C1","AA37C2","AA37C3","AA1","AA2","AA3","CA37C1","CA37C2","CA37C3","CA1","CA2","CA3")

dds_vst_out_mean <- data.frame(
  AA     = rowMeans(dds_vst_out[, c("AA1", "AA2", "AA3")]),
  AA37C  = rowMeans(dds_vst_out[, c("AA37C1", "AA37C2", "AA37C3")]),
  CA     = rowMeans(dds_vst_out[, c("CA1", "CA2", "CA3")]),
  CA37C  = rowMeans(dds_vst_out[, c("CA37C1", "CA37C2", "CA37C3")])
)

dds_vst_out_mean_long<- dds_vst_out_mean %>% 
  t() %>% scale() %>% t() %>% 
  as.data.frame() %>%
  rownames_to_column("Feature") %>%
  pivot_longer(
    cols = c(AA37C, AA, CA37C, CA),
    names_to = "Group",
    values_to = "Value"
  )

ggplot() +
  geom_boxplot(data = dds_vst_out_mean_long, 
               aes(x = Group, y = Value, fill = Group),
               position = position_dodge(1), width = 0.4) +
  scale_fill_manual(values = c(
    "AA" = "lightblue3",
    "AA37C" = "blue",
    "CA" = "pink3",
    "CA37C" = "red4"
  )) +
  ylab("Expression (Z scores)") +
  xlab("Genotypes") +
  theme_bw() +
  theme(title = element_text(face = "bold",size = 25),
        legend.title = element_text(face = "bold",size = 18),
        legend.text = element_text(face = "bold",size = 16),
        legend.position = "none",
        axis.title = element_text(face = "bold",size = 18),
        axis.text = element_text(face = "bold",size=16),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,size=16),
        strip.text =element_text(face = "bold",size=16))

dds_vst_sig <- dds_vst[rownames(dds_vst) %in% SigTEs,]
dds_norm_sig <- dds_norm[rownames(dds_norm) %in% SigTEs,]


dds_vst_sig_scale <- dds_vst_sig %>% t() %>% scale() %>% t()

library(circlize)
library(grid)
colnames(dds_vst_sig_scale) <- c("AA_H1", "AA_H2", "AA_H3","AA_C1","AA_C2","AA_C3","CA_H1","CA_H2","CA_H3","CA_C1","CA_C2","CA_C3")
ComplexHeatmap::Heatmap(dds_vst_sig_scale,
                        show_row_names = FALSE,
                        column_names_gp = grid::gpar(fontsize = 15),
                        row_names_gp = grid::gpar(fontsize = 15),
                        col = colorRamp2(c(-2, 0, 2), c("blue","white","red")),
                        heatmap_legend_param = list(
                          title = "Z score",
                          title_gp = gpar(fontsize = 15),
                          labels_gp = gpar(fontsize = 15)
                        ))

col_fun <- colorRamp2(c(-2, 0, 2), c("#2A78C6","white","#C83B38"))

##convert to long
Transcripts_long <- dds_vst_sig_scale %>% as.data.frame() %>%
  mutate(gene.TE = rownames(.)) %>%
  pivot_longer(
    cols = -gene.TE,
    names_to = "SampleID",
    values_to = "Z_score"
  )

dds_norm_sig_long <- dds_norm_sig %>% as.data.frame() %>%
  mutate(gene.TE = rownames(.)) %>%
  pivot_longer(
    cols = -gene.TE,
    names_to = "SampleID",
    values_to = "norm_counts"
  )



Transcripts_long_label<- Transcripts_long%>% 
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

dds_norm_sig_long_label <-  dds_norm_sig_long%>% 
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


ggplot(dds_norm_sig_long_label %>% filter(gene.TE=="Copia-18_BD-I:Copia-18_BD-I:Copia-18_BD-I"), aes(x = Conditions, y = norm_counts, fill = Conditions)) +
  geom_boxplot() +
  facet_wrap(~ Genotype) +
  theme_minimal() +
  labs(
    title = "TE Expression by Condition and Genotype",
    x = "Condition",
    y = "Expression (normalized)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




#####
### moderate heat stress
TranscriptsModHS<- list(
  AA_TE=read.delim("/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/10_TE_script/Results_AA_moderate_HS/AA_moderate_HS_TEtranscripts_out_gene_TE_analysis.txt"),
  CA_TE=read.delim("/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/10_TE_script/Results_CA_moderate_HS/CA_moderate_HS_TEtranscripts_out_gene_TE_analysis.txt"))


TranscriptsModHS_TE_only <- list()
for (TE in names(TranscriptsModHS)) {
  
  data <- TranscriptsModHS[[TE]]
  
  clean_TE_names <- rownames(data)[!grepl("^\\d+:\\d+:\\d+$", rownames(data))]
  filtered <- clean_TE_names[!grepl("^\\([ACGT]+\\)n", clean_TE_names)]
  
  filtered <- filtered[!grepl("^Soltu.DM.", filtered)]
  
  data <- data[rownames(data) %in% filtered,]
  
  TranscriptsModHS_TE_only[[TE]] <- data
  
}

TranscriptsModHS_TE_only$AA_TE %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 1) %>% nrow()
TranscriptsModHS_TE_only$AA_TE %>% filter(padj <= 0.05) %>% filter(log2FoldChange <= -1) %>% nrow()

TranscriptsModHS_TE_only$CA_TE %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 1) %>% nrow()
TranscriptsModHS_TE_only$CA_TE %>% filter(padj <= 0.05) %>% filter(log2FoldChange <= -1) %>% nrow()

AA_sigModHS <- TranscriptsTE_only$AA_TEtranscripts_out_gene_TE_analysis %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 1 | log2FoldChange <= -1)
CA_sigModHS <- TranscriptsTE_only$CA_TEtranscripts_out_gene_TE_analysis %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 1 | log2FoldChange <= -1)

AA_sigModHS <- AA_sigModHS %>% rownames_to_column("TE")
CA_sigModHS <- CA_sigModHS %>% rownames_to_column("TE")

SigTEs_ModHS <- unique(c(AA_sig$TE,
                         CA_sig$TE))

####do pca 

#### Get the norm counts
AA_TEtranscriptsModHS_out <- read.delim("/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/10_TE_script/Results_AA_moderate_HS/AA_moderate_HS_TEtranscripts_out.cntTable")
CA_TEtranscriptsModHS_out <- read.delim("/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/10_TE_script/Results_CA_moderate_HS/CA_moderate_HS_TEtranscripts_out.cntTable")

transcripts_ModHS_out <- inner_join(AA_TEtranscriptsModHS_out,CA_TEtranscriptsModHS_out, by="gene.TE")

colnames(transcripts_ModHS_out) <- gsub("_Run2_Aligned.sortedByCoord.out.bam.T","",colnames(transcripts_ModHS_out))
colnames(transcripts_ModHS_out) <- gsub("_Run2_Aligned.sortedByCoord.out.bam.C","",colnames(transcripts_ModHS_out))
colnames(transcripts_ModHS_out) <- gsub("X.media.rna.BAIJI.Darren_Doktorand.Darren_Epipotato.04_STAR.","",colnames(transcripts_ModHS_out))


colnames(transcripts_ModHS_out)[c(4,7)] <- c("AA_9_H3","AA_2_C2") #

transcripts_ModHS_out <- transcripts_ModHS_out %>% column_to_rownames("gene.TE")

MetaData_mod<- data.frame(SampleID=colnames(transcripts_ModHS_out))

MetaData_mod<- MetaData_mod%>% 
  mutate(Conditions=case_when(str_detect(SampleID,"_C") ~ "Control",
                              str_detect(SampleID,"_H") ~ "Heat"),
         Genotype=case_when(str_detect(SampleID,"AA") ~ "AA",
                            str_detect(SampleID,"CA") ~ "CA"))

dds_mod<- DESeqDataSetFromMatrix(
  countData = transcripts_ModHS_out,
  colData = MetaData_mod,
  design = ~Conditions +Genotype
)

dds_mod <- DESeq(dds_mod)

dds_mod_vst <- assay(vst(dds_mod))

library(factoextra)
library(FactoMineR)

PCA <- PCA(t(dds_mod_vst), scale.unit = TRUE, graph = FALSE)

fviz_pca_ind(PCA, 
             addEllipses = FALSE, label = "var",
             col.var = "black", repel = TRUE) + 
  geom_point(aes(color=MetaData_mod$Conditions), size = 3) +
  geom_text(aes(label = MetaData_mod$Genotype), 
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

dds_mod_vst_sig <- dds_mod_vst[rownames(dds_mod_vst) %in% SigTEs_ModHS,]

dds_mod_vst_sig_scale <- dds_mod_vst_sig %>% t() %>% scale() %>% t()

ComplexHeatmap::Heatmap(dds_mod_vst_sig_scale %>% na.omit(),
                        show_row_names = FALSE,
                        column_names_gp = grid::gpar(fontsize = 15),
                        row_names_gp = grid::gpar(fontsize = 15),
                        col = colorRamp2(c(-2, 0, 2), c("blue","white","red")),
                        heatmap_legend_param = list(
                          title = "Z score",
                          title_gp = gpar(fontsize = 15),
                          labels_gp = gpar(fontsize = 15)
                        ))

##convert to long
Transcripts_ModHS_long <- dds_mod_vst_sig_scale %>% as.data.frame() %>%
  mutate(gene.TE = rownames(.)) %>%
  pivot_longer(
    cols = -gene.TE,
    names_to = "SampleID",
    values_to = "Z_score"
  )

##inlcude label to long table
Transcripts_ModHS_long_label<- Transcripts_ModHS_long%>% 
  mutate(Conditions=case_when(str_detect(SampleID,"_C") ~ "Control",
                              str_detect(SampleID,"_H") ~ "Heat"),
         Genotype=case_when(str_detect(SampleID,"AA") ~ "AA",
                            str_detect(SampleID,"CA") ~ "CA"))


ggplot(Transcripts_ModHS_long_label %>% filter(gene.TE=="ATLINE1_7:ATLINE1_7:ATLINE1_7"), aes(x = Conditions, y = Z_score, fill = Conditions)) +
  geom_boxplot() +
  facet_wrap(~ Genotype) +
  theme_minimal() +
  labs(
    title = "TE Expression by Condition and Genotype",
    x = "Condition",
    y = "Expression (normalized)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

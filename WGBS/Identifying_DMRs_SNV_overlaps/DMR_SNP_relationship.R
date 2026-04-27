library(tidyverse)
library(rtracklayer)
library(readxl)
library(vcfR)

vcf <- read.vcfR("/media/rna/Epipotato16TB/WGS/05_VARIANT_CALL/filtered_output_parallel_NEWQUAL30.vcf")
PrepareVCFfile_GetimportantInfo<- function(VCFR_file){
  
  #convert to tidy file
  my.vcf.tidy <- vcfR2tidy(VCFR_file)
  
  my.vcf.df<- my.vcf.tidy
  
  my.vcf.df_subset<- my.vcf.df %>% as.data.frame() %>% select(CHROM,ID,POS,QUAL,TYPE)
  
  return(my.vcf.df_subset)
  
}

vcf_edit<- PrepareVCFfile_GetimportantInfo(vcf)

vcf@gt

Enrichment_GB_candidates<- Enrichment_GB %>% dplyr::select(locusName)

DMR_DEGs_regions_BT_AA_comparisonGeno<- list(
  Upstream = read_excel("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/AAvsCA_Controls_Summary_1st2nd_EXP/DMR_DEGs_Controls_BT/DMR_DEGs_regions_BT_AA_comparisonGeno.xlsx", sheet = "Upstream"),
  Genebody = read_excel("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/AAvsCA_Controls_Summary_1st2nd_EXP/DMR_DEGs_Controls_BT/DMR_DEGs_regions_BT_AA_comparisonGeno.xlsx", sheet = "Genebody"),
  Downstream = read_excel("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/AAvsCA_Controls_Summary_1st2nd_EXP/DMR_DEGs_Controls_BT/DMR_DEGs_regions_BT_AA_comparisonGeno.xlsx", sheet = "Downstream")
)

Enrichment_GB_candidates_DF<- inner_join(DMR_DEGs_regions_BT_AA_comparisonGeno %>% dplyr::select(chr,start_DMR,end_DMR,locusName,Strand,Methylation_direction,v6.1.Description)
,Enrichment_GB_candidates,
by="locusName")

Enrichment_GB_AACompare_candidates_GR <- GRanges(
  seqnames = Enrichment_GB_candidates_DF$chr,
  ranges = IRanges(start = Enrichment_GB_candidates_DF$start_DMR, 
                   end = Enrichment_GB_candidates_DF$end_DMR),
  locusName = Enrichment_GB_candidates_DF$locusName,
  score = Enrichment_GB_candidates_DF$Methylation_direction,
  Description = Enrichment_GB_candidates_DF$v6.1.Description
)

export(Enrichment_GB_AACompare_candidates_GR, "/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/SNP_DMRs/Enrichment_GB_AACompare_candidates_GR.bed", 
       format = "bed")


### SNP ####
DMR_SNP_overlap <- read.delim("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/SNP_DMRs/DMR_SNP_overlap.tsv", header=FALSE)

###
DMR_SNP_overlap_edit<- DMR_SNP_overlap[,c(1:7,10:16)]
colnames(DMR_SNP_overlap_edit) <- c("chr","POS","V3","REF","ALT","QUAL","V7","AA","CA","Chr_dmr","start_DMR","end_DMR","V15","Meth_Direction")


extract_dosage <- function(gt_string) {
  gt <- strsplit(gt_string, ":")[[1]][1]       # Get only genotype part
  alleles <- unlist(strsplit(gt, "/"))         # Split to alleles
  sum(as.numeric(alleles))                     # Count ALT alleles
}

DMR_SNP_overlap_edit$AA_DOS <- sapply(DMR_SNP_overlap_edit$AA, extract_dosage)
DMR_SNP_overlap_edit$CA_DOS <- sapply(DMR_SNP_overlap_edit$CA, extract_dosage)


### clean table 

DMR_SNP_overlap_edit_clean<- DMR_SNP_overlap_edit %>% select(chr,POS,REF,ALT,AA_DOS,CA_DOS,start_DMR,end_DMR,Meth_Direction)

###annotate before calculate so group by genes instead

Enrichment_GB_candidates_DF_edit <- Enrichment_GB_candidates_DF %>% select(chr,start_DMR,end_DMR, locusName)

##Because of bedtools stupid 0 start 
Enrichment_GB_candidates_DF_edit$start_DMR <- Enrichment_GB_candidates_DF_edit$start_DMR-1


DMR_SNP_overlap_edit_clean <- inner_join(DMR_SNP_overlap_edit_clean,
                                         Enrichment_GB_candidates_DF_edit,
                                         by=c("chr","start_DMR","end_DMR"))


DMR_SNP_overlap_CT_GA <- DMR_SNP_overlap_edit_clean %>% filter(
  (DMR_SNP_overlap_edit_clean$REF == "C" & DMR_SNP_overlap_edit_clean$ALT =="T") |
  (DMR_SNP_overlap_edit_clean$REF == "G" & DMR_SNP_overlap_edit_clean$ALT =="A") 
)

##SUM up total Alt allele
DMR_SNP_overlap_CT_GA_SUM <- DMR_SNP_overlap_CT_GA %>% 
  group_by(chr, start_DMR, end_DMR, locusName, Meth_Direction) %>%
  summarise(
    n_snps = n(),
    sum_AA_DOS = sum(AA_DOS),
    sum_CA_DOS = sum(CA_DOS)
  ) %>%
  ungroup()

##Add a column to classify alt allele difference
DMR_SNP_overlap_CT_GA_SUM_tagged <- DMR_SNP_overlap_CT_GA_SUM %>%
  mutate(
    ALT_comparison = case_when(
      sum_AA_DOS > sum_CA_DOS ~ "AA_more_ALT",
      sum_AA_DOS < sum_CA_DOS ~ "AA_less_ALT",
      sum_AA_DOS == sum_CA_DOS ~ "Equal_ALT"
    )
  )

### finally sum the tags with the meth direction C--T G -- A
DMR_SNP_overlap_CT_GA_finalSummary <- DMR_SNP_overlap_CT_GA_SUM_tagged %>%
  group_by(Meth_Direction, ALT_comparison) %>%
  summarise(
    n = n()
  ) %>%
  ungroup()


##plot the summary
plot1 <- ggplot(DMR_SNP_overlap_CT_GA_finalSummary, aes(x = factor(Meth_Direction),
                                y = n,
                                fill = ALT_comparison)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    x = "Methylation Direction",
    y = "Number of alternative allele",
    fill = "ALT Allele Enrichment",
    title = "ALT Allele Dosage Comparison by Methylation Direction C->T, G->A"
  ) +
  geom_text(aes(label = n),
            position = position_dodge(width = 0.9),
            vjust = -0.3, size = 7, fontface="bold")+
  scale_fill_manual(values = c(
    "AA_less_ALT" = "black",
    "AA_more_ALT" = "grey30",
    "Equal_ALT" = "grey70"
  ))+
  scale_x_discrete(labels = c("-1" = "Hypo", "1" = "Hyper")) +
  theme_bw() +
  theme(
    title = element_text(face = "bold",size = 25),
    legend.title = element_text(face = "bold",size = 18),
    legend.text = element_text(face = "bold",size = 16),
    axis.title = element_text(face = "bold",size = 18),
    axis.text = element_text(face = "bold",size=16),
    axis.text.x = element_text(size=16),
    strip.text =element_text(face = "bold",size=16))


##### now consider all snps together
DMR_SNP_overlap_edit_clean_SUM<- DMR_SNP_overlap_edit_clean %>% 
  group_by(chr, start_DMR, end_DMR, Meth_Direction) %>%
  summarise(
    n_snps = n(),
    sum_AA_DOS = sum(AA_DOS),
    sum_CA_DOS = sum(CA_DOS),
    AA_avg = mean(AA_DOS),
    CA_avg = mean(CA_DOS)
  ) %>%
  ungroup()


##Add a column to classify alt allele difference classify by sum
DMR_SNP_overlap_edit_clean_SUM_tagged <- DMR_SNP_overlap_edit_clean_SUM %>%
  mutate(
    ALT_comparison = case_when(
      sum_AA_DOS > sum_CA_DOS ~ "AA_more_ALT",
      sum_AA_DOS < sum_CA_DOS ~ "AA_less_ALT",
      sum_AA_DOS == sum_CA_DOS ~ "Equal_ALT"
    )
  )

##Add a column to classify alt allele difference classify by mean
DMR_SNP_overlap_edit_clean_mean_tagged <- DMR_SNP_overlap_edit_clean_SUM%>%
  mutate(
    ALT_comparison = case_when(
      sum_AA_DOS > sum_CA_DOS ~ "AA_more_ALT",
      sum_AA_DOS < sum_CA_DOS ~ "AA_less_ALT",
      sum_AA_DOS == sum_CA_DOS ~ "Equal_ALT"
    )
  )





### finally sum the tags with the meth direction, ALL SNPS
DMR_SNP_overlap_edit_clean_finalSummary <- DMR_SNP_overlap_edit_clean_mean_tagged %>%
  group_by(Meth_Direction, ALT_comparison) %>%
  summarise(
    n = n()
  ) %>%
  ungroup()



plot2 <- ggplot(DMR_SNP_overlap_edit_clean_finalSummary, aes(x = factor(Meth_Direction),
                                                        y = n,
                                                        fill = ALT_comparison)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    x = "Methylation Direction",
    y = "Number of alternative allele",
    fill = "ALT Allele Enrichment",
    title = "ALT Allele Dosage Comparison by Methylation Direction (All SNPs)"
  ) +   
  geom_text(aes(label = n),
    position = position_dodge(width = 0.9),
    vjust = -0.3, size = 7, fontface="bold")+
  scale_fill_manual(values = c(
    "AA_less_ALT" = "black",
    "AA_more_ALT" = "grey30",
    "Equal_ALT" = "grey70"
  ))+
  scale_x_discrete(labels = c("-1" = "Hypo", "1" = "Hyper")) +
  theme_bw() +
  theme(
    title = element_text(face = "bold",size = 25),
    legend.title = element_text(face = "bold",size = 18),
    legend.text = element_text(face = "bold",size = 16),
    axis.title = element_text(face = "bold",size = 18),
    axis.text = element_text(face = "bold",size=16),
    axis.text.x = element_text(size=16),
    strip.text =element_text(face = "bold",size=16))


gridExtra::grid.arrange(plot1,plot2)

##consider for T-->C and A-->C SNPs
DMR_SNP_overlap_edit_clean_TC_AG<- DMR_SNP_overlap_edit_clean %>% filter(
  (DMR_SNP_overlap_edit_clean$REF == "T" & DMR_SNP_overlap_edit_clean$ALT =="C") |
    (DMR_SNP_overlap_edit_clean$REF == "A" & DMR_SNP_overlap_edit_clean$ALT =="G") 
)

##function to automate
SumMeanDosagesEachDMR<- function(DMR_SNP_TABLE){
  
  ##calculate SUM and Mean for all SNPs within each DMR
  DMR_SNP_TABLE_SUM_MEAN<- DMR_SNP_TABLE %>% 
    group_by(chr, start_DMR, end_DMR, Meth_Direction) %>%
    summarise(
      n_snps = n(),
      sum_AA_DOS = sum(AA_DOS),
      sum_CA_DOS = sum(CA_DOS),
      AA_avg = mean(AA_DOS),
      CA_avg = mean(CA_DOS)
    ) %>%
    ungroup()
  
  #tag them
  DMR_SNP_TABLE_SUM_tagged <- DMR_SNP_TABLE_SUM_MEAN%>%
    mutate(
      ALT_comparison = case_when(
        sum_AA_DOS > sum_CA_DOS ~ "AA > CA",
        sum_AA_DOS < sum_CA_DOS ~ "AA < CA",
        sum_AA_DOS == sum_CA_DOS ~ "AA = CA"
      )
    )
  
  DMR_SNP_TABLE_MEAN_tagged <- DMR_SNP_TABLE_SUM_MEAN%>%
    mutate(
      ALT_comparison = case_when(
        AA_avg > CA_avg ~ "AA > CA",
        AA_avg < CA_avg ~ "AA < CA",
        AA_avg == CA_avg ~ "AA = CA"
      )
    )
  
  
  return(list(
    SUM_TAGGED= DMR_SNP_TABLE_SUM_tagged,
    MEAN_TAGGED= DMR_SNP_TABLE_MEAN_tagged
  ))
  
}
DMR_SNP_overlap_TC_AG_SUM_MEAN<- SumMeanDosagesEachDMR(DMR_SNP_overlap_edit_clean_TC_AG)


### finally sum the tags with the meth direction, T-->C A-->G
DMR_SNP_overlap_TC_AG_finalSummary <- DMR_SNP_overlap_TC_AG_SUM_MEAN$MEAN_TAGGED %>%
  group_by(Meth_Direction, ALT_comparison) %>%
  summarise(
    n = n()
  ) %>%
  ungroup()



##consider only for SP6A

SP6A_RegionDMR<- rbind(
DMR_DEGs_regions_BT_AA_comparisonGeno$Upstream %>% dplyr::filter(aktualisierte.annotation..bitte.ergänzen. == "StPEBP9, SP6A"),
DMR_DEGs_regions_BT_AA_comparisonGeno$Downstream %>% dplyr::filter(aktualisierte.annotation..bitte.ergänzen. == "StPEBP9, SP6A"))


SP6A_RegionDMR_GR <- GRanges(
  seqnames = SP6A_RegionDMR$chr,
  ranges = IRanges(start = SP6A_RegionDMR$start_DMR, 
                   end = SP6A_RegionDMR$end_DMR),
  locusName = SP6A_RegionDMR$locusName,
  score = SP6A_RegionDMR$Methylation_direction,
  Description = SP6A_RegionDMR$v6.1.Description
)


export(SP6A_RegionDMR_GR, "/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/SNP_DMRs/SP6A_RegionDMR_GR.bed", 
       format = "bed")

##### consider all SNPs now and all DMRs
ALL_DMR_SNP_overlap <- read.delim("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/SNP_DMRs/ALL_DMR_SNP_overlap.tsv", header=FALSE)

ALL_DMR_SNP_overlap
ALL_DMR_SNP_overlap_edit<- ALL_DMR_SNP_overlap[,c(1:7,10:16)]
colnames(ALL_DMR_SNP_overlap_edit) <- c("chr","POS","V3","REF","ALT","QUAL","V7","AA","CA","Chr_dmr","start_DMR","end_DMR","V15","Meth_Direction")

ALL_DMR_SNP_overlap_edit %>% select(start_DMR, end_DMR) %>% unique() %>% nrow()

#DMR_SNP_overlap_edit_clean <- inner_join(DMR_SNP_overlap_edit_clean,
#                                         Enrichment_GB_candidates_DF_edit,
#                                         by=c("chr","start_DMR","end_DMR"))
####Consider now for All upstream downstream and Genebody DMR DEGs
DMR_DEGs_regions_BT_AA_comparisonGeno<- list(
  Upstream = read_excel("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/AAvsCA_Controls_Summary_1st2nd_EXP/DMR_DEGs_Controls_BT/DMR_DEGs_regions_BT_AA_comparisonGeno.xlsx", sheet = "Upstream"),
  Genebody = read_excel("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/AAvsCA_Controls_Summary_1st2nd_EXP/DMR_DEGs_Controls_BT/DMR_DEGs_regions_BT_AA_comparisonGeno.xlsx", sheet = "Genebody"),
  Downstream = read_excel("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/AAvsCA_Controls_Summary_1st2nd_EXP/DMR_DEGs_Controls_BT/DMR_DEGs_regions_BT_AA_comparisonGeno.xlsx", sheet = "Downstream")
)

#Create a granges list to save, consider all methylation context....
GR_list_DMR_DEGs_AA_Comparison <- list()
for (regions in names(DMR_DEGs_regions_BT_AA_comparisonGeno)) {
  
  GR_list_DMR_DEGs_AA_Comparison[[regions]] <- GRanges(
    seqnames = DMR_DEGs_regions_BT_AA_comparisonGeno[[regions]]$chr,
    ranges = IRanges(start = DMR_DEGs_regions_BT_AA_comparisonGeno[[regions]]$start_DMR, 
                     end = DMR_DEGs_regions_BT_AA_comparisonGeno[[regions]]$end_DMR),
    locusName = DMR_DEGs_regions_BT_AA_comparisonGeno[[regions]]$locusName,
    score = DMR_DEGs_regions_BT_AA_comparisonGeno[[regions]]$Methylation_direction,
    Description = DMR_DEGs_regions_BT_AA_comparisonGeno[[regions]]$v6.1.Description
  )
  
}

##expore as bedfile
for (regions in names(GR_list_DMR_DEGs_AA_Comparison)) {
  
  export(GR_list_DMR_DEGs_AA_Comparison[[regions]], 
         paste0("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/SNP_DMRs/BED_files_DMR_DEGs/GR_list_DMR_DEGs_AA_Comparison_",regions,".bed"), 
         format = "bed")
}
 
##overlapping done by bedtools intersect
###load the overlap table
DMR_DEG_SNP_overlap_list<- list(
  Upstream = read.delim("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/SNP_DMRs/BED_files_DMR_DEG_overlaps/GR_list_DMR_DEGs_AA_Comparison_Upstream_DMRDEG_overlaps.tsv", header=FALSE),
  Genebody = read.delim("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/SNP_DMRs/BED_files_DMR_DEG_overlaps/GR_list_DMR_DEGs_AA_Comparison_Genebody_DMRDEG_overlaps.tsv", header=FALSE),
  Downstream = read.delim("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/SNP_DMRs/BED_files_DMR_DEG_overlaps/GR_list_DMR_DEGs_AA_Comparison_Downstream_DMRDEG_overlaps.tsv", header=FALSE)
)

sum(is.na(DMR_DEG_SNP_overlap_list$Genebody))
View(DMR_DEG_SNP_overlap_list$Genebody)
##function to get clean tables
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

test <- DMR_DEG_SNP_overlap_list$Genebody


test_edit <- test[,c(1:7,10:16)]
colnames(test_edit) <- c("chr","POS","V3","REF","ALT","QUAL","V7","AA","CA","Chr_dmr","start_DMR","end_DMR","V15","Meth_Direction")

test_edit$AA_DOS <- sapply(test_edit$AA, extract_dosage)
test_edit$CA_DOS <- sapply(test_edit$CA, extract_dosage)

# Find indices where NAs occur

###NA occurs due to cant determin dosage for either AA or CA
na_aa_indices <- which(is.na(test_edit$AA_DOS))
na_ca_indices <- which(is.na(test_edit$CA_DOS))

test_edit[na_aa_indices, ]
test_edit[na_ca_indices, ]


##sum dosages.. and clean table
DMR_DEG_SNP_overlap_clean_list<- lapply(DMR_DEG_SNP_overlap_list, CleanTables)

DMR_DEG_SNP_overlap_clean_SUM_MEAN_list<- lapply(DMR_DEG_SNP_overlap_clean_list, SumMeanDosagesEachDMR)

SUMMARY_table_func<- function(SUM_MEAN_table){
  
  SUM_MEAN_tablefinalSummary <- SUM_MEAN_table %>% na.omit() %>%
    group_by(Meth_Direction, ALT_comparison) %>%
    summarise(
      n = n()
    ) %>%
    ungroup()
  
  return(SUM_MEAN_tablefinalSummary)
  
}
DMR_DEG_SNP_overlap_clean_SUM_MEAN_list$Upstream

##Summarise
DMR_DEG_summary_Upstream<- lapply(DMR_DEG_SNP_overlap_clean_SUM_MEAN_list$Upstream, SUMMARY_table_func)
DMR_DEG_summary_Genebody<- lapply(DMR_DEG_SNP_overlap_clean_SUM_MEAN_list$Genebody, SUMMARY_table_func)
DMR_DEG_summary_Downstream<- lapply(DMR_DEG_SNP_overlap_clean_SUM_MEAN_list$Downstream, SUMMARY_table_func)

##Check values...
DMR_DEG_summary_Upstream$SUM_TAGGED
DMR_DEG_summary_Upstream$MEAN_TAGGED

DMR_DEG_summary_Genebody$SUM_TAGGED
DMR_DEG_summary_Genebody$MEAN_TAGGED

DMR_DEG_summary_Downstream$SUM_TAGGED
DMR_DEG_summary_Downstream$MEAN_TAGGED

rbind(DMR_DEG_summary_Upstream$SUM_TAGGED %>% mutate(Region="Upstream"),
DMR_DEG_summary_Genebody$SUM_TAGGED %>% mutate(Region="Genebody"),
DMR_DEG_summary_Downstream$SUM_TAGGED %>% mutate(Region="Downstream"))

### Now consider DEGs in the mix
##first need to match the respective dmrs to gene locus
DMR_DEGs_regions_BT_AA_comparisonGeno
DMR_DEG_SNP_overlap_clean_SUM_MEAN_list$Upstream
DMR_DEG_SNP_overlap_clean_SUM_MEAN_list$Genebody
DMR_DEG_SNP_overlap_clean_SUM_MEAN_list$Downstream

AnnotateDMRsSNPdosages2GeneLocus<- function(SnpDMR_Table_overlap, Original_DMR_DEG_Table){
  
  ### bedtools have an 0 starting system..
  Original_DMR_DEG_Table$start_DMR <- Original_DMR_DEG_Table$start_DMR - 1
  
  ### select only the key columns from the original table for the merge
  Original_DMR_DEG_Table_select <- Original_DMR_DEG_Table %>% select(chr, start_DMR, end_DMR, locusName, 
                                                                     Context, Methylation_direction, 
                                                                     log2FoldChange, regulation) %>% dplyr::rename(
                                                                       Meth_Direction = Methylation_direction)
 
   SnpDMR_Table_overlap_annot <- inner_join(SnpDMR_Table_overlap, Original_DMR_DEG_Table_select,
                                           by=c("chr","start_DMR","end_DMR","Meth_Direction"))
  
  
   return(SnpDMR_Table_overlap_annot)
}

## create list just for the SUM table
DMR_DEG_SNP_overlap_clean_SUM_list <- list(
  Upstream = DMR_DEG_SNP_overlap_clean_SUM_MEAN_list$Upstream$SUM_TAGGED,
  Genebody = DMR_DEG_SNP_overlap_clean_SUM_MEAN_list$Genebody$SUM_TAGGED,
  Downstream = DMR_DEG_SNP_overlap_clean_SUM_MEAN_list$Downstream$SUM_TAGGED
)

which(is.na(DMR_DEG_SNP_overlap_clean_SUM_MEAN_list$Genebody$SUM_TAGGED))
DMR_DEG_SNP_overlap_clean_SUM_list$Genebody[which(is.na(DMR_DEG_SNP_overlap_clean_SUM_list$Genebody)), ]


#Annotation 
DMR_DEG_SNP_overlap_clean_SUM_list_annot <- list()
for (regions in names(DMR_DEG_SNP_overlap_clean_SUM_list)) {

  ##run the annotation function, like lapply
  DMR_DEG_SNP_overlap_clean_SUM_list_annot[[regions]] <-  AnnotateDMRsSNPdosages2GeneLocus(
    DMR_DEG_SNP_overlap_clean_SUM_list[[regions]],
    DMR_DEGs_regions_BT_AA_comparisonGeno[[regions]]
  )
  
}

###NAs are present in the genebody tables because one sample for a given SNP, 
#you cant determine the SNP GT, thus cant determine the dosages
test <- DMR_DEG_SNP_overlap_clean_SUM_list_annot$Genebody[which(is.na(DMR_DEG_SNP_overlap_clean_SUM_list_annot$Genebody)), ]


### some DMRs are in the vicinity of two genes...  note...
nrow(DMR_DEG_SNP_overlap_clean_SUM_list$Upstream)
nrow(DMR_DEG_SNP_overlap_clean_SUM_list$Genebody)
nrow(DMR_DEG_SNP_overlap_clean_SUM_list$Downstream)

nrow(DMR_DEG_SNP_overlap_clean_SUM_list_annot$Upstream)
nrow(DMR_DEG_SNP_overlap_clean_SUM_list_annot$Genebody)
nrow(DMR_DEG_SNP_overlap_clean_SUM_list_annot$Downstream)

###Sum all, without regulation....
DMR_DEG_SNP_overlap_clean_SUM_list_summary <- list()
for (regions in names(DMR_DEG_SNP_overlap_clean_SUM_list)) {
  
  DMR_DEG_SNP_overlap_clean_SUM_list_summary[[regions]] <-  DMR_DEG_SNP_overlap_clean_SUM_list[[regions]] %>% 
    na.omit() %>%
    group_by(Meth_Direction, ALT_comparison) %>%
    summarise(
      n = n()
    ) %>%
    ungroup()
}

###Sum all
DMR_DEG_SNP_overlap_clean_SUM_list_annot_summary2 <- list()
for (regions in names(DMR_DEG_SNP_overlap_clean_SUM_list_annot)) {
  
  DMR_DEG_SNP_overlap_clean_SUM_list_annot_summary2[[regions]] <-  DMR_DEG_SNP_overlap_clean_SUM_list_annot[[regions]] %>% 
    na.omit() %>%
    group_by(ALT_comparison,regulation) %>%
    summarise(
      n = n()
    ) %>%
    ungroup()
}


glm(n ~ Meth_Direction * ALT_comparison * regulation, 
    family = poisson, data = DMR_DEG_SNP_overlap_clean_SUM_list_annot_summary$Upstream)

##summary
#Downregulation
#hypomethylation --> more ALT allele
#Hypermethylation --> less AlT allele

#Upregulation
##Hypermethyation --> less alt allele

###plot everything, without regulation
DMR_DEG_SNP_plot<- rbind(DMR_DEG_SNP_overlap_clean_SUM_list_summary$Upstream %>% mutate(Region="Upstream"),
                          DMR_DEG_SNP_overlap_clean_SUM_list_summary$Genebody %>% mutate(Region="Genebody"),
                          DMR_DEG_SNP_overlap_clean_SUM_list_summary$Downstream %>% mutate(Region="Downstream"))

DMR_DEG_SNP_plot$ALT_comparison <- factor(DMR_DEG_SNP_plot$ALT_comparison, levels = c("AA < CA", "AA > CA", "AA = CA"))

plot1<- ggplot(DMR_DEG_SNP_plot, 
       aes(x = factor(Meth_Direction),
          y = n,
          fill = ALT_comparison)) +
  facet_grid(~Region)+
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    x = "Methylation Direction",
    y = "Number of DMRs (unique regions)",
    fill = "ALT Allele Enrichment",
    title = "ALT Allele Dosage Comparison by Methylation Direction (All SNPs)"
  ) +   
  geom_text(aes(label = n),
            position = position_dodge(width = 0.9),
            vjust = -0.3, size = 7, fontface="bold")+
  scale_fill_manual(values = c(
    "AA < CA" = "black",
    "AA > CA" = "grey30",
    "AA = CA" = "grey70"))+
  scale_x_discrete(labels = c("-1" = "Hypo", "1" = "Hyper")) +
  theme_bw() +
  theme(
    title = element_text(face = "bold",size = 25),
    legend.title = element_text(face = "bold",size = 18),
    legend.text = element_text(face = "bold",size = 16),
    axis.title = element_text(face = "bold",size = 18),
    axis.text = element_text(face = "bold",size=16),
    axis.text.x = element_text(size=16),
    strip.text =element_text(face = "bold",size=16))

ggsave("DMR_DEG_SNP_plot.svg",
  plot = plot1,
  dpi = 300,
  path = "/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/SNP_DMRs/SNP_DMR/DMR_DEG_SNP_summary/",
  device = "svg", 
  width = 18, 
  height = 10
)

###plot Alt comparison with regulation
ggplot(rbind(DMR_DEG_SNP_overlap_clean_SUM_list_annot_summary2$Upstream %>% mutate(Region="Upstream"),
             DMR_DEG_SNP_overlap_clean_SUM_list_annot_summary2$Genebody %>% mutate(Region="Genebody"),
             DMR_DEG_SNP_overlap_clean_SUM_list_annot_summary2$Downstream %>% mutate(Region="Downstream")), 
       aes(x = factor(regulation),
           y = n,
           fill = ALT_comparison)) +
  facet_grid(~Region)+
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    x = "Regulation",
    y = "Number of alternative allele",
    fill = "ALT Allele Enrichment",
    title = "ALT Allele Dosage Comparison by Methylation Direction (All SNPs)"
  ) +   
  geom_text(aes(label = n),
            position = position_dodge(width = 0.9),
            vjust = -0.3, size = 7, fontface="bold")+
  scale_fill_manual(values = c(
    "AA_less_ALT" = "black",
    "AA_more_ALT" = "grey30",
    "Equal_ALT" = "grey70"
  ))+
  theme_bw() +
  theme(
    title = element_text(face = "bold",size = 25),
    legend.title = element_text(face = "bold",size = 18),
    legend.text = element_text(face = "bold",size = 16),
    axis.title = element_text(face = "bold",size = 18),
    axis.text = element_text(face = "bold",size=16),
    axis.text.x = element_text(size=16),
    strip.text =element_text(face = "bold",size=16))

##TRUE for all regions, effect is strongest in upstream and downstream regions


### consider on Context specific methylation....
DMR_DEG_SNP_overlap_clean_SUM_list_annot
DMR_DEG_SNP_overlap_clean_SUM_list_annot_summary <- list()
for (regions in names(DMR_DEG_SNP_overlap_clean_SUM_list_annot)) {
  
  DMR_DEG_SNP_overlap_clean_SUM_list_annot_summary[[regions]] <-  DMR_DEG_SNP_overlap_clean_SUM_list_annot[[regions]] %>% 
    na.omit() %>%
    group_by(Meth_Direction, Context, ALT_comparison, regulation) %>%
    summarise(
      n = n()
    ) %>%
    ungroup()
}

ggplot(DMR_DEG_SNP_overlap_clean_SUM_list_annot_summary$Genebody, 
       aes(x = factor(Meth_Direction),
           y = n,
           fill = ALT_comparison)) +
  facet_grid(regulation~Context)+
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    x = "Methylation Direction",
    y = "Number of alternative allele",
    fill = "ALT Allele Enrichment",
    title = "ALT Allele Dosage Comparison by Methylation Direction (All SNPs)"
  ) +   
  geom_text(aes(label = n),
            position = position_dodge(width = 0.9),
            vjust = -0.3, size = 7, fontface="bold")+
  scale_fill_manual(values = c(
    "AA_less_ALT" = "black",
    "AA_more_ALT" = "grey30",
    "Equal_ALT" = "grey70"
  ))+
  scale_x_discrete(labels = c("-1" = "Hypo", "1" = "Hyper")) +
  theme_bw() +
  theme(
    title = element_text(face = "bold",size = 25),
    legend.title = element_text(face = "bold",size = 18),
    legend.text = element_text(face = "bold",size = 16),
    axis.title = element_text(face = "bold",size = 18),
    axis.text = element_text(face = "bold",size=16),
    axis.text.x = element_text(size=16),
    strip.text =element_text(face = "bold",size=16))


### consider on alt allele differences and regulation only....
DMR_DEG_SNP_overlap_clean_SUM_list_annot
DMR_DEG_SNP_overlap_clean_SUM_list_annot_summary <- list()
for (regions in names(DMR_DEG_SNP_overlap_clean_SUM_list_annot)) {
  
  DMR_DEG_SNP_overlap_clean_SUM_list_annot_summary[[regions]] <-  DMR_DEG_SNP_overlap_clean_SUM_list_annot[[regions]] %>% 
    na.omit() %>%
    group_by(ALT_comparison, regulation) %>%
    summarise(
      n = n()
    ) %>%
    ungroup()
}

ggplot(rbind(DMR_DEG_SNP_overlap_clean_SUM_list_annot_summary$Upstream %>% mutate(Region="Upstream"),
             DMR_DEG_SNP_overlap_clean_SUM_list_annot_summary$Genebody %>% mutate(Region="Genebody"),
             DMR_DEG_SNP_overlap_clean_SUM_list_annot_summary$Downstream %>% mutate(Region="Downstream")),
       aes(x = factor(ALT_comparison),
           y = n,
           fill = ALT_comparison)) +
  facet_grid(Region~regulation)+
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    x = "Methylation Direction",
    y = "Number of alternative allele",
    fill = "ALT Allele Enrichment",
    title = "ALT Allele Dosage Comparison by Methylation Direction (All SNPs)"
  ) +   
  geom_text(aes(label = n),
            position = position_dodge(width = 0.9),
            vjust = -0.3, size = 7, fontface="bold")+
  scale_fill_manual(values = c(
    "AA_less_ALT" = "black",
    "AA_more_ALT" = "grey30",
    "Equal_ALT" = "grey70"
  ))+
  scale_x_discrete(labels = c("-1" = "Hypo", "1" = "Hyper")) +
  scale_y_continuous(
    breaks = c(0, 100, 200, 300, 400, 500),
    limits = c(0, 600),
    expand = c(0, 0) # removes extra whitespace at top/bottom
  )+
  theme_bw() +
  theme(
    title = element_text(face = "bold",size = 25),
    legend.title = element_text(face = "bold",size = 18),
    legend.text = element_text(face = "bold",size = 16),
    axis.title = element_text(face = "bold",size = 18),
    axis.text = element_text(face = "bold",size=16),
    axis.text.x = element_text(size=16),
    strip.text =element_text(face = "bold",size=16))


View(DMR_DEG_SNP_overlap_clean_SUM_list_annot$Upstream)

sum(is.na(DMR_DEG_SNP_overlap_clean_SUM_list_annot$Upstream))
sum(is.na(DMR_DEG_SNP_overlap_clean_SUM_list_annot$Genebody))
sum(is.na(DMR_DEG_SNP_overlap_clean_SUM_list_annot$Downstream))


####Consider DMR GENES and not DEGs.... #####
###load DMR all GENEs

#CG DMR GENES
CG_DMR_Genes <- list(
Upstream = read.delim("/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/WGBS_HIGH_HEAT/ANALYSIS/CpG_DMRs/FISHER/BEDTOOLS_Controls/CG_UpstreamBEDTOOLS.bed", header=FALSE),
Genebody = read.delim("/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/WGBS_HIGH_HEAT/ANALYSIS/CpG_DMRs/FISHER/BEDTOOLS_Controls/CG_GenebodyBEDTOOLS.bed", header=FALSE),
Downstream = read.delim("/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/WGBS_HIGH_HEAT/ANALYSIS/CpG_DMRs/FISHER/BEDTOOLS_Controls/CG_DownstreamBEDTOOLS.bed", header=FALSE))

#CHG DMR GENES
CHG_DMR_Genes <- list(
  Upstream = read.delim("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/CHG_DMRs/FISHER/BEDTOOLS_Controls/CHG_UpstreamBEDTOOLS.bed", header=FALSE),
  Genebody = read.delim("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/CHG_DMRs/FISHER/BEDTOOLS_Controls/CHG_GenebodyBEDTOOLS.bed", header=FALSE),
  Downstream = read.delim("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/CHG_DMRs/FISHER/BEDTOOLS_Controls/CHG_DownstreamBEDTOOLS.bed", header=FALSE))

#CHH DMR GENES
CHH_DMR_Genes <- list(
  Upstream = read.delim("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/CHH_DMRs/FISHER/BEDTOOLS_Controls/CHH_UpstreamBEDTOOLS.bed", header=FALSE),
  Genebody = read.delim("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/CHH_DMRs/FISHER/BEDTOOLS_Controls/CHH_GenebodyBEDTOOLS.bed", header=FALSE),
  Downstream = read.delim("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/CHH_DMRs/FISHER/BEDTOOLS_Controls/CHH_DownstreamBEDTOOLS.bed", header=FALSE))

RenameEditColumns <- function(DMR_GENE_Table, CONTEXT){
  
  DMR_GENE_Table_select <- DMR_GENE_Table %>% select(
    V1, V2, V3, V5, V10, V11, V13, V15) %>% dplyr::rename(
     chr = V1, 
     start_DMR = V2, 
     end_DMR = V3, 
     Meth_Direction = V5, 
     start = V10, 
     end = V11, 
     strand = V13, 
     locusName = V15)
  
  DMR_GENE_Table_select_edit <- DMR_GENE_Table_select %>% separate(
    locusName, into = c("ID", "locusName"), sep=";"
  ) %>% select(-ID) %>% select(chr, start_DMR, end_DMR, Meth_Direction, Context, start, end, strand, locusName)
  
  DMR_GENE_Table_select_edit$locusName <- gsub("Name=","",DMR_GENE_Table_select_edit$locusName)
  
  ###reverse the meth direction so AA becomes the comparison and not CA
  DMR_GENE_Table_select_edit<-  DMR_GENE_Table_select_edit %>%  
    mutate(Meth_Direction = Meth_Direction * -1)
  
  
  return(DMR_GENE_Table_select_edit)
}
?separate

CG_DMR_Genes <- lapply(CG_DMR_Genes, RenameEditColumns, CONTEXT = "CG")
CHG_DMR_Genes <- lapply(CHG_DMR_Genes, RenameEditColumns, CONTEXT = "CHG")
CHH_DMR_Genes <- lapply(CHH_DMR_Genes, RenameEditColumns, CONTEXT = "CHH")

##rearrange to regions
Upstream_AA_compare_DMR_GENEs <- rbind(
  CG_DMR_Genes$Upstream,
  CHG_DMR_Genes$Upstream,
  CHH_DMR_Genes$Upstream)

Genebody_AA_compare_DMR_GENEs <- rbind(
  CG_DMR_Genes$Genebody,
  CHG_DMR_Genes$Genebody,
  CHH_DMR_Genes$Genebody)

Downstream_AA_compare_DMR_GENEs <- rbind(
  CG_DMR_Genes$Downstream,
  CHG_DMR_Genes$Downstream,
  CHH_DMR_Genes$Downstream)

nrow(Upstream_AA_compare_DMR_GENEs)
nrow(Genebody_AA_compare_DMR_GENEs)
nrow(Downstream_AA_compare_DMR_GENEs)

AA_compare_DMR_GENEs<- list(
  Upstream = Upstream_AA_compare_DMR_GENEs,
  Genebody = Genebody_AA_compare_DMR_GENEs,
  Downstream = Downstream_AA_compare_DMR_GENEs
)

library(rtracklayer)
#Create a granges list to save, consider all methylation context....
GR_list_DMR_GENEs_AA_Comparison <- list()
for (regions in names(AA_compare_DMR_GENEs)) {
  
  GR_list_DMR_GENEs_AA_Comparison[[regions]] <- GRanges(
    seqnames = AA_compare_DMR_GENEs[[regions]]$chr,
    ranges = IRanges(start = AA_compare_DMR_GENEs[[regions]]$start_DMR, 
                     end = AA_compare_DMR_GENEs[[regions]]$end_DMR),
    locusName = AA_compare_DMR_GENEs[[regions]]$locusName,
    score = AA_compare_DMR_GENEs[[regions]]$Meth_Direction
  )
  
}

##expore as bedfile
for (regions in names(GR_list_DMR_GENEs_AA_Comparison)) {
  
  export(GR_list_DMR_GENEs_AA_Comparison[[regions]], 
         paste0("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/SNP_DMRs/BED_files_DMR_GENEs/GR_list_DMR_GENEs_AA_Comparison_",regions,".bed"), 
         format = "bed")
}

DMR_GENEs_SNP_overlap_list<- list(
Upstream = read.delim("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/SNP_DMRs/BED_files_DMR_GENEs_Overlaps/DMR_GENEs_AA_Comparison_Upstream_SNP_overlaps.tsv", header=FALSE),
Genebody = read.delim("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/SNP_DMRs/BED_files_DMR_GENEs_Overlaps/DMR_GENEs_AA_Comparison_Genebody_SNP_overlaps.tsv", header=FALSE),
Downstream = read.delim("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/SNP_DMRs/BED_files_DMR_GENEs_Overlaps/DMR_GENEs_AA_Comparison_Downstream_SNP_overlaps.tsv", header=FALSE))

### total DMR GENEs Upstream --> 7779, 87.3% of DMRs consists of SNPs
DMR_GENEs_SNP_overlap_list$Upstream %>% select(V12,V13,V14) %>% unique() %>% nrow()

### total DMR GENEs Genebody --> 19557, 89.7% of DMRs consists of SNPs
DMR_GENEs_SNP_overlap_list$Genebody %>% select(V12,V13,V14) %>% unique() %>% nrow()

### total DMR GENEs Downstream --> 9465,  84.1% of DMRs consists of SNPs
DMR_GENEs_SNP_overlap_list$Downstream %>% select(V12,V13,V14) %>% unique() %>% nrow()


##sum dosages.. and clean table
DMR_GENEs_SNP_overlap_clean_list<- lapply(DMR_GENEs_SNP_overlap_list, CleanTables)

DMR_GENEs_SNP_overlap_clean_SUM_MEAN_list<- lapply(DMR_GENEs_SNP_overlap_clean_list, SumMeanDosagesEachDMR)
write.table(DMR_GENEs_SNP_overlap_clean_SUM_MEAN_list$Upstream$SUM_TAGGED,
            "/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/SNP_DMRs/SNP_DMR/DMR_GENEs_SNP_overlap_clean_SUM_MEAN_upstream.txt",
            col.names = TRUE,row.names = FALSE, quote =FALSE, sep= "\t")

write.table(DMR_GENEs_SNP_overlap_clean_SUM_MEAN_list$Genebody$SUM_TAGGED,
            "/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/SNP_DMRs/SNP_DMR/DMR_GENEs_SNP_overlap_clean_SUM_MEAN_GB.txt",
            col.names = TRUE,row.names = FALSE, quote =FALSE, sep= "\t")

write.table(DMR_GENEs_SNP_overlap_clean_SUM_MEAN_list$Downstream$SUM_TAGGED,
            "/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/SNP_DMRs/SNP_DMR/DMR_GENEs_SNP_overlap_clean_SUM_MEAN_downstream.txt",
            col.names = TRUE,row.names = FALSE, quote =FALSE, sep= "\t")


##Summarise
DMR_GENEs_summary_Upstream<- lapply(DMR_GENEs_SNP_overlap_clean_SUM_MEAN_list$Upstream, SUMMARY_table_func)
DMR_GENEs_summary_Genebody<- lapply(DMR_GENEs_SNP_overlap_clean_SUM_MEAN_list$Genebody, SUMMARY_table_func)
DMR_GENEs_summary_Downstream<- lapply(DMR_GENEs_SNP_overlap_clean_SUM_MEAN_list$Downstream, SUMMARY_table_func)

#
DMR_GENEs_summary_Upstream_sum <- DMR_GENEs_summary_Upstream$SUM_TAGGED %>% mutate(Region="Upstream")
DMR_GENEs_summary_Genebody_sum <- DMR_GENEs_summary_Genebody$SUM_TAGGED %>% mutate(Region="Genebody")
DMR_GENEs_summary_Downstream_sum <- DMR_GENEs_summary_Downstream$SUM_TAGGED %>% mutate(Region="Downstream")

##make DF for plot
DMR_GENES_SNP_plot <- rbind(DMR_GENEs_summary_Upstream_sum, 
                            DMR_GENEs_summary_Genebody_sum, 
                            DMR_GENEs_summary_Downstream_sum)

DMR_GENES_SNP_plot$ALT_comparison <- factor(DMR_GENES_SNP_plot$ALT_comparison, 
                                            levels = c("AA < CA", "AA > CA", "AA = CA"))


plot2<- ggplot(DMR_GENES_SNP_plot, 
               aes(x = factor(Meth_Direction),
                   y = n,
                   fill = ALT_comparison)) +
  facet_grid(~Region)+
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    x = "Methylation Direction",
    y = "Number of DMRs (unique regions)",
    fill = "ALT Allele Enrichment",
    title = "ALT Allele Dosage Comparison by Methylation Direction (All SNPs)"
  ) +   
  geom_text(aes(label = n),
            position = position_dodge(width = 0.9),
            vjust = -0.3, size = 7, fontface="bold")+
  scale_fill_manual(values = c(
    "AA < CA" = "black",
    "AA > CA" = "grey30",
    "AA = CA" = "grey70"))+
  scale_x_discrete(labels = c("-1" = "Hypo", "1" = "Hyper")) +
  theme_bw() +
  theme(
    title = element_text(face = "bold",size = 25),
    legend.title = element_text(face = "bold",size = 18),
    legend.text = element_text(face = "bold",size = 16),
    axis.title = element_text(face = "bold",size = 18),
    axis.text = element_text(face = "bold",size=16),
    axis.text.x = element_text(size=16),
    strip.text =element_text(face = "bold",size=16))

ggsave("DMR_GENEs_SNP_plot.svg",
       plot = plot2,
       dpi = 300,
       path = "/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/SNP_DMRs/SNP_DMR/DMR_DEG_SNP_summary/",
       device = "svg", 
       width = 18, 
       height = 10
)


#model <- glm(n ~ Meth_Direction * ALT_comparison, family = "poisson", 
#             data = DMR_GENEs_summary_Upstream_sum)
#summary(model)


##CHisq
pivot_wider(DMR_GENEs_summary_Upstream_sum)
df_wide <- DMR_GENEs_summary_Upstream_sum %>%
  pivot_wider(
    names_from = ALT_comparison,
    values_from = n
  ) %>% select(-Region) %>% column_to_rownames("Meth_Direction")

chisq.test(df_wide)


### CONSIDER C-->T and G-->A
###COnsider also T-->C and A-->G, because if Reference v6.1 is T, but the entire potato genetic library population is C, and say CA has a T and AA has a C allele,
## this will be reported in the VCF file as REF==T and ALT==C


####LOOOK INTO DMRS with CT GA TC and AG SNPS ONLY not both

# Allowed SNP types
allowed_snps <- function(df) {
  (df$REF == "C" & df$ALT == "T") |
    (df$REF == "G" & df$ALT == "A") |
    (df$REF == "T" & df$ALT == "C") |
    (df$REF == "A" & df$ALT == "G")
}

FilterRespectiveSNPs<- function(df){
  
  
  # Step 1: Find all unique DMRs
  all_dmrs <- df %>% distinct(start_DMR, end_DMR)
  
  # Step 2: For each DMR, check if all SNPs are allowed types
  dmr_snps_check <- df %>%
    group_by(start_DMR, end_DMR) %>%
    summarise(all_allowed = all(allowed_snps(cur_data_all())), .groups = 'drop')
  
  # Step 3: Keep only DMRs where all SNPs are allowed
  dmrs_only_allowed <- dmr_snps_check %>%
    filter(all_allowed) %>%
    select(start_DMR, end_DMR)
  
  # Step 4: Filter original data to keep only these DMRs
  filtered_df <- df %>%
    semi_join(dmrs_only_allowed, by = c("start_DMR", "end_DMR")) %>%
    # And keep only SNPs of allowed types
    filter(allowed_snps(.))
  
  return(filtered_df)
}

DMR_GENEs_SNP_CG_GA_TC_AG_overlap_clean_list<- lapply(DMR_GENEs_SNP_overlap_clean_list, 
                                                      FilterRespectiveSNPs)

##verify for some examples
#View(DMR_GENEs_SNP_CG_GA_TC_AG_overlap_clean_list$Genebody)
#View(DMR_GENEs_SNP_overlap_clean_list$Genebody)


## Sum the dosages again
DMR_GENEs_SNP_CG_GA_TC_AG_overlap_clean_SUM_MEAN_list<- lapply(DMR_GENEs_SNP_CG_GA_TC_AG_overlap_clean_list, 
                                                   SumMeanDosagesEachDMR)


##Summarise
DMR_GENEs_summary_Upstream<- lapply(DMR_GENEs_SNP_CG_GA_TC_AG_overlap_clean_SUM_MEAN_list$Upstream, SUMMARY_table_func)
DMR_GENEs_summary_Genebody<- lapply(DMR_GENEs_SNP_CG_GA_TC_AG_overlap_clean_SUM_MEAN_list$Genebody, SUMMARY_table_func)
DMR_GENEs_summary_Downstream<- lapply(DMR_GENEs_SNP_CG_GA_TC_AG_overlap_clean_SUM_MEAN_list$Downstream, SUMMARY_table_func)

#
DMR_GENEs_CT_GA_TC_AG_plot<- rbind(DMR_GENEs_summary_Upstream$SUM_TAGGED %>% mutate(Region="Upstream"),
      DMR_GENEs_summary_Genebody$SUM_TAGGED %>% mutate(Region="Genebody"),
      DMR_GENEs_summary_Downstream$SUM_TAGGED %>% mutate(Region="Downstream") )

DMR_GENEs_CT_GA_TC_AG_plot$ALT_comparison <- factor(DMR_GENEs_CT_GA_TC_AG_plot$ALT_comparison,
                                                    levels = c("AA < CA",
                                                               "AA > CA",
                                                               "AA = CA"
                                                               ))

plot3<- ggplot(DMR_GENEs_CT_GA_TC_AG_plot, 
               aes(x = factor(Meth_Direction),
                   y = n,
                   fill = ALT_comparison)) +
  facet_grid(~Region)+
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    x = "Methylation Direction",
    y = "Number of DMRs (unique regions)",
    fill = "ALT Allele Enrichment",
    title = "ALT Allele Dosage Comparison by Methylation Direction (C--T, G-->A, T-->C, A-->G)"
  ) +   
  geom_text(aes(label = n),
            position = position_dodge(width = 0.9),
            vjust = -0.3, size = 7, fontface="bold")+
  scale_fill_manual(values = c(
    "AA < CA" = "black",
    "AA > CA" = "grey30",
    "AA = CA" = "grey70"))+
  scale_x_discrete(labels = c("-1" = "Hypo", "1" = "Hyper")) +
  theme_bw() +
  theme(
    title = element_text(face = "bold",size = 25),
    legend.title = element_text(face = "bold",size = 18),
    legend.text = element_text(face = "bold",size = 16),
    axis.title = element_text(face = "bold",size = 18),
    axis.text = element_text(face = "bold",size=16),
    axis.text.x = element_text(size=16),
    strip.text =element_text(face = "bold",size=16))

ggsave("DMR_GENEs_CT_GA_TC_AG_SNP_plot.svg",
       plot = plot3,
       dpi = 300,
       path = "/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/SNP_DMRs/SNP_DMR/DMR_DEG_SNP_summary/",
       device = "svg", 
       width = 18, 
       height = 10
)

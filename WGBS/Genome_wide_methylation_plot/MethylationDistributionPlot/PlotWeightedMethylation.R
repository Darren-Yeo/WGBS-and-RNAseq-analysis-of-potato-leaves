library(tidyverse)

WeightedMethylation<- list(
  CG  = read.csv("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/MethylationPopStruct_1st2nd_Exp/MethylationDistributionPlot/tile200bp_CpG_weightMeth_na_omit.txt", sep=""),
  CHG  = read.csv("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/MethylationPopStruct_1st2nd_Exp/MethylationDistributionPlot/tile200bp_CHG_weightMeth_na_omit.txt", sep=""),
  CHH  = read.csv("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/MethylationPopStruct_1st2nd_Exp/MethylationDistributionPlot/tile200bp_CHH_weightMeth_na_omit.txt", sep="")
)

##initialize the dataframe with the column names from one of the dataframe
WeightedMethylation_DF <- data.frame(matrix(ncol = length(colnames(WeightedMethylation$CG))-1, nrow = 0))
colnames(WeightedMethylation_DF) <- colnames(WeightedMethylation$CG[2:ncol(WeightedMethylation$CG)])
##loop all dataframes to add contexts
for (context in names(WeightedMethylation)) {
  
  WeightedMethylation_in_loop_context<- WeightedMethylation[[context]] %>% select(-tile) %>% mutate(Context = context)
  
  ##bind all context df together
  WeightedMethylation_DF <-  rbind(WeightedMethylation_DF,WeightedMethylation_in_loop_context)
}

#tail(WeightedMethylation_DF)
#head(WeightedMethylation_DF)

#summary(WeightedMethylation$CHG$AA_28dap_C)
#summary(WeightedMethylation$CHG$AA_42dap_H)
#summary(WeightedMethylation$CHG$CA_28dap_C)
#summary(WeightedMethylation$CHG$CA_42dap_H)

#compare results with computemethylationprofile for DMRcaller 
##for CHG (most variable), median was about the same...
#summary(CHG_MethylationProfile_AA_noNA$Proportion)
#summary(CHG_MethylationProfile_AA_H_noNA$Proportion)
#summary(CHG_MethylationProfile_CA_noNA$Proportion)
#summary(CHG_MethylationProfile_CA_H_noNA$Proportion)


WeightedMethylation_DF_long <- WeightedMethylation_DF %>% 
  pivot_longer(cols = -Context,
               names_to = "Samples",
               values_to = "Weighted_methylation")

WeightedMethylation_DF_long<- WeightedMethylation_DF_long %>% mutate(
  sample_label=case_when(str_detect(Samples,"AA_28dap_C")~"AA",
                         str_detect(Samples,"AA_42dap_H")~"AA37C",
                         str_detect(Samples,"CA_28dap_C")~"CA",
                         str_detect(Samples,"CA_42dap_H")~"CA37C")
)

Plot<- ggplot(WeightedMethylation_DF_long,
       aes(x = sample_label, y = Weighted_methylation, fill = sample_label)) +
geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  scale_fill_manual(values = c(
    "AA" = "lightblue3",
    "AA37C" = "blue",
    "CA" = "pink3",
    "CA37C" = "red4"
  ))+
  facet_grid(. ~ Context) +
  theme_bw() +
  xlab("") +
  ylab("Methylation level %") +
  ggtitle("Weighted methylation % (200bp windows)") +
  labs(fill = "Samples")+
  theme(
    title = element_text(face = "bold",size = 25),
    legend.title = element_text(face = "bold",size = 18),
    legend.text = element_text(face = "bold",size = 16),
    axis.title = element_text(face = "bold",size = 18),
    axis.text = element_text(face = "bold",size=16),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,size=16),
    strip.text =element_text(face = "bold",size=16))


ggsave("WeightedMethylation200bp_window_boxplot.svg",
       path = "/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/MethylationPopStruct_1st2nd_Exp/MethylationDistributionPlot/",
       plot = Plot, 
       width = 14, 
       height = 12, 
       dpi = 300,
       device = "svg")


###wilcoxon ranked sum test, mann whitney test... heat vs control... before and after

#2870331 CG 200bp sites
nrow(WeightedMethylation$CG)

#2988118 CHG 200bp sites
nrow(WeightedMethylation$CHG)

#3308465 CHH 200bp sites
nrow(WeightedMethylation$CHH)

#### normality? 
##CG is rather bimodal, skewed towards 0 and 100

gridExtra::grid.arrange(
  plot1,
  plot2,
  plot3,
  plot4
)



plot1 <- hist(WeightedMethylation$CG$AA_28dap_C, breaks = 20, probability = TRUE, main = "CG_density")
plot2 <- hist(WeightedMethylation$CG$AA_42dap_H, breaks = 20, probability = TRUE, main = "CG_density")
plot3 <- hist(WeightedMethylation$CG$CA_28dap_C, breaks = 20, probability = TRUE, main = "CG_density")
plot4 <- hist(WeightedMethylation$CG$CA_42dap_H, breaks = 20, probability = TRUE, main = "CG_density")

##CHG also rather bimodal, althout more skewed towards 0
hist(WeightedMethylation$CHG$AA_28dap_C, breaks = 20, probability = TRUE, main = "CHG_density")
hist(WeightedMethylation$CHG$AA_42dap_H, breaks = 20, probability = TRUE, main = "CHG_density")
hist(WeightedMethylation$CHG$CA_28dap_C, breaks = 20, probability = TRUE, main = "CHG_density")
hist(WeightedMethylation$CHG$CA_42dap_H, breaks = 20, probability = TRUE, main = "CHG_density")

#CHH is skewed towards 0
hist(WeightedMethylation$CHH$AA_28dap_C, breaks = 20, probability = TRUE, main = "CHH_density")
hist(WeightedMethylation$CHH$AA_42dap_H, breaks = 20, probability = TRUE, main = "CHH_density")
hist(WeightedMethylation$CHH$CA_28dap_C, breaks = 20, probability = TRUE, main = "CHH_density")
hist(WeightedMethylation$CHH$CA_42dap_H, breaks = 20, probability = TRUE, main = "CHH_density")

### wilcoxon rank sum
ConvertDFtolong <- function(Wide_DF){
  
 Long_DF<-  Wide_DF %>% 
    select(AA_28dap_C, AA_42dap_H, CA_28dap_C, CA_42dap_H,Context) %>% 
    pivot_longer(
      cols = c("AA_28dap_C", "AA_42dap_H", "CA_28dap_C", "CA_42dap_H"),
      names_to = c("genotype", "timepoint", "condition"),
      names_pattern = "(AA|CA)_(\\d+dap)_(C|H)"
    )
  
  return(Long_DF)
}
library(rstatix)
WeightedMethylation_DF_long<- ConvertDFtolong(WeightedMethylation_DF)

histPlotlist <- list()
for (samples in unique(WeightedMethylation_DF_long$Samples)) {
  
  plot <- ggplot(WeightedMethylation_DF_long %>% filter(Samples==samples),
         aes(x = Weighted_methylation, color=Context)) + 
    geom_freqpoly(bins=50, position="identity")+
    scale_color_manual(values = c(
      "CG" = "#E69F00",   # orange
      "CHG" = "#56B4E9",  # blue
      "CHH" = "#009E73"   # green
    )) + 
    theme_bw() +
    theme(
      axis.title = element_text(size = 15, face = "bold"),
      axis.text = element_text(size = 15, face = "bold"),
      legend.title = element_text(size = 15, face = "bold"),
      legend.text = element_text(size = 14, face = "bold"),
      strip.text = element_text(size = 14, face = "bold")
    )
  
  histPlotlist[[samples]] <- plot
  
}

gridExtra::grid.arrange(histPlotlist$AA_28dap_C,
                        histPlotlist$AA_42dap_H,
                        histPlotlist$CA_28dap_C,
                        histPlotlist$CA_42dap_H)

results_wilcox <- WeightedMethylation_DF_long %>%
  group_by(Context,genotype) %>%
  wilcox_test(value ~ condition, exact = FALSE) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()


library(effsize)
results_cliffDelta <- WeightedMethylation_DF_long %>%
  group_by(Context, genotype) %>%
  summarise(cliff = list(cliff.delta(value, condition)), .groups = "drop") %>%
  unnest_wider(cliff)

### separate by context
WeightedMethylation_DF_CG<- WeightedMethylation_DF %>% dplyr::filter(Context=="CG")
WeightedMethylation_DF_CHG<- WeightedMethylation_DF %>% dplyr::filter(Context=="CHG")
WeightedMethylation_DF_CHH<- WeightedMethylation_DF %>% dplyr::filter(Context=="CHH")

####compute summary
summary(WeightedMethylation_DF_CG)
summary(WeightedMethylation_DF_CHG)
summary(WeightedMethylation_DF_CHH)

abs(median(WeightedMethylation_DF_CG$AA_42dap_H)- median(WeightedMethylation_DF_CG$AA_28dap_C))
abs(median(WeightedMethylation_DF_CG$CA_42dap_H)- median(WeightedMethylation_DF_CG$CA_28dap_C))

abs(median(WeightedMethylation_DF_CHG$AA_42dap_H)- median(WeightedMethylation_DF_CHG$AA_28dap_C))
abs(median(WeightedMethylation_DF_CHG$CA_42dap_H)- median(WeightedMethylation_DF_CHG$CA_28dap_C))

abs(median(WeightedMethylation_DF_CHH$AA_42dap_H)- median(WeightedMethylation_DF_CHH$AA_28dap_C))
abs(median(WeightedMethylation_DF_CHH$CA_42dap_H)- median(WeightedMethylation_DF_CHH$CA_28dap_C))


cliff.delta(WeightedMethylation_DF_CG_sample$AA_42dap_H,WeightedMethylation_DF_CG_sample$AA_28dap_C)
cliff.delta(WeightedMethylation_DF_CG_sample$CA_42dap_H,WeightedMethylation_DF_CG_sample$CA_28dap_C)

cliff.delta(WeightedMethylation_DF_CHG_sample$AA_42dap_H,WeightedMethylation_DF_CHG_sample$AA_28dap_C)
cliff.delta(WeightedMethylation_DF_CHG_sample$CA_42dap_H,WeightedMethylation_DF_CHG_sample$CA_28dap_C)

cliff.delta(WeightedMethylation_DF_CHH_sample$AA_42dap_H,WeightedMethylation_DF_CHH_sample$AA_28dap_C)
cliff.delta(WeightedMethylation_DF_CHH_sample$CA_42dap_H,WeightedMethylation_DF_CHH_sample$CA_28dap_C)
##save as df
write.table(results_wilcox,
            "/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/MethylationPopStruct_1st2nd_Exp/MethylationDistributionPlot/results_wilcox.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

### sample 10000 values to test
WeightedMethylation_DF_sample<- ConvertDFtolong(
  sample_n(WeightedMethylation_DF, size=100000, replace = FALSE)
)

results_sample_wilcox_approx_30000 <- WeightedMethylation_DF_sample %>%
  group_by(Context,genotype) %>%
  wilcox_test(value ~ condition, exact = FALSE) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()

write.table(results_sample_wilcox_approx_30000,
            "/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/MethylationPopStruct_1st2nd_Exp/MethylationDistributionPlot/results_sample_wilcox_approx_30000.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

#### Methylation within GBregions 100bp window 200bp within gene
Metagene_AA <- read.delim("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/MethylationPopStruct_1st2nd_Exp/BISMARKPLOT/merged/GB/Metagene_AA.tsv")
Metagene_AA37C <- read.delim("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/MethylationPopStruct_1st2nd_Exp/BISMARKPLOT/merged/GB/Metagene_AA37C.tsv")
Metagene_CA <- read.delim("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/MethylationPopStruct_1st2nd_Exp/BISMARKPLOT/merged/GB/Metagene_CA.tsv")
Metagene_CA37C <- read.delim("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/MethylationPopStruct_1st2nd_Exp/BISMARKPLOT/merged/GB/Metagene_CA37C.tsv")

metagene_list <- list(
  AA = Metagene_AA,
  AA37C = Metagene_AA37C,
  CA = Metagene_CA,
  CA37C = Metagene_CA37C
)
library(tidyverse)

# Clean and compute methylation level

CleanComputeMeanMeth<- function(METAGENE_LIST){
  
  METAGENE_summary_list <- list()
  for (samples in names(METAGENE_LIST)) {
    
    METAGENE_summary_list[[samples]] <- METAGENE_LIST[[samples]] %>%
      filter(count > 0) %>%
      group_by(fragment, context) %>%
      summarise(
        total_sum = sum(sum, na.rm = TRUE),
        total_count = sum(count, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(!!samples := total_sum / total_count) %>%
      select(fragment, context, !!samples)
    
  }
  return(METAGENE_summary_list)
}

metagene_summary_list<- CleanComputeMeanMeth(metagene_list)

merged_metagene <- reduce(metagene_summary_list, function(x, y) inner_join(x, y, by = c("fragment", "context")))

merged_metagene_long<- merged_metagene %>% pivot_longer(
  
  cols = c(-fragment,-context),
  names_to = "Samples",
  values_to = "Methylation"
  
) %>% mutate(
  Methylation= Methylation *100
)


plot_1 <-  ggplot(merged_metagene_long, 
                aes(x = fragment, y = Methylation, color = Samples)) +
  geom_line(alpha = 1) +
  scale_x_continuous(
    breaks = c(0, 100, 200, 300, 400),
    labels = c("", "TSS", "Gene body", "TES", "")
  ) +
  scale_color_manual(values = c(
    "AA" = "lightblue3",
    "AA37C" = "blue",
    "CA" = "pink3",
    "CA37C" = "red4"
  )) +
  facet_wrap(~ context, ncol = 3, scales = "free_y") +  # one column, separate y-axis
  theme_bw() +
  labs(
    title = "",
    x = "",
    y = "Methylation (%)",
    color = "Samples"
  ) +
  theme(
    axis.title = element_text(size = 15, face = "bold"),
    axis.text = element_text(size = 15, face = "bold"),
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold")
  )

ggsave(filename = paste0("All_context_GB_LandScape.svg"),
       path = "/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/MethylationPopStruct_1st2nd_Exp/BISMARKPLOT/merged/GB/",
       device = "svg",
       plot = plot_1, 
       width = 16, 
       height = 7, 
       dpi = 300)

for (CONTEXT in c("CG", "CHG", "CHH")) {
  
  ###plot
 plot <-  ggplot(merged_metagene_long %>% filter(context==CONTEXT), 
                 aes(x = fragment, y = Methylation, color = Samples)) +
    geom_line(alpha = 1) +
    scale_x_continuous(
      breaks = c(0, 100, 200,300, 400),        # change to your bin positions
      labels = c("Upstream", "TSS", "Gene body","TES", "Downstream")
    ) +
    scale_color_manual(values = c(
      "AA" = "lightblue3",  # greenish
      "AA37C" = "blue",  # orange
      "CA" = "pink3",  # purple
      "CA37C" = "red4"   # pink
    )) +
    theme_bw() +
    labs(
      title = "",
      x = "",
      y = "Methylation (%)",
      color = "Samples"
    ) +theme(
      axis.title = element_text(size = 15,face = "bold"),
      axis.text = element_text(size = 15,face = "bold"),
      legend.title = element_text(size = 15,face = "bold"),
      legend.text = element_text(size = 14,face = "bold")
    )
  
 ggsave(filename = paste0(CONTEXT, "_GB.svg"),
        path = "/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/MethylationPopStruct_1st2nd_Exp/BISMARKPLOT/merged/GB/",
        device = "svg",
        plot = plot, 
        width = 18, 
        height = 14, 
        dpi = 300)
 
}


plot<- ggplot(merged_metagene_long, aes(x = Samples, y = Methylation, fill = Samples)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  facet_wrap(~ context, nrow = 1) +
  scale_fill_manual(values = c(
    "AA" = "lightblue3",
    "AA37C" = "blue",
    "CA" = "pink3",
    "CA37C" = "red4"
  )) +
  labs(
    x = "",
    y = "Methylation (%)"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 15,face = "bold"),
    axis.text = element_text(size = 15,face = "bold"),
    legend.title = element_text(size = 15,face = "bold"),
    legend.text = element_text(size = 14,face = "bold"),
    strip.text = element_text(size = 14, face = "bold")
  ) 

ggsave(filename = paste0("Boxplot_GB.svg"),
       path = "/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/MethylationPopStruct_1st2nd_Exp/BISMARKPLOT/merged/GB/",
       device = "svg",
       plot = plot, 
       width = 15, 
       height = 13, 
       dpi = 300)


merged_metagene_CHH<- merged_metagene %>% filter(context=="CHH")
merged_metagene_CG<- merged_metagene %>% filter(context=="CG")
merged_metagene_CHG<- merged_metagene %>% filter(context=="CHG")

wilcox.test(merged_metagene_CHG$CA,
            merged_metagene_CHG$CA37C, paired = TRUE, exact = TRUE
            )
library(effsize)

cohen.d(merged_metagene_CHG$CA37C,
        merged_metagene_CHG$CA,
        paired = TRUE)

qqnorm(merged_metagene_CHH$AA)
qqline(merged_metagene_CHH$AA, col = "red")
### DO same thing for transposable elements
### Aggregation and summary methylation done in bash.. because too many rows..

metagene_TE_list <- list(
  AA = read.delim("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/MethylationPopStruct_1st2nd_Exp/BISMARKPLOT/merged/TE/Summary_Metagene_TE_AA.tsv"),
  AA37C = read.delim("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/MethylationPopStruct_1st2nd_Exp/BISMARKPLOT/merged/TE/Summary_Metagene_TE_AA37C.tsv"),
  CA = read.delim("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/MethylationPopStruct_1st2nd_Exp/BISMARKPLOT/merged/TE/Summary_Metagene_TE_CA.tsv"),
  CA37C = read.delim("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/MethylationPopStruct_1st2nd_Exp/BISMARKPLOT/merged/TE/Summary_Metagene_TE_CA37C.tsv")
)

for (df in names(metagene_TE_list)) {
  
  colnames(metagene_TE_list[[df]])[3] <- df
  
}

merged_metagene_TE <- reduce(metagene_TE_list, function(x, y) inner_join(x, y, by = c("fragment", "context")))

###convert ot long
merged_metagene_TE_long<- merged_metagene_TE %>% pivot_longer(
  
  cols = c(-fragment,-context),
  names_to = "Samples",
  values_to = "Methylation"
  
)

merged_metagene_TE_long <- merged_metagene_TE_long %>% 
  mutate(Methylation = Methylation * 100)

for (CONTEXT in c("CG", "CHG", "CHH")) {
  
  ###plot
  plot <-  ggplot(merged_metagene_TE_long %>% filter(context==CONTEXT), 
                  aes(x = fragment, y = Methylation, color = Samples)) +
    geom_line(alpha = 1) +
    scale_x_continuous(
      breaks = c(0, 100, 200,300, 400),        # change to your bin positions
      labels = c("Upstream", "TSS", "Repeat body","TES", "Downstream")
    ) +
    scale_color_manual(values = c(
      "AA" = "lightblue3",  # greenish
      "AA37C" = "blue",  # orange
      "CA" = "pink3",  # purple
      "CA37C" = "red4"   # pink
    )) +
    theme_bw() +
    labs(
      title = "",
      x = "",
      y = "Methylation (%)",
      color = "Samples"
    ) +theme(
      axis.title = element_text(size = 15,face = "bold"),
      axis.text = element_text(size = 15,face = "bold"),
      legend.title = element_text(size = 15,face = "bold"),
      legend.text = element_text(size = 14,face = "bold")
    )
  
  ggsave(filename = paste0(CONTEXT, "_TE.svg"),
         path = "/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/MethylationPopStruct_1st2nd_Exp/BISMARKPLOT/merged/TE",
         device = "svg",
         plot = plot, 
         width = 18, 
         height = 14, 
         dpi = 300)
  
}

plot_2 <-  ggplot(merged_metagene_TE_long, 
                aes(x = fragment, y = Methylation, color = Samples)) +
  geom_line(alpha = 1) +
  scale_x_continuous(
    breaks = c(0, 100, 200, 300, 400),
    labels = c("", "TSS", "Repeat body", "TES", "")
  ) +
  scale_color_manual(values = c(
    "AA" = "lightblue3",
    "AA37C" = "blue",
    "CA" = "pink3",
    "CA37C" = "red4"
  )) +
  facet_wrap(~ context, ncol = 3, scales = "free_y") +  # one column, separate y-axis
  theme_bw() +
  labs(
    title = "",
    x = "",
    y = "Methylation (%)",
    color = "Samples"
  ) +
  theme(
    axis.title = element_text(size = 15, face = "bold"),
    axis.text = element_text(size = 15, face = "bold"),
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold")
  )

ggsave(filename = paste0("All_context_TE_LandScape.svg"),
       path = "/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/MethylationPopStruct_1st2nd_Exp/BISMARKPLOT/merged/TE/",
       device = "svg",
       plot = plot_2, 
       width = 16, 
       height = 7, 
       dpi = 300)



###boxplot TEs
plot <- ggplot(merged_metagene_TE_long, aes(x = Samples, y = Methylation, fill = Samples)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  facet_wrap(~ context, nrow = 1) +
  scale_fill_manual(values = c(
    "AA" = "lightblue3",
    "AA37C" = "blue",
    "CA" = "pink3",
    "CA37C" = "red4"
  )) +
  labs(
    x = "",
    y = "Methylation (%)"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 15,face = "bold"),
    axis.text = element_text(size = 15,face = "bold"),
    legend.title = element_text(size = 15,face = "bold"),
    legend.text = element_text(size = 14,face = "bold"),
    strip.text = element_text(size = 14, face = "bold")
  ) 

ggsave(filename = paste0("Boxplot_TE.svg"),
       path = "/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/MethylationPopStruct_1st2nd_Exp/BISMARKPLOT/merged/TE/",
       device = "svg",
       plot = plot, 
       width = 15, 
       height = 13, 
       dpi = 300)


###save the metagene plot of GB and TE together
##GB plot
plot_1 <-  ggplot(merged_metagene_long, 
                  aes(x = fragment, y = Methylation, color = Samples)) +
  geom_line(alpha = 1) +
  scale_x_continuous(
    breaks = c(0, 100, 200, 300, 400),
    labels = c("Upstream", "TSS", "Gene body", "TES", "Downstream")
  ) +
  scale_color_manual(values = c(
    "AA" = "lightblue3",
    "AA37C" = "blue",
    "CA" = "pink3",
    "CA37C" = "red4"
  )) +
  facet_wrap(~ context, ncol = 1, scales = "free_y") +  # one column, separate y-axis
  theme_bw() +
  labs(
    title = "",
    x = "",
    y = "Methylation (%)",
    color = "Samples"
  ) +
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 17, face = "bold"),
    legend.position="none",
    legend.title = element_blank(),
    legend.text = element_blank(),
    strip.text = element_text(size = 18, face = "bold")
  )

#TE plot
plot_2 <-  ggplot(merged_metagene_TE_long, 
                  aes(x = fragment, y = Methylation, color = Samples)) +
  geom_line(alpha = 1) +
  scale_x_continuous(
    breaks = c(0, 100, 200, 300, 400),
    labels = c("Upstream", "TSS", "Repeat body", "TES", "Downstream")
  ) +
  scale_color_manual(values = c(
    "AA" = "lightblue3",
    "AA37C" = "blue",
    "CA" = "pink3",
    "CA37C" = "red4"
  )) +
  facet_wrap(~ context, ncol = 1, scales = "free_y") +  # one column, separate y-axis
  theme_bw() +
  labs(
    title = "",
    x = "",
    y = "Methylation (%)",
    color = "Samples"
  ) +
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    axis.title.y  = element_blank(),
    axis.text = element_text(size = 17, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 17, face = "bold"),
    strip.text = element_text(size = 18, face = "bold")
  )

gridExtra::grid.arrange(
  plot_1,plot_2,
  ncol=2
)

ggsave(filename = paste0("GB_TE_metagenePlot.svg"),
       path = "/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/MethylationPopStruct_1st2nd_Exp/BISMARKPLOT/merged/TE/",
       device = "svg",
       plot = gridExtra::grid.arrange(
         plot_1,plot_2,
         ncol=2
       ), 
       width = 16, 
       height = 18, 
       dpi = 300)





merged_metagene_TE_CHH<- merged_metagene_TE %>% filter(context=="CHH")
merged_metagene_TE_CG<- merged_metagene_TE %>% filter(context=="CG")
merged_metagene_TE_CHG<- merged_metagene_TE %>% filter(context=="CHG")

wilcox.test(merged_metagene_TE_CHG$AA,
            merged_metagene_TE_CHG$AA37C,
            paired = TRUE)

cohen.d(merged_metagene_TE_CHH$AA37C,
        merged_metagene_TE_CHH$AA,
        paired = TRUE)


effsize::cliff.delta(merged_metagene_TE_CHG$CA37C,
                     merged_metagene_TE_CHG$CA)
merged_metagene_CG
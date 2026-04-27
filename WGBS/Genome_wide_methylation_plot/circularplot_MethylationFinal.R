# Load the packages ----
library(here)
library(karyoploteR)
library(circlize)
library(rtracklayer)
library(vcfR)
library(tidyverse)
library(Rsamtools)
library(gridExtra)
library(readxl)




# Define the potato genome 6.1
soltub_genome_frame <- data.frame(chr = c("chr01", "chr02", "chr03", "chr04", "chr05", "chr06", 
                                          "chr07", "chr08", "chr09", "chr10", "chr11", "chr12"),
                                  start = rep(0, times = 12),
                                  end = c(88591686-1, 46102915-1, 60707570-1, 69236331-1, 55599697-1, 59091578-1,
                                          57639317-1, 59226000-1, 67600300-1, 61044151-1, 46777387-1, 59670755-1))
soltub_genome <- toGRanges(soltub_genome_frame)
# Define the cytobands
soltub_cytobands <- data.frame(chr = rep(c("chr01", "chr02", "chr03", "chr04", "chr05", "chr06", 
                                           "chr07", "chr08", "chr09", "chr10", "chr11", "chr12"), each = 3),
                               start = c(0, 33, 33.3+1e-6,
                                         0, 0, 1e-6, 
                                         0, 11.9, 12.2+1e-6,
                                         0, 26.3, 27.7+1e-6,
                                         0, 27.0, 27.5+1e-6,
                                         0, 15.5, 16.8+1e-6,
                                         0, 16.8, 17.5+1e-6,
                                         0, 18.6, 18.7+1e-6,
                                         0, 18.3, 20.3+1e-6,
                                         0, 24.5, 25.1+1e-6,
                                         0, 21.0, 22.6+1e-6,
                                         0, 18.4, 19.5+1e-6)*1e6,
                               end = c(33*1e6 - 1, 33.3*1e6, 88591686,
                                       0, 0, 46102915,
                                       11.9*1e6 - 1, 12.2*1e6, 60707570,
                                       26.3*1e6 - 1, 27.7*1e6, 69236331,
                                       27.0*1e6 - 1, 27.5*1e6, 55599697,
                                       15.5*1e6 - 1, 16.8*1e6, 59091578,
                                       16.8*1e6 - 1, 17.5*1e6, 57639317,
                                       18.6*1e6 - 1, 18.7*1e6, 59226000,
                                       18.3*1e6 - 1, 20.3*1e6, 67600300,
                                       24.5*1e6 - 1, 25.1*1e6, 61044151,
                                       21.0*1e6 - 1, 22.6*1e6, 46777387,
                                       18.4*1e6 - 1, 19.5*1e6, 59670755),
                               name = c("chr01_pre-cen", "chr01_cen", "chr01_post-cen",
                                        "chr02_pre-cen", "chr02_cen", "chr02_post-cen",
                                        "chr03_pre-cen", "chr03_cen", "chr03_post-cen",
                                        "chr04_pre-cen", "chr04_cen", "chr04_post-cen",
                                        "chr05_pre-cen", "chr05_cen", "chr05_post-cen",
                                        "chr06_pre-cen", "chr06_cen", "chr06_post-cen",
                                        "chr07_pre-cen", "chr07_cen", "chr07_post-cen",
                                        "chr08_pre-cen", "chr08_cen", "chr08_post-cen",
                                        "chr09_pre-cen", "chr09_cen", "chr09_post-cen",
                                        "chr10_pre-cen", "chr10_cen", "chr10_post-cen",
                                        "chr11_pre-cen", "chr11_cen", "chr11_post-cen",
                                        "chr12_pre-cen", "chr12_cen", "chr12_post-cen"),
                               gieStain = rep(c("gpos25", "acen", "gpos25"), times = 12))

##Create the gene track from gene GFF3 files
Genes_gff3 <- import.gff3("/media/rna/Epipotato16TB/00_META_DATA/v6.1_DM_Phytozome/download.20230626.090008/Phytozome/PhytozomeV13/Stuberosum/v6.1/annotation/Stuberosum_686_v6.1.gene_exons.gff3", 
                          version = "3")

Genes_gff3 <- import.gff3("/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/00_META_DATA/Potato_DM_v6.1_reference/download.20230626.090008/Phytozome/PhytozomeV13/Stuberosum/v6.1/annotation/Stuberosum_686_v6.1.gene_exons.gff3", 
                          version = "3")

Genes_gff3 <- Genes_gff3[mcols(Genes_gff3)$type == "gene"]

Genes_gff3_chr <- Genes_gff3[seqnames(Genes_gff3) %in% c("chr01", "chr02", "chr03", "chr04", "chr05", "chr06", 
                                                         "chr07", "chr08", "chr09", "chr10", "chr11", "chr12")]

# Circular karyogram ----
# Create a data frame from the granges object for the karyogram
genes_df <- data.frame(chr = seqnames(Genes_gff3_chr), grange = ranges(Genes_gff3_chr), ID = mcols(Genes_gff3_chr)$ID)


##CREATE THE track from REPEAT ELEMENTD/ transposons gff
TE_GFF<- rtracklayer::import("/media/rna/Epipotato16TB/00_META_DATA/v6.1_DM_Phytozome/DM_1-3_516_R44_potato_genome_assembly.v6.1.hm.out.gff")
TE_GFF<- rtracklayer::import("/media/rna/BAIJI/Darren_Doktorand/Darren_Epipotato/00_META_DATA/Potato_DM_v6.1_reference/DM_1-3_516_R44_potato_genome_assembly.v6.1.hm.out.gff")
TE_GFF_chr <- TE_GFF[seqnames(TE_GFF) %in% c("chr01", "chr02", "chr03", "chr04", "chr05", "chr06", 
                                             "chr07", "chr08", "chr09", "chr10", "chr11", "chr12")]


tile(TE_GFF_chr + 100, width=10)

TE_df <- data.frame(chr = seqnames(TE_GFF_chr), 
                    grange = ranges(TE_GFF_chr), 
                    Target = mcols(TE_GFF_chr)$Target)


Soltu_Chr<- c("chr01", "chr02", "chr03", "chr04", "chr05", "chr06", 
              "chr07", "chr08", "chr09", "chr10", "chr11", "chr12")

soltu_Chr_len <- c(88591686-1, 46102915-1, 60707570-1, 69236331-1, 55599697-1, 59091578-1,
                   57639317-1, 59226000-1, 67600300-1, 61044151-1, 46777387-1, 59670755-1)

names(soltu_Chr_len) <- Soltu_Chr

bin_100kb_chr<- tileGenome(seqlengths = soltu_Chr_len,
                           tilewidth = 100000,
                           cut.last.tile.in.chrom = TRUE)

TE_counts_100kb_bins<- countOverlaps(bin_100kb_chr,TE_GFF_chr)

TE_df_100kb <- data.frame(
  chr = as.character(seqnames(bin_100kb_chr)),
  start = start(bin_100kb_chr) - 1,  # Circos uses 0-based start
  end = end(bin_100kb_chr),
  TE_count = TE_counts_100kb_bins
)



####load DMRs..
AA28_CvsCA28_C<- read.csv("/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/WGBS_HIGH_HEAT/ANALYSIS/CpG_DMRs/FISHER/AA28_CvsCA28_C.txt", sep="")
CHG_AA28_CvsCA28_C <- read.csv("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/CHG_DMRs/FISHER/CHG_AA28_CvsCA28_C.txt", sep="")
CHH_AA28_CvsCA28_C <- read.csv("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/CHH_DMRs/FISHER/CHH_AA28_CvsCA28_C.txt", sep="")

DMR_list <- list(
  CpG = AA28_CvsCA28_C %>% select(seqnames,start,end,width,direction,regionType) %>% mutate(value=1),
  CHG = CHG_AA28_CvsCA28_C %>% select(seqnames,start,end,width,direction,regionType) %>% mutate(value=1),
  CHH = CHH_AA28_CvsCA28_C %>% select(seqnames,start,end,width,direction,regionType) %>% mutate(value=1))


DMR_list_hyperhypo <- list()
for (DMR_Context in names(DMR_list)) {
  
  ##CAMEL WAS USED THE COMPARISON IN THE ANALYSES, BUT I WANT TO MAKE ANNABELLE THE COMPARISON.... hence swwitching the direction
  Hypo <- DMR_list[[DMR_Context]] %>% mutate(regionType=case_when(str_detect(regionType,"gain") ~ "hypomethylated")) %>% filter(regionType == "hypomethylated")
  Hyper <- DMR_list[[DMR_Context]] %>% mutate(regionType=case_when(str_detect(regionType,"loss") ~ "hypermethylated")) %>% filter(regionType == "hypermethylated")
  
  Hypo$value <- 1
  Hyper$value <- 1
  
  DMR_list_hyperhypo[[DMR_Context]]<- list(
    Hypomethylated = Hypo,
    Hypermethylated = Hyper
  )
  
}

path <-"/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/WGBS_HIGH_HEAT/03_BISMARK/CX_FILES/Single_Chr_level_RDS/"
###sinlge sample GW meth 500kb, high heat
Meth_500kb_GW_single<- list(
CG_AA = readRDS(paste0(path,"CG_levels_AA_1.rds")),
CHG_AA = readRDS(paste0(path,"CHG_levels_AA_1.rds")),
CHH_AA = readRDS(paste0(path,"CHH_levels_AA_1.rds")),

CG_CA = readRDS(paste0(path,"CG_levels_CA_1.rds")),
CHG_CA = readRDS(paste0(path,"CHG_levels_CA_1.rds")),
CHH_CA = readRDS(paste0(path,"CHH_levels_CA_1.rds"))
)

##filter for only canonimal chromosomes
Meth_500kb_GW_single_chr <- list()
for (dataframe in names(Meth_500kb_GW_single)) {
  
  Meth_500kb_GW_single_chr[[dataframe]] <- Meth_500kb_GW_single[[dataframe]] %>% dplyr::filter(str_detect(chr,"chr"))
  
}

##ensure only canonical chromosomes left
for (dataframe in names(Meth_500kb_GW_single_chr)) {
  
   print(unique(Meth_500kb_GW_single_chr[[dataframe]]$chr))
  
}
###merged replicates each parent 500kb
Meth_500kb_GW_merged<- list(
  CG_AA = readRDS("/media/rna/Epipotato16TB/Epipotato_methylome/04_BISMARK/METHYLATION_CX_REPORT/Unfiltered/CX_PARENTS_MERGED/CG_AA_chr_level.rds"),
  CHG_AA = readRDS("/media/rna/Epipotato16TB/Epipotato_methylome/04_BISMARK/METHYLATION_CX_REPORT/Unfiltered/CX_PARENTS_MERGED/CHG_AA_chr_level.rds"),
  CHH_AA = readRDS("/media/rna/Epipotato16TB/Epipotato_methylome/04_BISMARK/METHYLATION_CX_REPORT/Unfiltered/CX_PARENTS_MERGED/CHH_AA_chr_level.rds"),
  
  CG_CA = readRDS("/media/rna/Epipotato16TB/Epipotato_methylome/04_BISMARK/METHYLATION_CX_REPORT/Unfiltered/CX_PARENTS_MERGED/CG_CA_chr_level.rds"),
  CHG_CA = readRDS("/media/rna/Epipotato16TB/Epipotato_methylome/04_BISMARK/METHYLATION_CX_REPORT/Unfiltered/CX_PARENTS_MERGED/CHG_CA_chr_level.rds"),
  CHH_CA = readRDS("/media/rna/Epipotato16TB/Epipotato_methylome/04_BISMARK/METHYLATION_CX_REPORT/Unfiltered/CX_PARENTS_MERGED/CHH_CA_chr_level.rds")
)

PrepPlotStartEnd<- function(Meth_DF){
  
  Meth_DF$start <- Meth_DF$chr_pos
  Meth_DF$end <- Meth_DF$chr_pos +500000 -1
  Meth_DF_edit <- Meth_DF %>% dplyr::select(chr,start,end,density)
  
  return(Meth_DF_edit)
  
}

Meth_500kb_GW_single_edit<- lapply(Meth_500kb_GW_single_chr, PrepPlotStartEnd)

circos.genomicTrack(CG_AA_edit, 
                    panel.fun = function(region, value, ...) {
                      circos.genomicLines(region, value, area = TRUE)
                    })
col_fun <- colorRamp2(
  c(min(CG_AA_edit$density),  max(CG_AA_edit$density)),
  c("white","red3")
)


test <- (inner_join(Meth_500kb_GW_single_edit$CG_AA,Meth_500kb_GW_single_edit$CG_CA,by=c("chr","start","end")))

circos.genomicHeatmap(Meth_500kb_GW_single_edit$CG_AA,
                      col= colorRamp2(
                        c(min(Meth_500kb_GW_single_edit$CG_AA$density),  
                          max(Meth_500kb_GW_single_edit$CG_AA$density)),
                        c("white","red3")
                      ), 
                      connection_height = NULL,
                      heatmap_height = 0.07)

# Add track box (border)
circos.track(
  track.index = get.current.track.index(),  # gets the last added track
  panel.fun = function(x, y) {
    circos.rect(
      xleft = CELL_META$xlim[1],
      ybottom = CELL_META$ylim[1],
      xright = CELL_META$xlim[2],
      ytop = CELL_META$ylim[2],
      border = "black",
      col = NA,
      lwd = 1
    )
  },
  bg.border = NA
)

##long reads... ONT
AA_500kb_mean_cov <- read.delim("/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/ONT_READS/05_GA_CORRECTED_MINIMAP2/Primary/AA_500kb_mean_cov.bed", header=FALSE)
CA_500kb_mean_cov <- read.delim("/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/ONT_READS/05_GA_CORRECTED_MINIMAP2/Primary/CA_500kb_mean_cov.bed", header=FALSE)

CA_500kb_mean_cov<- CA_500kb_mean_cov %>% select(V1,V3,V4) %>% rename(
  chr=V1,
  Position=V3,
  Density=V4
)

AA_500kb_mean_cov<- AA_500kb_mean_cov %>% select(V1,V3,V4) %>% rename(
  chr=V1,
  Position=V3,
  Density=V4
)

circos.clear()
circos.par("track.height" = 0.15, start.degree = 80, gap.degree = c(rep(2,11),25), cell.padding=c(0.003,0,0.003,0))
circos.genomicInitialize(soltub_cytobands, 
                         sector.names = c("Chr01\n", "Chr02\n", "Chr03\n", "\nChr04", "\nChr05", "\nChr06",
                                          "\nChr07", "\nChr08", "\nChr09", "Chr10\n", "Chr11\n", "Chr12\n"),
                         labels.cex = 2.3, major.by = 2e7,
                         axis.labels.cex = 0.7*par("cex"))
# circos.genomicInitialize(soltub_cytobands %>% filter(chr == "chr01"), sector.names = c("Chr01"),
#                          labels.cex = 2, major.by = 2e7)
# Karyogram
circos.genomicIdeogram(soltub_cytobands, track.height = 0.03)
# Gene density
circos.genomicDensity(genes_df, col = c("gray20"), track.height = 0.1)
#text(-0.015, 0.75, "Gene density", cex = 1.6)
#TE density
circos.genomicDensity(TE_df, col = c("grey70"), track.height = 0.1)
circos.genomicDensity(DMR_list$CpG, col = c("purple4"), track.height = 0.12)

## long read coverage
circos.genomicHeatmap(AA_500kb_mean_cov, col= colorRamp2(
  c(min(AA_500kb_mean_cov$V4), median(AA_500kb_mean_cov$V4) ,max(AA_500kb_mean_cov$V4)),
  c("white","red3", "red4")
), connection_height = NULL,heatmap_height = 0.07)

circos.genomicHeatmap(CA_500kb_mean_cov, col= colorRamp2(
  c(min(CA_500kb_mean_cov$V4), median(CA_500kb_mean_cov$V4) ,max(CA_500kb_mean_cov$V4)),
  c("white","red3", "red4")
), connection_height = NULL,heatmap_height = 0.07)


DMR_list$CpG$value1 <- 1
circos.genomicTrackPlotRegion(data = DMR_list$CHG %>% select(seqnames,start,end,value), ylim = c(0,1), track.height = 0.1,track.margin = c(0, 0),
                              panel.fun = function(region, value, ...) {
                                circos.genomicLines(region, value, type = "h",
                                                    col = add_transparency("skyblue2", 0.9))
                              })


AddTrackBorder <- function(){
  # Add track box (border)
  circos.track(
    track.index = get.current.track.index(),  # gets the last added track
    panel.fun = function(x, y) {
      circos.rect(
        xleft = CELL_META$xlim[1],
        ybottom = CELL_META$ylim[1],
        xright = CELL_META$xlim[2],
        ytop = CELL_META$ylim[2],
        border = "black",
        col = NA,
        lwd = 1
      )
    },
    bg.border = NA
  )
}

###### plot circos function
PlotCircosMethlation<- function(path,
                                GENE_Track,
                                TE_TRACK,
                                AA_meth_Track,
                                CA_methTrack,
                                DMR_track,
                                Context){
  
  # Open SVG device
  svg_filename <- paste0(path,"Circos_DMR_Highheat_test", Context, ".svg")
  svg(svg_filename, width = 13, height = 13, bg = "white")  # Adjust size as needed
  
  circos.clear()
  circos.par("track.height" = 0.15, start.degree = 80, gap.degree = c(rep(2,11),25), cell.padding=c(0.003,0,0.003,0))
  circos.genomicInitialize(soltub_cytobands, 
                           sector.names = c("Chr01\n", "Chr02\n", "Chr03\n", "\nChr04", "\nChr05", "\nChr06",
                                            "\nChr07", "\nChr08", "\nChr09", "Chr10\n", "Chr11\n", "Chr12\n"),
                           labels.cex = 2.3, major.by = 2e7,
                           axis.labels.cex = 0.7*par("cex"))
  # Karyogram
  circos.genomicIdeogram(soltub_cytobands, track.height = 0.03)
  # Gene density
  circos.genomicDensity(GENE_Track, col = c("gray20"), track.height = 0.1)
  #TE density
  circos.genomicDensity(TE_TRACK, col = c("gray70"), track.height = 0.1)
  
  ##plot annabelle track
  circos.genomicHeatmap(AA_meth_Track, col= colorRamp2(
    c(min(AA_meth_Track$density), max(AA_meth_Track$density)),
    c("white", "red3")
  ), connection_height = NULL,heatmap_height = 0.07)
  
  ##add track border for annabelle meth heatmap
  AddTrackBorder()
  
  ##plot Camel track
  circos.genomicHeatmap(CA_methTrack,col= colorRamp2(
    c(min(CA_methTrack$density), max(CA_methTrack$density)),
    c("white", "red3")
  ), connection_height = NULL,heatmap_height = 0.07)
  
  ##add track border for Camel meth heatmap
  AddTrackBorder()
  
  ##DMR track
  circos.genomicTrackPlotRegion(data = DMR_track %>% select(seqnames,start,end,value), ylim = c(0,1), track.height = 0.07,track.margin = c(0, 0),
                                panel.fun = function(region, value, ...) {
                                  circos.genomicLines(region, value, type = "h",
                                                      col = add_transparency("purple4", 0.9))
                                })
  
  dev.off()
  
}

PATH <- "/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/WGBS_HIGH_HEAT/ANALYSIS/CircosPlotGenome/"

##plot CG circos
PlotCircosMethlation(path = PATH,
                     GENE_Track=genes_df,
                     TE_TRACK = TE_df,
                     AA_meth_Track = Meth_500kb_GW_single_edit$CG_AA,
                     CA_methTrack = Meth_500kb_GW_single_edit$CG_CA,
                     DMR_track = DMR_list$CpG,
                     Context = "CG")

##plot CHG circos
PlotCircosMethlation(path = PATH,
                     GENE_Track=genes_df,
                     TE_TRACK = TE_df,
                     AA_meth_Track = Meth_500kb_GW_single_edit$CHG_AA,
                     CA_methTrack = Meth_500kb_GW_single_edit$CHG_CA,
                     DMR_track = DMR_list$CHG,
                     Context = "CHG")

##plot CHH circos
PlotCircosMethlation(path = PATH,
                     GENE_Track=genes_df,
                     TE_TRACK = TE_df,
                     AA_meth_Track = Meth_500kb_GW_single_edit$CHH_AA,
                     CA_methTrack = Meth_500kb_GW_single_edit$CHH_CA,
                     DMR_track = DMR_list$CHH,
                     Context = "CHH")

###create the legends
library(circlize)
plot.new()
lgd_ = rep(NA, 11)
lgd_[c(0.2,0.95)] = c(0.2,0.95)
legend(x = 0.4, y = 0.5, #coordinates of legend
       legend = lgd_,
       fill = colorRampPalette(colors = c('white','purple4'))(15),
       border = NA,
       y.intersp = 0.5,
       cex = 2, text.font = 2,
       bty = "n")
?legend
library(circlize)
plot.new()
lgd_ = rep(NA, 11)
lgd_[c(0.2,0.95)] = c(0.2,0.95)
legend(x = 0.4, y = 0.5, #coordinates of legend
       legend = lgd_,
       fill = colorRampPalette(colors = c('white','red3'))(15),
       border = NA,
       y.intersp = 0.5,
       cex = 2, text.font = 2,
       bty = "n")

##### Convert to bed

####load DMRs..
AA28_CvsCA28_C<- read.csv("/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/WGBS_HIGH_HEAT/ANALYSIS/CpG_DMRs/FISHER/AA28_CvsCA28_C.txt", sep="")
CHG_AA28_CvsCA28_C <- read.csv("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/CHG_DMRs/FISHER/CHG_AA28_CvsCA28_C.txt", sep="")
CHH_AA28_CvsCA28_C <- read.csv("/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/CHH_DMRs/FISHER/CHH_AA28_CvsCA28_C.txt", sep="")



# Ensure start, end, and seqnames are correctly typed
DMR_list_GR <- list()
for (contexts in names(DMR_list)) {
  
  DMR_list_GR[[contexts]]<- GRanges(
    seqnames = DMR_list[[contexts]]$seqnames,
    ranges = IRanges(start = DMR_list[[contexts]]$start, end = DMR_list[[contexts]]$end),
    strand = DMR_list[[contexts]]$strand,
    score = DMR_list[[contexts]]$direction
  )
  
  
  
}

DMR_list_GR$CpG
export(DMR_list_GR$CpG, "/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/WGBS_HIGH_HEAT/ANALYSIS/CpG_DMRs/FISHER/AA28_CvsCA28_C.bed", format = "BED")
export(DMR_list_GR$CHG, "/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/CHG_DMRs/FISHER/CHG_AA28_CvsCA28_C.bed", format = "BED")
export(DMR_list_GR$CHH, "/media/rna/VAQUITA/DARREN_WGS_HIGH_HEAT/Analysis/CHH_DMRs/FISHER/CHH_AA28_CvsCA28_C.bed", format = "BED")


####plot vcf density
##read vcf file 
vcf <- read.vcfR("/media/rna/POLLUX/Darren_EPIPOTATO_WGS_WGBS/WGBS_HIGH_HEAT/ANALYSIS/filtered_output_parallel_NEWQUAL30.vcf")

# Create a data frame from the fixed data of the VCF file
vcf_fix <- getFIX(vcf)
#vcf_fix <- getFIX(vcf_file)
vcf_df <- as.data.frame(vcf_fix[,1:2])
vcf_df$POS <- as.numeric(vcf_df$POS)
colnames(vcf_df) <- c("chr", "end")
vcf_df$start <- vcf_df$end - 1
vcf_df <- vcf_df[,c(1,3,2)]
vcf_df$value1 <- 1



circos.clear()
circos.par("track.height" = 0.15, start.degree = 80, gap.degree = c(rep(2,11),25), cell.padding=c(0.003,0,0.003,0))
circos.genomicInitialize(soltub_cytobands, 
                         sector.names = c("Chr01\n", "Chr02\n", "Chr03\n", "\nChr04", "\nChr05", "\nChr06",
                                          "\nChr07", "\nChr08", "\nChr09", "Chr10\n", "Chr11\n", "Chr12\n"),
                         labels.cex = 2.3, major.by = 2e7,
                         axis.labels.cex = 0.7*par("cex"))
# Karyogram
circos.genomicIdeogram(soltub_cytobands, track.height = 0.03)
# Gene density
circos.genomicDensity(genes_df, col = c("gray20"), track.height = 0.1)
#TE density
circos.genomicDensity(TE_df, col = c("gray70"), track.height = 0.1)

##DMR track
circos.genomicTrackPlotRegion(data = DMR_list$CpG %>% select(seqnames,start,end,value), ylim = c(0,1), track.height = 0.07,track.margin = c(0, 0),
                              panel.fun = function(region, value, ...) {
                                circos.genomicLines(region, value, type = "h",
                                                    col = add_transparency("purple4", 0.9))
                              })

vcf_df_sampled <- vcf_df %>%
  sample_n(10000)
circos.genomicTrackPlotRegion(data=vcf_df_sampled, ylim = c(0,1), track.height = 0.12,
                              panel.fun = function(region, value, ...) {
                                circos.genomicLines(region, value, type = "h",
                                                    col = add_transparency(col = "black", transparency = 0.9))
                              })

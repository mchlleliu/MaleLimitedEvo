###################################
#
#                             Grieshop et al. 2023
#             DsRed experimental evolution - transcriptomics analysis
#                   SSAV candidates vs other SA candidates
# 
# 
###################################

rm(list=ls()) # Clears the environment
setwd("~/Desktop/UofT/SSAV_RNA/")

# Packages
#########
library(readxl)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggblend)
library(ggpubr)
library(ggstatsplot)
#########

# Get chromosome locations  
##########

all.genes <- read.delim(file="~/Desktop/UofT/SSAV_RNA/Data/all.genes.tsv", sep="\t", header=FALSE)
colnames(all.genes) = c("FlyBaseID")

Xchr <- read.delim(file="~/Desktop/UofT/SSAV_RNA/Data/X.chromosome.genes.tsv", sep="\t", header=TRUE)
colnames(Xchr) = c("FlyBaseID")
Xchr$Chr <- rep("X", dim(Xchr)[1])

Ychr <- read.delim(file="~/Desktop/UofT/SSAV_RNA/Data/Y.chromosome.genes.tsv", sep="\t", header=TRUE)
colnames(Ychr) = c("FlyBaseID")
Ychr$Chr <- rep("Y", dim(Ychr)[1])


chr2L <- read.delim(file="~/Desktop/UofT/SSAV_RNA/Data/2L.chromosome.genes.tsv", sep="\t", header=FALSE)
colnames(chr2L) = c("FlyBaseID")
chr2L$Chr <- rep("2L", dim(chr2L)[1])
chr2R <- read.delim(file="~/Desktop/UofT/SSAV_RNA/Data/2R.chromosome.genes.tsv", sep="\t", header=FALSE)
colnames(chr2R) = c("FlyBaseID")
chr2R$Chr <- rep("2R", dim(chr2R)[1])
# 
chr2 <- rbind(chr2L, chr2R)
chr2$Chr <- rep("2", dim(chr2)[1])

chr3L <- read.delim(file="~/Desktop/UofT/SSAV_RNA/Data/3L.chromosome.genes.tsv", sep="\t", header=FALSE)
colnames(chr3L) = c("FlyBaseID")
chr3L$Chr <- rep("3L", dim(chr3L)[1])
chr3R <- read.delim(file="~/Desktop/UofT/SSAV_RNA/Data/3R.chromosome.genes.tsv", sep="\t", header=FALSE)
colnames(chr3R) = c("FlyBaseID")
chr3R$Chr <- rep("3R", dim(chr3R)[1])
#
chr3 <- rbind(chr3L, chr3R)
chr3$Chr <- rep("3", dim(chr3)[1])

Chrs <- rbind(Xchr, Ychr, chr2, chr3) # Not all genes; just X, Y, 2, and 3.
Chrs$Chr <- as.factor(Chrs$Chr)

Chrs_All <- rbind(Xchr, Ychr, chr2L, chr2R, chr3L, chr3R)
Chrs_All$Chr <- as.factor(Chrs_All$Chr)
##########

# Prepare dataset
########
# load results if not loaded in env.
SSAV.geno <- read.delim("Results/All.geno_candidates.tsv", header = T, sep = "\t")

SSAV.geno <- merge(SSAV.geno, Chrs, by = "FlyBaseID", all = TRUE)
SSAV.geno <- SSAV.geno[!is.na(SSAV.geno$Sig) & !is.na(SSAV.geno$Chr),]

dim(SSAV.geno[!is.na(SSAV.geno$Sig),])

# AS candidates
jseq.All.geno <- read.delim(file="Results/jseq.All.geno.txt", sep = "\t", header = T)
jseq.All.geno <- merge(jseq.All.geno, Chrs, by = "FlyBaseID") 


# combined all candidates
test_df <- merge(jseq.All.geno[,c("FlyBaseID", "Sig")], 
                 SSAV.geno[,c("FlyBaseID", "Sig")], 
                 by ="FlyBaseID", all = T)
str(test_df)
test_df <- test_df[!(is.na(test_df$Sig.x) & is.na(test_df$Sig.y)),]
test_df <- test_df %>%
  dplyr::mutate(Sig = ifelse((!is.na(Sig.x) & Sig.x), TRUE,
                             ifelse(!is.na(Sig.y) & Sig.y, TRUE, FALSE)))
colnames(test_df) <- c("FlyBaseID", "Sig.DS", "Sig.DE", "Sig")
str(test_df)
#########


TestEnrichment <- function(data, df1 = "Sig", df2, xlab = "DE", ylab){
  test <- fisher.test(x = data[[df1]], y = data[[df2]])
  print(test)
  
  # kinda impractical but can't figure out how to 
  #     get ggbarstats to work with strings... alas
  colnames(data)[colnames(data) == paste0(df1)] = "df1"
  colnames(data)[colnames(data) == paste0(df2)] = "df2"
  
  mos_plot <- ggbarstats(
    data, df2, df1,
    results.subtitle = FALSE, # include fisher's exact test results
    subtitle = paste0(
      "Fisher's exact test", ", p-value = ",
      ifelse(test$p.value < 0.001, "< 0.001", round(test$p.value, 3))
    ),
    label = "both", # use both percentage and count as labels, 
    perc.k = 2, # include 2 decimal points in perc
    xlab = NULL
  ) +
    # some plot and theme settings
    scale_fill_manual(labels = c(paste0(ylab), paste0("not ",ylab)),
                      values = c("purple3", "darkgrey")) + # "Chr-2", "Chr-3", "X-Chr"
    scale_x_discrete(labels = c(paste0("not",xlab), paste0(xlab))) +
    theme(plot.title.position = c("panel"),
          legend.title = element_blank(),
          legend.position = c("bottom"),
          legend.box.background = element_rect(),
          legend.text = element_text(size = 20, color = "black"),
          axis.text.x = element_text(size=20, margin = margin(5,0,0,0), color = "black"),
          axis.text.y = element_text(size=10, margin = margin(0,5,0,0), color = "black"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=40, margin = margin(0,10,0,0), color = "black"),
          plot.title = element_text(size=40, margin = margin(0,0,0,0), color = "black"),
          plot.margin = margin(6,6,6,6)
    )
  
  # get the contingency table
  plot.data <- mos_plot$data[,1:4]
  # rename columns appropriately
  colnames(plot.data)[colnames(plot.data) == "df1"] = paste0(xlab)
  colnames(plot.data)[colnames(plot.data) == "df2"] = paste0(df2)
  
  # assign new object in parent env
  assign(paste0(xlab,".vs.",ylab,".data"), plot.data, envir = parent.frame())
  
  # return the plot
  return(mos_plot)
}



# Innocenti and Morrow (2010)
##########
# Get Innocenti & Morrow (2010) SA candidates
InnoMorrow_SBGE <- read.delim("Data/InnoMorrow.SBGE.tsv")
str(InnoMorrow_SBGE)

# combine InnoMorrow status with SSAV dataset
SSAV.geno <- SSAV.geno %>%
  dplyr::mutate(IsInnoMorr = FlyBaseID %in% InnoMorrow_SBGE$FlyBaseID[InnoMorrow_SBGE$Sig])
jseq.All.geno <- jseq.All.geno %>% 
  dplyr::mutate(IsInnoMorr = FlyBaseID %in% InnoMorrow_SBGE$FlyBaseID[InnoMorrow_SBGE$Sig])
test_df <- test_df %>% 
  mutate(IsInnoMorr = FlyBaseID %in% InnoMorrow_SBGE$FlyBaseID[InnoMorrow_SBGE$Sig])
rm(InnoMorrow_SBGE)

# Only DE data
DE.vs.InnoMorr <- TestEnrichment(SSAV.geno, df2 = "IsInnoMorr", xlab = "DE", ylab = "InnoMorr")
DE.vs.InnoMorr.data

# Only DS data
DS.vs.InnoMorr <- TestEnrichment(jseq.All.geno, df2 = "IsInnoMorr", xlab = "DS", ylab = "InnoMorr")
DS.vs.InnoMorr.data

# all DE and DS results
SA.vs.InnoMorr <- TestEnrichment(test_df, df2 = "IsInnoMorr", xlab = "SA", ylab = "InnoMorr")
SA.vs.InnoMorr.data

##########


# Ruzicka et al. (2019)
##########
# Get Ruzicka et al. (2019) SA candidates
Ruzicka <- read_excel("~/Desktop/UofT/SSAV_RNA/Data/Ruzicka_SA_genes.xlsx")
colnames(Ruzicka) <- c("FlyBaseID", "FlyBaseName", "Chr", "SA_miss", 
                       "SA_nonmiss", "NSA_miss", "NSA_nonmiss")

SSAV.geno <- SSAV.geno %>% 
  mutate(IsRuz = FlyBaseID %in% Ruzicka$FlyBaseID)
jseq.All.geno <- jseq.All.geno %>% 
  mutate(IsRuz = FlyBaseID %in% Ruzicka$FlyBaseID)
test_df <- test_df %>% 
  mutate(IsRuz = FlyBaseID %in% Ruzicka$FlyBaseID)

rm(Ruzicka)

# Only DE data
DE.vs.Ruz <- TestEnrichment(SSAV.geno, df2 = "IsRuz", xlab = "DE", ylab = "Ruz")
DE.vs.Ruz.data

# Only DS data
DS.vs.Ruz <- TestEnrichment(jseq.All.geno, df2 = "IsRuz", xlab = "DS", ylab = "Ruz")
DS.vs.Ruz.data

# all DE and DS results
SA.vs.Ruz <- TestEnrichment(test_df, df2 = "IsRuz", xlab = "SA", ylab = "Ruz")
SA.vs.Ruz.data
##########



# Wong and Holman (2023)
# MESSY. NOT CLEANED BCS NOT INCLUDED IN FINAL RESULT
#########
# Get Wong & Holman TWAS SA candidate genes
Wong_Holman <- read.csv(file="Data/Wong_Holman_TWAS.csv", header = TRUE)
colnames(Wong_Holman)[1] <- "FlyBaseID"
Wong_Holman <- Wong_Holman %>% mutate(SA = ifelse((Female.early.effect < 0 & Female.late.effect < 0) &
                                                    (Male.early.effect > 0 & Male.late.effect > 0), TRUE, 
                                                  ifelse((Female.early.effect > 0 & Female.late.effect > 0) &
                                                           (Male.early.effect < 0 & Male.late.effect < 0), TRUE, FALSE)))

SSAV.geno <- SSAV.geno %>% 
  mutate(IsWong = FlyBaseID %in% Wong_Holman[Wong_Holman$SA,]$FlyBaseID)
jseq.All.geno.tmp <- jseq.All.geno.tmp %>% 
  mutate(IsWong = FlyBaseID %in% Wong_Holman[Wong_Holman$SA,]$FlyBaseID)

SSAV.geno_Chr2 <- SSAV.geno_Chr2 %>% 
  mutate(IsWong = FlyBaseID %in% Wong_Holman[Wong_Holman$SA,]$FlyBaseID)
jseq.All.geno.tmp_Chr2 <- jseq.All.geno.tmp_Chr2 %>% 
  mutate(IsWong = FlyBaseID %in% Wong_Holman[Wong_Holman$SA,]$FlyBaseID)

test <- fisher.test(x = SSAV.geno$Sig, y = SSAV.geno$IsWong)

mos_plot_Wong
mos_plot_Wong <- ggbarstats(
  SSAV.geno, IsWong, Sig,
  results.subtitle = FALSE,
  subtitle = paste0(
    "Fisher's exact test", ", p-value = ",
    ifelse(test$p.value < 0.001, "< 0.001", round(test$p.value, 3))
  ), 
  xlab = NULL
) +
  scale_fill_manual(labels = c("in_Wong", "not_in_Wong"),
                    values = c("darkorchid4", "darkgrey")) + # "Chr-2", "Chr-3", "X-Chr"
  scale_x_discrete(labels = c("Background", "Candidates")) +
  theme(plot.title.position = c("panel"),
        legend.title = element_blank(),
        legend.position = c("bottom"),
        #legend.justification = c("right", "bottom"),
        #legend.box.just = "left",
        #legend.box.background = element_rect(),
        legend.box.background = element_rect(),
        #legend.box.margin = margin(4, 6, 6, 6),
        legend.text = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size=20, margin = margin(5,0,0,0), color = "black"),
        axis.text.y = element_text(size=10, margin = margin(0,5,0,0), color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=40, margin = margin(0,10,0,0), color = "black"),
        plot.title = element_text(size=40, margin = margin(0,0,0,0), color = "black"),
        plot.margin = margin(6,6,6,6)
  )


#########




# Chloe's genomics candidates
# MESSY. NOT CLEANED BCS NOT INCLUDED IN FINAL RESULT
##########
# # Get SSAV candidates from genomics data
# SSAV_chloe <- read.delim("Data/SSAV_genomics_candidates_list.txt", header = TRUE)

# New candidates
SSAV_chloe <- read.csv("Data/SSAV_cand_NEW.csv", header = TRUE)
SSAV_chloe <- SSAV_chloe[SSAV_chloe$new.FDR5 & !is.na(SSAV_chloe$geneID),] %>%
  dplyr::distinct(geneID, .keep_all = TRUE) %>%
  dplyr::rename(FlyBaseID = geneID)
dim(SSAV_chloe) 

# keep only genes in both genomics and RNA-seq data
SSAV.geno <- SSAV.geno %>% 
  mutate(IsChloe = FlyBaseID %in% SSAV_chloe[SSAV_chloe$new.FDR5,]$FlyBaseID)
# jseq.All.geno.tmp <- jseq.All.geno.tmp %>% 
#   mutate(IsChloe = FlyBaseID %in% SSAV_chloe[SSAV_chloe$new.FDR5,]$FlyBaseID)

SSAV.geno_Chr2 <- SSAV.geno_Chr2 %>% 
  mutate(IsChloe = FlyBaseID %in% SSAV_chloe[SSAV_chloe$new.FDR5,]$FlyBaseID)
# jseq.All.geno.tmp_Chr2 <- jseq.All.geno.tmp_Chr2 %>% 
#   mutate(IsChloe = FlyBaseID %in% SSAV_chloe[SSAV_chloe$new.FDR5,]$FlyBaseID)


# only genes in both Genomics & Transcriptomics
SSAV_chloe_ALL <- read.csv("Data/SSAV_cand_NEW.csv", header = TRUE)
SSAV_chloe_ALL <- unique(SSAV_chloe_ALL$geneID)
SSAV_geno_trans <- SSAV.geno[SSAV.geno$FlyBaseID %in% SSAV_chloe_ALL,]
dim(SSAV_geno_trans)
SSAV_geno_trans_Chr2 <- SSAV_geno_trans[SSAV_geno_trans$Chr == "2",]
dim(SSAV_geno_trans_Chr2)

test <- fisher.test(x = SSAV_geno_trans$Sig, y = SSAV_geno_trans$IsChloe)

mos_plot_Chloe # the genomics candidates only consists of Chr2 genes. So filtering out the data below...


mos_plot_Chloe <- ggbarstats(
  SSAV_geno_trans[SSAV_geno_trans$Chr == "2",], IsChloe, Sig,
  results.subtitle = FALSE,
  subtitle = paste0(
    "Fisher's exact test", ", p-value = ",
    ifelse(test$p.value < 0.001, "< 0.001", round(test$p.value, 3))
  ), 
  xlab = NULL
) +
  scale_fill_manual(labels = c("in_SSAV_geno", "not_in_SSAV_geno"),
                    values = c("darkorchid4", "darkgrey")) + # "Chr-2", "Chr-3", "X-Chr"
  scale_x_discrete(labels = c("Background", "Candidates")) +
  theme(plot.title.position = c("panel"),
        legend.title = element_blank(),
        legend.position = c("bottom"),
        #legend.justification = c("right", "bottom"),
        #legend.box.just = "left",
        #legend.box.background = element_rect(),
        legend.box.background = element_rect(),
        #legend.box.margin = margin(4, 6, 6, 6),
        legend.text = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size=20, margin = margin(5,0,0,0), color = "black"),
        axis.text.y = element_text(size=10, margin = margin(0,5,0,0), color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=40, margin = margin(0,10,0,0), color = "black"),
        plot.title = element_text(size=40, margin = margin(0,0,0,0), color = "black"),
        plot.margin = margin(6,6,6,6)
  )

##########


pdf(file = "~/Desktop/UofT/SSAV_RNA/Plots/Enrichment_Tests/Chloe_Chr2.pdf",  # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 7) # The height of the plot in inches
mos_plot_Chloe_ovl2 # mos_plot_InnoMorr_Chr2, mos_plot_Ruz, mos_plot_Ruz_Chr2, mos_plot_Chloe_Chr2, mos_plot_Wong, mos_plot_Wong_Chr2
dev.off()


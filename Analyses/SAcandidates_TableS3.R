###################################
#
#                             Grieshop et al. 2025
#                             Author: Michelle Liu
#             DsRed experimental evolution - transcriptomics analysis
#                   SSAV candidates vs other SA candidates
#                                   Table S3
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
source("Chromosome_df.R")

# Prepare dataset
########
# load results if not loaded in env.
SSAV.geno <- read.delim("Results/All.geno_DE.candidates.tsv", header = T, sep = "\t")

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

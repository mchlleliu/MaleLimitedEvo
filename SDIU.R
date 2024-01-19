###################################
#
#                             Grieshop et al. 2023
#                             Author: Michelle Liu
#             DsRed experimental evolution - transcriptomics analysis
#                             Sex-specific splicing
# 
# 
###################################

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

# Load Singh & Agrawal 2023 dataset
#######
# load Singh & Agrawal 2023 dataset
SDIU <- read.csv(file="~/Desktop/UofT/SSAV_RNA/Data/SBGEandSSSdataForMBE.csv", sep=",", header=TRUE)
colnames(SDIU)[2] <- "FlyBaseID"
# change to order for plotting
SDIU <- mutate(SDIU, SBGEcat.body.Osada = case_when(
  SBGEcat.body.Osada == "extFB"   ~ "a.ext.fbg",
  SBGEcat.body.Osada == "sFB"  ~ "b.more.fbg",
  SBGEcat.body.Osada == "FB" ~ "c.fbg",
  SBGEcat.body.Osada == "UB" ~ "d.ubg",
  SBGEcat.body.Osada == "MB" ~ "e.mbg",
  SBGEcat.body.Osada == "sMB" ~ "f.more.mbg",
  SBGEcat.body.Osada == "extMB" ~ "g.ext.mbg",
  TRUE              ~ SBGEcat.body.Osada  # Keep other values unchanged
))

# change to order for plotting
SDIU <- mutate(SDIU, SBGEcat.head.Osada = case_when(
  SBGEcat.head.Osada == "extFB"   ~ "a.ext.fbg",
  SBGEcat.head.Osada == "sFB"  ~ "b.more.fbg",
  SBGEcat.head.Osada == "FB" ~ "c.fbg",
  SBGEcat.head.Osada == "UB" ~ "d.ubg",
  SBGEcat.head.Osada == "MB" ~ "e.mbg",
  SBGEcat.head.Osada == "sMB" ~ "f.more.mbg",
  SBGEcat.head.Osada == "extMB" ~ "g.ext.mbg",
  TRUE              ~ SBGEcat.head.Osada  # Keep other values unchanged
))
str(SDIU)
#######

# Prepare dataset
#######
# load results if not loaded in env.
A.f.geno <- read.delim("Results/A.f.geno_candidates.tsv")
A.m.geno <- read.delim("Results/A.m.geno_candidates.tsv")

# list of all genes that could be quantified in either or both SSAV male and female samples
all.genes <- merge(A.f.geno, A.m.geno, by = "FlyBaseID", all = T)
all.genes <- data.frame(FlyBaseID = unlist(all.genes[, "FlyBaseID"]))


#######


## combine SDIU status with SSAV dataset
#######
# all_genes_SDIU <- tibble(all.genes)
# all_genes_SDIU <- all_genes_SDIU %>% 
#   mutate(Sig.Af = FlyBaseID %in% A.f.geno[A.f.geno$Sig,]$FlyBaseID,
#          Sig.Am = FlyBaseID %in% A.m.geno[A.m.geno$Sig,]$FlyBaseID,
#          Sig = Sig.Af | Sig.Am,
#          SDIU.body = FlyBaseID %in% SDIU[!is.na(SDIU$SDIU.body.sig) & 
#                                             SDIU$SDIU.body.sig,]$FlyBaseID,
#          SDIU.head = FlyBaseID %in% SDIU[!is.na(SDIU$SDIU.head.sig) & 
#                                            SDIU$SDIU.head.sig,]$FlyBaseID)
# 
# # only Chr 2 genes
# all_genes_SDIU_Chr2 <- merge(all_genes_SDIU, Chrs, by = "FlyBaseID", all = T)
# all_genes_SDIU_Chr2 <- all_genes_SDIU_Chr2[!is.na(all_genes_SDIU_Chr2$Sig) &
#                                              all_genes_SDIU_Chr2$Chr == "2",]

## another approach to perform the analysis with genes in both studies
SSAV.geno <- merge(A.m.geno, A.f.geno, by = "FlyBaseID", all = TRUE)
colnames(SSAV.geno) <- c("FlyBaseID", "A.m.exp_geno", "A.m.se_geno", "A.m.padj", "A.m.TopSig", "A.m.Sig",
                         "A.f.exp_geno", "A.f.se_geno", "A.f.padj", "A.f.Sig")
# column denotes genes that are candidates in males or females
SSAV.geno <- SSAV.geno %>% mutate(Sig = ifelse(!is.na(A.m.Sig) & A.m.Sig, TRUE, 
                                               ifelse(!is.na(A.f.Sig) & A.f.Sig, TRUE, FALSE))) 
SSAV.geno <- merge(SSAV.geno, SDIU, by = "FlyBaseID", all = TRUE)
SSAV.geno <- SSAV.geno[!is.na(SSAV.geno$Sig) & 
                         (!is.na(SSAV.geno$SDIU.body.sig) | !is.na(SSAV.geno$SDIU.head.sig)),]
SSAV.geno_Chr2 <- SSAV.geno[SSAV.geno$chrm == "2L" | 
                              SSAV.geno$chrm == "2R",]
######



#######

test <- fisher.test(x = SSAV.geno_Chr2$Sig, 
                    y = SSAV.geno_Chr2$SDIU.head.sig)

mos_plot_SDIU <- ggbarstats(
  SSAV.geno, SDIU.head.sig, Sig,
  results.subtitle = FALSE,
  subtitle = paste0(
    "Fisher's exact test", ", p-value = ",
    ifelse(test$p.value < 0.001, "< 0.001", round(test$p.value, 3))
  ), 
  xlab = NULL
) +
  scale_fill_manual(labels = c("SDIU", "not_SDIU"),
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




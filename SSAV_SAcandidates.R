###################################
#
#                             Grieshop et al. 2023
#             DsRed experimental evolution - transcriptomics analysis
#                   SSAV candidates vs other SA candidates
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

# Prepare dataset
#########
# load results if not loaded in env.
A.f.geno <- read.delim("Results/A.f.geno_candidates.tsv")
A.m.geno <- read.delim("Results/A.m.geno_candidates.tsv")

# list of all genes that could be quantified in either or both SSAV male and female samples
all.genes <- merge(A.f.geno, A.m.geno, by = "FlyBaseID", all = T)
all.genes <- data.frame(FlyBaseID = unlist(all.genes[, "FlyBaseID"]))

SSAV.geno <- merge(A.m.geno, A.f.geno, by = "FlyBaseID", all = TRUE)
colnames(SSAV.geno) <- c("FlyBaseID", "A.m.exp_geno", "A.m.se_geno", "A.m.padj", "A.m.TopSig", "A.m.Sig",
                         "A.f.exp_geno", "A.f.se_geno", "A.f.padj", "A.f.Sig")
# column denotes genes that are candidates in males or females
SSAV.geno <- SSAV.geno %>% mutate(Sig = ifelse(!is.na(A.m.Sig) & A.m.Sig, TRUE, 
                                               ifelse(!is.na(A.f.Sig) & A.f.Sig, TRUE, FALSE))) 
SSAV.geno <- merge(SSAV.geno, Chrs, by = "FlyBaseID", all = TRUE)
SSAV.geno <- SSAV.geno[!is.na(SSAV.geno$Sig) & !is.na(SSAV.geno$Chr),]
SSAV.geno_Chr2 <- SSAV.geno[SSAV.geno$Chr == "2",]

jseq.All.geno.tmp <- read.delim("Results/jseq.All.geno.txt", sep = "\t", header=T)
jseq.All.geno.tmp <-merge(jseq.All.geno.tmp, Chrs, by = "FlyBaseID", all = TRUE)
# jseq.All.geno.tmp <- jseq.All.geno.tmp[!is.na(jseq.All.geno.tmp$Sig) & !is.na(jseq.All.geno.tmp$Chrs),]
jseq.All.geno.tmp_Chr2 <- jseq.All.geno.tmp[jseq.All.geno.tmp$Chr=="2",]
#########


# Innocenti and Morrow (2010)
##########
# Get Innocenti & Morrow (2010) SA candidates
Inno_Morrow <- read.delim(file="~/Desktop/UofT/SSAV_RNA/Data/FlyBaseID_InnoMorrow.txt", sep="\t", header=TRUE)
colnames(Inno_Morrow) <- c("Symbol", "FlyBaseID", "Related_Record")

# some of the Symbols are associated with the same FlyBaseID gene. In that case, I just retained the FlyBaseID info
# i.e., these features will be regarded as one.
Innocenti_Morrow_SA_genes <- read_excel("~/Desktop/UofT/SSAV_RNA/Data/Innocenti_Morrow_SA_genes.xls", sheet = "Antagonistic genes")
Innocenti_Morrow_SA_genes <- merge(Innocenti_Morrow_SA_genes, Inno_Morrow, by = "Symbol", all = T)

# combine InnoMorrow status with SSAV dataset
SSAV.geno <- SSAV.geno %>% 
  mutate(IsInnoMorr = FlyBaseID %in% Innocenti_Morrow_SA_genes$FlyBaseID)
jseq.All.geno.tmp <- jseq.All.geno.tmp %>% 
  mutate(IsInnoMorr = FlyBaseID %in% Innocenti_Morrow_SA_genes$FlyBaseID)

# only Chr 2 genes
SSAV.geno_Chr2 <- SSAV.geno_Chr2 %>% 
  mutate(IsInnoMorr = FlyBaseID %in% Innocenti_Morrow_SA_genes$FlyBaseID)
jseq.All.geno.tmp_Chr2 <- jseq.All.geno.tmp_Chr2 %>% 
  mutate(IsInnoMorr = FlyBaseID %in% Innocenti_Morrow_SA_genes$FlyBaseID)

test <- fisher.test(x = jseq.All.geno.tmp_Chr2$Sig, y = jseq.All.geno.tmp_Chr2$IsInnoMorr)

# enrichment for all genes in SSAV males and females, vs Innocenti & Morrow all genes
mos_plot_InnoMorr

# enrichment for Chr2 genes in SSAV males and females, vs Innocenti & Morrow Chr2 genes
mos_plot_InnoMorr_Chr2 <- ggbarstats(
  SSAV.geno_Chr2, IsInnoMorr, Sig,
  results.subtitle = FALSE,
  subtitle = paste0(
    "Fisher's exact test", ", p-value = ",
    ifelse(test$p.value < 0.001, "< 0.001", round(test$p.value, 3))
  ), 
  xlab = NULL
) +
  scale_fill_manual(labels = c("SA_InnoMorr", "notSA_InnoMorr"),
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


# Ruzicka et al. (2019)
##########
# Get Ruzicka et al. (2019) SA candidates
Ruzicka <- read_excel("~/Desktop/UofT/SSAV_RNA/Data/Ruzicka_SA_genes.xlsx")
colnames(Ruzicka) <- c("FlyBaseID", "FlyBaseName", "Chr", "SA_miss", 
                       "SA_nonmiss", "NSA_miss", "NSA_nonmiss")


SSAV.geno <- SSAV.geno %>% 
  mutate(IsRuz = FlyBaseID %in% Ruzicka$FlyBaseID)
jseq.All.geno.tmp <- jseq.All.geno.tmp %>%
  mutate(IsRuz = FlyBaseID %in% Ruzicka$FlyBaseID)

# no difference if only looking at Chr2 either.
SSAV.geno_Chr2 <- SSAV.geno_Chr2 %>% 
  mutate(IsRuz = FlyBaseID %in% Ruzicka$FlyBaseID)
jseq.All.geno.tmp_Chr2 <- jseq.All.geno.tmp_Chr2 %>% 
  mutate(IsRuz = FlyBaseID %in% Ruzicka$FlyBaseID)

test <- fisher.test(x = jseq.All.geno.tmp_Chr2$Sig, y = jseq.All.geno.tmp_Chr2$IsRuz)

# enrichment for all genes in SSAV males and females, vs Ruzicka all genes
mos_plot_Ruz

# enrichment for Chr 2 genes in SSAV males and females, vs Ruzicka Chr 2 genes
mos_plot_Ruz_Chr2 <- ggbarstats(
  jseq.All.geno.tmp, IsRuz, Sig,
  results.subtitle = FALSE,
  subtitle = paste0(
    "Fisher's exact test", ", p-value = ",
    ifelse(test$p.value < 0.001, "< 0.001", round(test$p.value, 3))
  ), 
  xlab = NULL
) +
  scale_fill_manual(labels = c("in_Ruz", "not_in_Ruz"),
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



# Wong and Holman (2023)
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

test <- fisher.test(x = jseq.All.geno.tmp_Chr2$Sig, y = jseq.All.geno.tmp_Chr2$IsWong)

mos_plot_Wong
mos_plot_Wong_Chr2 <- ggbarstats(
  SSAV.geno_Chr2, IsWong, Sig,
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
##########
# Get SSAV candidates from genomics data
SSAV_chloe <- read.delim("Data/SSAV_genomics_candidates_list.txt", header = TRUE)

# New candidates
SSAV_chloe <- read.csv("Data/SSAV_cand_NEW.csv", header = TRUE)
SSAV_chloe <- SSAV_chloe[SSAV_chloe$new.FDR5 & !is.na(SSAV_chloe$geneID),] %>%
  dplyr::distinct(geneID, .keep_all = TRUE) %>%
  dplyr::rename(FlyBaseID = geneID)
dim(SSAV_chloe) 

# keep only genes in both genomics and RNA-seq data
SSAV.geno <- SSAV.geno %>% 
  mutate(IsChloe = FlyBaseID %in% SSAV_chloe[SSAV_chloe$new.FDR5,]$FlyBaseID)
jseq.All.geno.tmp <- jseq.All.geno.tmp %>% 
  mutate(IsChloe = FlyBaseID %in% SSAV_chloe[SSAV_chloe$new.FDR5,]$FlyBaseID)

SSAV.geno_Chr2 <- SSAV.geno_Chr2 %>% 
  mutate(IsChloe = FlyBaseID %in% SSAV_chloe[SSAV_chloe$new.FDR5,]$FlyBaseID)
jseq.All.geno.tmp_Chr2 <- jseq.All.geno.tmp_Chr2 %>% 
  mutate(IsChloe = FlyBaseID %in% SSAV_chloe[SSAV_chloe$new.FDR5,]$FlyBaseID)


# only genes in both Genomics & Transcriptomics
SSAV_chloe_ALL <- read.csv("Data/SSAV_cand_NEW.csv", header = TRUE)
SSAV_chloe_ALL <- unique(SSAV_chloe_ALL$geneID)
SSAV_geno_trans <- SSAV.geno[SSAV.geno$FlyBaseID %in% SSAV_chloe_ALL,]
dim(SSAV_geno_trans)
SSAV_geno_trans_Chr2 <- SSAV_geno_trans[SSAV_geno_trans$Chr == "2",]
dim(SSAV_geno_trans_Chr2)

test <- fisher.test(x = jseq.All.geno.tmp_Chr2$Sig, y = jseq.All.geno.tmp_Chr2$IsChloe)

mos_plot_Chloe # the genomics candidates only consists of Chr2 genes. So filtering out the data below...


mos_plot_Chloe_ovl2 <- ggbarstats(
  jseq.All.geno.tmp_Chr2, IsChloe, Sig,
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


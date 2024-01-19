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
all_genes_InnoMorrow <- tibble(all.genes)
all_genes_InnoMorrow <- all_genes_InnoMorrow %>% 
  mutate(Sig.Af = FlyBaseID %in% A.f.geno[A.f.geno$Sig,]$FlyBaseID,
         Sig.Am = FlyBaseID %in% A.m.geno[A.m.geno$Sig,]$FlyBaseID,
         Sig = Sig.Af | Sig.Am,
         IsInnoMorr = FlyBaseID %in% Innocenti_Morrow_SA_genes$FlyBaseID)

# only Chr 2 genes
all_genes_InnoMorrow_Chr2 <- merge(all_genes_InnoMorrow, Chrs, by = "FlyBaseID", all = T)
all_genes_InnoMorrow_Chr2 <- all_genes_InnoMorrow_Chr2[!is.na(all_genes_InnoMorrow_Chr2$Sig) &
                                                         all_genes_InnoMorrow_Chr2$Chr == "2",]

test <- fisher.test(x = all_genes_InnoMorrow$Sig, y = all_genes_InnoMorrow$IsInnoMorr)

# enrichment for all genes in SSAV males and females, vs Innocenti & Morrow all genes
mos_plot_InnoMorr

# enrichment for Chr2 genes in SSAV males and females, vs Innocenti & Morrow Chr2 genes
mos_plot_InnoMorr_Chr2 <- ggbarstats(
  all_genes_InnoMorrow_Chr2, IsInnoMorr, Sig,
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


all_genes_Ruz <- tibble(all.genes)
all_genes_Ruz <- all_genes_Ruz %>% 
  mutate(Sig.Af = FlyBaseID %in% A.f.geno[A.f.geno$Sig,]$FlyBaseID,
         Sig.Am = FlyBaseID %in% A.m.geno[A.m.geno$Sig,]$FlyBaseID,
         Sig = Sig.Af | Sig.Am,
         IsRuz = FlyBaseID %in% Ruzicka$FlyBaseID)

# no difference if only looking at Chr2 either.
all_genes_Ruz_Chr2 <- merge(all_genes_Ruz, Chrs, by = "FlyBaseID", all = T)
all_genes_Ruz_Chr2 <- all_genes_Ruz_Chr2[!is.na(all_genes_Ruz_Chr2$Sig) & 
                                           all_genes_Ruz_Chr2$Chr == "2",]

test <- fisher.test(x = all_genes_Ruz_Chr2$Sig, y = all_genes_Ruz_Chr2$IsRuz)

# enrichment for all genes in SSAV males and females, vs Ruzicka all genes
mos_plot_Ruz

# enrichment for Chr 2 genes in SSAV males and females, vs Ruzicka Chr 2 genes
mos_plot_Ruz_Chr2 <- ggbarstats(
  all_genes_Ruz_Chr2[!is.na(all_genes_Ruz_Chr2$Sig) & 
                       !is.na(all_genes_Ruz_Chr2$IsRuz),], IsRuz, Sig,
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

# function modified from Wong & Holman suppl. mat
# under Plotting and modelling the evidence for antagonism
# https://lukeholman.github.io/fitnessGWAS/plot_models_variant_effects.html#Frequencies_of_antagonistic_transcripts
# don't really have the correct data for this...
get_antagonism_ratios <- function(dat){
  dat %>%
    
    # Convert the LFSR to the probability that the effect size is positive
    mutate(pp_female_early = ifelse(Female.early.effect > 0, Female.early.pval, 1 - Female.early.pval),
           pp_female_late  = ifelse(Female.late.effect > 0, Female.late.pval, 1 - Female.late.pval),
           pp_male_early   = ifelse(Male.early.effect > 0, Male.early.pval, 1 - Male.early.pval),
           pp_male_late    = ifelse(Male.late.effect > 0, Male.late.pval, 1 - Male.late.pval)) %>%
    
    # Calculate the probabilities that beta_i and beta_j have the same/opposite signs
    mutate(p_sex_concord_early = pp_female_early * pp_male_early + 
             (1 - pp_female_early) * (1 - pp_male_early),
           p_sex_antag_early = pp_female_early * (1 - pp_male_early) + 
             (1 - pp_female_early) * pp_male_early,
           p_sex_concord_late  = pp_female_late * pp_male_late + 
             (1 - pp_female_late) * (1 - pp_male_late),
           p_sex_antag_late = pp_female_late * (1 - pp_male_late) + 
             (1 - pp_female_late) * pp_male_late,
           p_age_concord_females = pp_female_early * pp_female_late + 
             (1 - pp_female_early) * (1 - pp_female_late),
           p_age_antag_females = pp_female_early * (1 - pp_female_late) + 
             (1 - pp_female_early) * pp_female_late,
           p_age_concord_males = pp_male_early * pp_male_late + (1 - pp_male_early) * (1 - pp_male_late),
           p_age_antag_males = pp_male_early * (1 - pp_male_late) + (1 - pp_male_early) * pp_male_late) %>%
    
    # Find the ratios of some of these two probabilities (i.e. the "evidence ratios")
    mutate(inter_sex_early = p_sex_concord_early / p_sex_antag_early,
           inter_sex_late = p_sex_concord_late / p_sex_antag_late,
           inter_age_females = p_age_concord_females / p_age_antag_females,
           inter_age_males = p_age_concord_males / p_age_antag_males) 
}
test <- Wong_Holman %>% get_antagonism_ratios()
head(test)
rm(test)

all_genes_Wong <- tibble(merge(all.genes, Chrs, by = "FlyBaseID", all = T))
all_genes_Wong <- all_genes_Wong %>% 
  mutate(Sig.Af = FlyBaseID %in% A.f.geno[A.f.geno$Sig,]$FlyBaseID,
         Sig.Am = FlyBaseID %in% A.m.geno[A.m.geno$Sig,]$FlyBaseID,
         Sig = Sig.Af | Sig.Am,
         IsWong = FlyBaseID %in% Wong_Holman[Wong_Holman$SA,]$FlyBaseID)

all_genes_Wong_Chr2 <- all_genes_Wong[all_genes_Wong$Chr == "2",]

test <- fisher.test(x = all_genes_Wong_Chr2$Sig, y = all_genes_Wong_Chr2$IsWong)

mos_plot_Wong
mos_plot_Wong_Chr2 <- ggbarstats(
  all_genes_Wong_Chr2[!is.na(all_genes_Wong_Chr2$Sig) & 
                       !is.na(all_genes_Wong_Chr2$IsWong),], IsWong, Sig,
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

all_genes_Chloe <- tibble(merge(all.genes, Chrs, by = "FlyBaseID", all = T))
all_genes_Chloe <- all_genes_Chloe %>% 
  mutate(Sig.Af = FlyBaseID %in% A.f.geno[A.f.geno$Sig,]$FlyBaseID,
         Sig.Am = FlyBaseID %in% A.m.geno[A.m.geno$Sig,]$FlyBaseID,
         Sig = Sig.Af | Sig.Am,
         IsChloe = FlyBaseID %in% SSAV_chloe$FlyBaseID)


all_genes_Chloe_Chr2 <- all_genes_Chloe[all_genes_Chloe$Chr == "2",]

test <- fisher.test(x = all_genes_Chloe_Chr2$Sig, y = all_genes_Chloe_Chr2$IsChloe)

mos_plot_Chloe # the genomics candidates only consists of Chr2 genes. So filtering out the data below...

mos_plot_Chloe_Chr2 <- ggbarstats(
  all_genes_Chloe_Chr2[!is.na(all_genes_Chloe_Chr2$Sig) & 
                         !is.na(all_genes_Chloe_Chr2$IsChloe),], IsChloe, Sig,
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


pdf(file = "~/Desktop/UofT/SSAV_RNA/Plots/Enrichment_Tests/Ruz_Chr2.pdf",  # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 7) # The height of the plot in inches
mos_plot_Ruz_Chr2 # mos_plot_InnoMorr_Chr2, mos_plot_Ruz, mos_plot_Ruz_Chr2, mos_plot_Chloe_Chr2, mos_plot_Wong, mos_plot_Wong_Chr2
dev.off()


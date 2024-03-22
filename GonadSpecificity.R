###################################
#
#                             Grieshop et al. 2023
#                             Author: Michelle Liu
#             DsRed experimental evolution - transcriptomics analysis
#         Comparing gonad specificity between candidates vs non-candidates
# 
# 
###################################

# packages 
#########
library(tidyr)
library(plyr)
library(dplyr)
library(broom)
library(ggplot2)
library(ggblend)
library(ggpubr)
library(RColorBrewer)
########

## prepare dataset
########
# load results if not loaded in env.
SSAV.geno <- read.delim("Results/All.geno_candidates.tsv")

# load Gonad specificity data
Gonad <- read.csv(file="~/Desktop/UofT/SSAV_RNA/Data/dtGonadSpecifityVals.csv", header=TRUE)
colnames(Gonad)[1] <- "FlyBaseID"

# combine with Gonad specificity data
SSAV.geno <- merge(SSAV.geno, Gonad, by = "FlyBaseID", all = TRUE)
# include SBGE categories (using Mishra et al. dataset. Look at External_data.R)
SSAV.geno <- merge(SSAV.geno, ASE, by = "FlyBaseID", all = TRUE)
SSAV.geno <- SSAV.geno[!is.na(SSAV.geno$Sig) & 
                         !is.na(SSAV.geno$geoMeanGSI) & 
                         !is.na(SSAV.geno$SBGE_comp),]
########


## plotting function 
pointSEplot <- function(boot_dat, perm_dat, x_col, SBGE_cat = NA){
  # set y-axis value above each error bar
  if(!is.na(SBGE_cat)){
    y_count_FALSE <- boot_dat[!boot_dat$Sig,] %>% 
      group_by(SBGE_comp) %>%
      summarise(max = q95 + 0.05)
    y_count_TRUE <- boot_dat[boot_dat$Sig,] %>%
      group_by(SBGE_comp) %>%
      summarise(max = q95 + 0.05)
  }
  else{
    y_count_FALSE <- boot_dat[!boot_dat$Sig,] %>% 
      summarise(max = q95 + 0.05)
    y_count_TRUE <- boot_dat[boot_dat$Sig,] %>% 
      summarise(max = q95 + 0.05)
    perm_dat <- as.data.frame(t(perm_dat))
    colnames(perm_dat) <- c("obs_diff", "n_TRUE", "n_FALSE", "pval", "Sig")
    str(perm_dat)
  }
  # join y-axis coordinate with text dataframes
  perm_dat$y_count_FALSE <- y_count_FALSE$max
  perm_dat$y_count_TRUE <- y_count_TRUE$max
  perm_dat$y_star <- perm_dat$y_count_FALSE + 0.1
  
  if(!is.na(SBGE_cat)){
    pointSE_plot <- ggplot(boot_dat) + 
      geom_point(aes_string(SBGE_cat, "q50", color = "Sig"), 
                 position = position_jitterdodge(jitter.width = 0.001), size = 6) +
      geom_errorbar(aes_string(SBGE_cat, "q50", color = "Sig"),
                    ymin = boot_dat$q05, ymax = boot_dat$q95, 
                    position = position_jitterdodge(jitter.width = 0.01), width = 0.5) +
      labs(x = "SBGE (ASE)", # "omegaA_MK" = expression(italic("\u03c9A")[MK]); "alpha_MK" = expression(italic("\u03b1")[MK])
           y = print(x_col)) +
      geom_text(data = perm_dat %>% mutate(sig1 = if_else(Sig==1, "*", "")),
                aes(x =SBGE_comp, y = y_star, label = sig1), color = "black", size = 10) +
      geom_text(data = perm_dat,aes(x =SBGE_comp, y = y_count_FALSE, label = n_FALSE), color = "black", size = 6,
                hjust = 1) +
      geom_text(data = perm_dat,aes(x =SBGE_comp, y = y_count_TRUE, label = n_TRUE), color = "black", 
                size = 5.5, hjust = -1) +
      scale_x_discrete(labels = c("Highly FB","Female-Biased", "Unbiased", "Male-Biased", "Highly MB")) 
  }
  else{
    pointSE_plot <- ggplot(boot_dat) +
      geom_point(aes(x = Sig, y = q50, color = Sig),
                 position = position_jitterdodge(jitter.width = 0.001), size = 6) +
      geom_errorbar(aes(x = Sig, y = q50, color = Sig),
                    ymin = boot_dat$q05, ymax = boot_dat$q95, 
                    position = position_jitterdodge(jitter.width = 0.01), width = 0.5) +
      labs(y = print(x_col)) +
      geom_text(data = perm_dat %>% mutate(sig1 = if_else(Sig==1, "*", "")),
                aes(x = 1.5, y = y_star, label = sig1), color = "black", size = 10) +
      geom_text(data = perm_dat,aes(x = 1.04, y = y_count_FALSE, label = n_FALSE), 
                color = "black", size = 6, hjust = 1) +
      geom_text(data = perm_dat,aes(x = 1.9, y = y_count_TRUE, label = n_TRUE), color = "black", 
                size = 5.5, hjust = -1) 
  }
  
  pointSE_plot <- pointSE_plot +
    scale_colour_manual(values = c("orange3", "forestgreen"), # "red3", "steelblue3", "#888888" # "purple3", "chartreuse3", "orange2"
                        labels = c("background", "candidates")) + # "Chr-2", "Chr-3", "X-Chr"
    guides(color = guide_legend(override.aes = list(shape = c(16, 16),
                                                    size = c(4, 4),
                                                    alpha = 1))) +
    scale_y_continuous(limits = c(0,1)) +
    theme_classic() +
    theme(plot.title.position = c("panel"),
          legend.title = element_blank(),
          legend.position = c("bottom"),
          legend.box.background = element_rect(),
          legend.text = element_text(size = 20, color = "black"),
          axis.text.x = element_text(size=20, margin = margin(5,0,0,0), color = "black"),
          axis.text.y = element_text(size=20, margin = margin(0,5,0,0), color = "black"),
          axis.title.x = element_text(size=40, margin = margin(10,0,0,0), color = "black"),
          axis.title.y = element_text(size=40, margin = margin(0,10,0,0), color = "black"),
          plot.title = element_text(size=40, margin = margin(0,0,0,0), color = "black"),
          plot.margin = margin(6,6,6,6)
    )
  
  return(pointSE_plot)
}


## Candidates vs Non-candidates
########
boot_testes_All
boot_ovaries_All <- TwoBoot_SBGE(SSAV.geno[!is.na(SSAV.geno$ovariesSpecificity),], 
                                x_col = "ovariesSpecificity", 
                                groupBy = "Sig",
                                SBGE_cat = "SBGE_comp")

perm_testes_All <- TwoPerm_SBGE(SSAV.geno[!is.na(SSAV.geno$testesSpecificity),], 
                                x_col = "testesSpecificity", 
                                groupBy = "Sig",
                                SBGE_cat = "SBGE_comp")
perm_ovaries_All <- TwoPerm_SBGE(SSAV.geno[!is.na(SSAV.geno$ovariesSpecificity),], 
                                x_col = "ovariesSpecificity", 
                                groupBy = "Sig",
                                SBGE_cat = "SBGE_comp")

# plot ovaries specificity
ovaries_All_SBGE <- pointSEplot(boot_dat = boot_ovaries_All, 
                                perm_dat = perm_ovaries_All,
                                x_col = "ovariesSpecificity", 
                                SBGE_cat = "SBGE_comp") + coord_cartesian(ylim = c(0,1)) 

# plot testes specificity
testes_All_SBGE <- pointSEplot(boot_dat = boot_testes_All, 
                               perm_dat = perm_testes_All,
                              x_col = "testesSpecificity", 
                              SBGE_cat = "SBGE_comp") + coord_cartesian(ylim = c(0,1)) 

########


## Proportion of most differentiated genes
# SSAV.geno_ASE <-  SSAV.geno_ASE %>% mutate(Most_A.m = ifelse((A.m.exp_geno > quantile(A.m.exp_geno, 0.875, na.rm = T) |
#                                               A.m.exp_geno < quantile(A.m.exp_geno, 0.125, na.rm = T)), 
#                                            TRUE, FALSE),
#                          Most_A.f = ifelse((A.f.exp_geno > quantile(A.f.exp_geno, 0.875, na.rm = T) |
#                                               A.f.exp_geno < quantile(A.f.exp_geno, 0.125, na.rm = T)), 
#                                            TRUE, FALSE))

## most upregulated
SSAV.geno_ASE <-  SSAV.geno_ASE %>% mutate(Most_A.m = ifelse(A.m.exp_geno > quantile(A.m.exp_geno, 0.75, na.rm = T), 
                                                             TRUE, FALSE),
                                           Most_A.f = ifelse(A.f.exp_geno > quantile(A.f.exp_geno, 0.75, na.rm = T), 
                                                             TRUE, FALSE))
# only Highly MB genes
High_MB.geno_ASE <- SSAV.geno_ASE[SSAV.geno_ASE$SBGE_comp == "e.more.mbg",]
High_MB.geno_ASE <-  High_MB.geno_ASE %>% mutate(t5_A.m = ifelse(A.m.exp_geno > quantile(A.m.exp_geno, 0.95, na.rm = T), 
                                                                 TRUE, FALSE),
                                                t10_A.m = ifelse(A.m.exp_geno > quantile(A.m.exp_geno, 0.90, na.rm = T), 
                                                                 TRUE, FALSE),
                                                t25_A.m = ifelse(A.m.exp_geno > quantile(A.m.exp_geno, 0.75, na.rm = T), 
                                                                 TRUE, FALSE))

# Test for Highly Mb genes
######
t5_HiMB_boot <- TwoBoot(High_MB.geno_ASE[!is.na(High_MB.geno_ASE$testesSpecificity) &
                                             !is.na(High_MB.geno_ASE$t5_A.m) &
                                             !is.na(High_MB.geno_ASE$exp_SBGE_ase),], 
                             x_col = "testesSpecificity", 
                             groupBy = "t5_A.m")
t5_HiMB_perm <- TwoPerm(High_MB.geno_ASE[!is.na(High_MB.geno_ASE$testesSpecificity) &
                                           !is.na(High_MB.geno_ASE$t5_A.m) &
                                           !is.na(High_MB.geno_ASE$exp_SBGE_ase),], 
                        x_col = "testesSpecificity", 
                        groupBy = "t5_A.m")

t10_HiMB_boot <- TwoBoot(High_MB.geno_ASE[!is.na(High_MB.geno_ASE$testesSpecificity) &
                                           !is.na(High_MB.geno_ASE$t10_A.m) &
                                           !is.na(High_MB.geno_ASE$exp_SBGE_ase),], 
                        x_col = "testesSpecificity", 
                        groupBy = "t10_A.m")
t10_HiMB_perm <- TwoPerm(High_MB.geno_ASE[!is.na(High_MB.geno_ASE$testesSpecificity) &
                                           !is.na(High_MB.geno_ASE$t10_A.m) &
                                           !is.na(High_MB.geno_ASE$exp_SBGE_ase),], 
                        x_col = "testesSpecificity", 
                        groupBy = "t10_A.m")
t25_HiMB_boot <- TwoBoot(High_MB.geno_ASE[!is.na(High_MB.geno_ASE$testesSpecificity) &
                                            !is.na(High_MB.geno_ASE$t25_A.m) &
                                            !is.na(High_MB.geno_ASE$exp_SBGE_ase),], 
                         x_col = "testesSpecificity", 
                         groupBy = "t25_A.m")
t25_HiMB_perm <- TwoPerm(High_MB.geno_ASE[!is.na(High_MB.geno_ASE$testesSpecificity) &
                                            !is.na(High_MB.geno_ASE$t25_A.m) &
                                            !is.na(High_MB.geno_ASE$exp_SBGE_ase),], 
                         x_col = "testesSpecificity", 
                         groupBy = "t25_A.m")
  
######


## ovaries
#######
boot_ovaries_Af_most <- TwoBoot_SBGE(SSAV.geno_ASE[!is.na(SSAV.geno_ASE$ovariesSpecificity) &
                                                    !is.na(SSAV.geno_ASE$Most_A.f) &
                                                     !is.na(SSAV.geno_ASE$exp_SBGE_ase),], 
                                    x_col = "ovariesSpecificity", 
                                    groupBy = "Most_A.f",
                                    SBGE_cat = "SBGE_comp")
perm_ovaries_Af_most <- TwoPerm_SBGE(SSAV.geno_ASE[!is.na(SSAV.geno_ASE$ovariesSpecificity) &
                                                    !is.na(SSAV.geno_ASE$Most_A.f) &
                                                     !is.na(SSAV.geno_ASE$exp_SBGE_ase),], 
                                    x_col = "ovariesSpecificity", 
                                    groupBy = "Most_A.f",
                                    SBGE_cat = "SBGE_comp")


boot_ovaries_Am_most <- TwoBoot_SBGE(SSAV.geno_ASE[!is.na(SSAV.geno_ASE$ovariesSpecificity) &
                                                    !is.na(SSAV.geno_ASE$Most_A.m) &
                                                     !is.na(SSAV.geno_ASE$exp_SBGE_ase),], 
                                    x_col = "ovariesSpecificity", 
                                    groupBy = "Most_A.m",
                                    SBGE_cat = "SBGE_comp")
perm_ovaries_Am_most <- TwoPerm_SBGE(SSAV.geno_ASE[!is.na(SSAV.geno_ASE$ovariesSpecificity) &
                                                    !is.na(SSAV.geno_ASE$Most_A.m) &
                                                     !is.na(SSAV.geno_ASE$exp_SBGE_ase),], 
                                    x_col = "ovariesSpecificity", 
                                    groupBy = "Most_A.m",
                                    SBGE_cat = "SBGE_comp")

# plot ovaries specificity for DE vs not-DE genes in SSAV females
ovaries_Af_most <- pointSEplot(boot_dat = boot_ovaries_Af_most, 
                              perm_dat = perm_ovaries_Af_most, 
                              x_col = "ovariesSpecificity", 
                              SBGE_cat = "SBGE_comp") + 
  coord_cartesian(ylim = c(0,1)) +
  scale_colour_manual(values = c("sienna", "olivedrab"), 
                      labels = c("not DE", "DE"))# "Chr-2", "Chr-3", "X-Chr"

# plot ovaries specificity for DE vs not-DE genes in SSAV males
ovaries_Am_most <- pointSEplot(boot_dat = boot_ovaries_Am_most, 
                              perm_dat = perm_ovaries_Am_most, 
                              x_col = "ovariesSpecificity",
                              SBGE_cat = "SBGE_comp") + 
  coord_cartesian(ylim = c(0,1)) +
  scale_colour_manual(values = c("sienna", "olivedrab"), 
                      labels = c("not DE", "DE"))# "Chr-2", "Chr-3", "X-Chr"

#######



## testes
#######
boot_testes_Af_most <- TwoBoot_SBGE(SSAV.geno_ASE[!is.na(SSAV.geno_ASE$testesSpecificity) &
                                                       !is.na(SSAV.geno_ASE$Most_A.f) &
                                                    !is.na(SSAV.geno_ASE$exp_SBGE_ase),], 
                                       x_col = "testesSpecificity", 
                                       groupBy = "Most_A.f",
                                       SBGE_cat = "SBGE_comp")
perm_testes_Af_most <- TwoPerm_SBGE(SSAV.geno_ASE[!is.na(SSAV.geno_ASE$testesSpecificity) &
                                                   !is.na(SSAV.geno_ASE$Most_A.f) &
                                                    !is.na(SSAV.geno_ASE$exp_SBGE_ase),], 
                                    x_col = "testesSpecificity", 
                                    groupBy = "Most_A.f",
                                    SBGE_cat = "SBGE_comp")


boot_testes_Am_most <- TwoBoot_SBGE(SSAV.geno_ASE[!is.na(SSAV.geno_ASE$testesSpecificity) &
                                           !is.na(SSAV.geno_ASE$Most_A.m) &
                                             !is.na(SSAV.geno_ASE$exp_SBGE_ase),], 
                                    x_col = "testesSpecificity", 
                                    groupBy = "Most_A.m",
                                    SBGE_cat = "SBGE_comp")
perm_testes_Am_most <- TwoPerm_SBGE(SSAV.geno_ASE[!is.na(SSAV.geno_ASE$testesSpecificity) &
                                                    !is.na(SSAV.geno_ASE$Most_A.m) &
                                                    !is.na(SSAV.geno_ASE$exp_SBGE_ase),], 
                                    x_col = "testesSpecificity", 
                                    groupBy = "Most_A.m",
                                    SBGE_cat = "SBGE_comp")

# plot testes specificity for DE vs not-DE genes in SSAV females
testes_Af_most <- pointSEplot(boot_dat = boot_testes_Af_most, 
                              perm_dat = perm_testes_Af_most, 
                              x_col = "testesSpecificity", 
                              SBGE_cat = "SBGE_comp") + 
  coord_cartesian(ylim = c(0,1)) +
  scale_colour_manual(values = c("sienna", "olivedrab"), 
                      labels = c("not DE", "DE"))# "Chr-2", "Chr-3", "X-Chr"

# plot testes specificity for DE vs not-DE genes in SSAV males
testes_Am_most <- pointSEplot(boot_dat = boot_testes_Am_most, 
                              perm_dat = perm_testes_Am_most, 
                              x_col = "testesSpecificity",
                              SBGE_cat = "SBGE_comp") + 
  coord_cartesian(ylim = c(0,1)) +
  scale_colour_manual(values = c("sienna", "olivedrab"), 
                      labels = c("not DE", "DE"))# "Chr-2", "Chr-3", "X-Chr"

#######





pdf(file = "~/Desktop/UofT/SSAV_RNA/Plots/Gonads_DEvNDE.pdf",   # The directory you want to save the file in
    width = 24, # 12, 24 The width of the plot in inches
    height = 20) # 10, 20 The height of the plot in inches
# ggarrange(testes_All_SBGE, NA, ovaries_All_SBGE,
#           labels = c("A)", NA, "B)"),
#           widths = c(1, 0.05, 1),
#           ncol = 3,
#           font.label = list(size = 30))
ggarrange(testes_Af_most, NA, ovaries_Af_most,
          NA,NA,NA,
          testes_Am_most, NA, ovaries_Am_most,
          labels = c("A)", NA, "B)",
                     NA, NA, NA,
                     "C", NA, "D)"),
          widths = c(1, 0.05, 1),
          heights = c(1, 0.05, 1),
          ncol = 3, nrow = 3,
          font.label = list(size = 30))
dev.off()





### -------------- general tissue specificity measures (tao) -----------

## prepare dataset using Singh & Agrawal 2023 data
##########
# load results if not loaded in env.
A.f.geno <- read.delim("Results/A.f.geno_candidates.tsv")
A.m.geno <- read.delim("Results/A.m.geno_candidates.tsv")
# Combine results into one data frame
# Genes present in both SSAV males and SSAV females data
SSAV.geno <- merge(A.m.geno, A.f.geno, by = "FlyBaseID", all = TRUE)
colnames(SSAV.geno) <- c("FlyBaseID", "A.m.exp_geno", "A.m.se_geno", "A.m.padj", "A.m.TopSig", "A.m.Sig",
                         "A.f.exp_geno", "A.f.se_geno", "A.f.padj", "A.f.Sig")
# column denotes genes that are candidates in males or females
SSAV.geno <- SSAV.geno %>% mutate(Sig = ifelse(!is.na(A.m.Sig) & A.m.Sig, TRUE, 
                                               ifelse(!is.na(A.f.Sig) & A.f.Sig, TRUE, FALSE))) 


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

SSAV.geno <- merge(SSAV.geno, SDIU, by = "FlyBaseID", all = TRUE)
##########

# this uses the SBGE estimations from Osada et al. 
perm_tao <- TwoPerm_SBGE(SSAV.geno[!is.na(SSAV.geno$Sig) & 
                                !is.na(SSAV.geno$tao.FML.FA2) & 
                                  !is.na(SSAV.geno$SBGEcat.body.Osada),], 
                         x_col = "tao.FML.FA2", groupBy = "Sig", SBGE_cat = "SBGEcat.body.Osada")
boot_tao <- TwoBoot_SBGE(SSAV.geno[!is.na(SSAV.geno$Sig) & 
                                !is.na(SSAV.geno$tao.FML.FA2) & 
                                  !is.na(SSAV.geno$SBGEcat.body.Osada),], 
                    x_col = "tao.FML.FA2", groupBy = "Sig", SBGE_cat = "SBGEcat.body.Osada")
tao_all <- pointSEplot(boot_dat = boot_tao, perm_dat = perm_tao,
                       x_col = "tao.FML.FA2", SBGE_cat = "SBGEcat.body.Osada") +
  coord_cartesian(ylim = c(0, 1.2)) + 
  scale_x_discrete(labels = c("Extreme FB","Strong FB","Female-Biased", 
                              "Unbiased", "Male-Biased", "Strong MB", "Extreme MB"))




###################################
#
#                             Grieshop et al. 2025
#                             Author: Michelle Liu
#             DsRed experimental evolution - transcriptomics analysis
#         Comparing gonad specificity between DE candidates vs non-candidates
#                                   Table S6
# 
# 
###################################

rm(list=ls())
setwd("~/Desktop/UofT/SSAV_RNA/")

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

# load Gonad specificity data
Gonad <- read.csv(file="~/Desktop/UofT/SSAV_RNA/Data/dtGonadSpecifityVals.csv", header=TRUE)
colnames(Gonad)[1] <- "FlyBaseID"

source("Mishra_et.al_SBGE.R")

## prepare dataset
########
# load results if not loaded in env.
SSAV.geno <- read.delim("Results/All.geno_DE.candidates.tsv")


# combine with Gonad specificity data
SSAV.geno <- merge(SSAV.geno, Gonad, by = "FlyBaseID", all = TRUE)
# include SBGE categories (using Mishra et al. dataset. Look at External_data.R)
SSAV.geno <- merge(SSAV.geno, Mishra, by = "FlyBaseID", all = TRUE)
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
      labs(x = "SBGE (Mishra)", # "omegaA_MK" = expression(italic("\u03c9A")[MK]); "alpha_MK" = expression(italic("\u03b1")[MK])
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
# bootstrap functions
source("boot_permute.R")


## most upregulated
SSAV.geno <-  SSAV.geno %>% mutate(Most_A.m = ifelse(A.m.exp_geno > quantile(A.m.exp_geno, 0.75, na.rm = T), 
                                              TRUE, FALSE),
                                   Most_A.f = ifelse(A.f.exp_geno > quantile(A.f.exp_geno, 0.75, na.rm = T), 
                                              TRUE, FALSE))
# only Highly MB genes
High_MB <- SSAV.geno[SSAV.geno$SBGE_comp == "e.more.mbg",]
High_MB <-  High_MB %>% mutate(t5_A.m = ifelse(A.m.exp_geno > quantile(A.m.exp_geno, 0.95, na.rm = T), 
                                                                 TRUE, FALSE),
                                t10_A.m = ifelse(A.m.exp_geno > quantile(A.m.exp_geno, 0.90, na.rm = T), 
                                                                 TRUE, FALSE),
                                 t25_A.m = ifelse(A.m.exp_geno > quantile(A.m.exp_geno, 0.75, na.rm = T), 
                                                                 TRUE, FALSE))

# Test for Highly Mb genes (log2FC > 5)
# Are the most upregulated between Red and NonRed genes testes-specific genes?
# If yes, might indicate that the change is due to differences in gonad size in Red males
# If no, then differences could be due to genuine upregulation in Red of genes of this SBGE category
######
# Only include assayed genes with complete data on sex-bias and testes specificity

# top 5% upregulated vs. others
# run bootstrapping function. 
t5_HiMB_boot <- TwoBoot(High_MB[!is.na(High_MB$testesSpecificity) &
                                 !is.na(High_MB$t5_A.m) &
                                  !is.na(High_MB$exp_SBGE_ase),], 
                             x_col = "testesSpecificity", 
                             groupBy = "t5_A.m")
# permute means
t5_HiMB_perm <- TwoPerm(High_MB[!is.na(High_MB$testesSpecificity) &
                                 !is.na(High_MB$t5_A.m) &
                                  !is.na(High_MB$exp_SBGE_ase),], 
                        x_col = "testesSpecificity", 
                        groupBy = "t5_A.m", n_perm = 1000)
# plot
t5_plot <- pointSEplot(boot_dat = t5_HiMB_boot %>% 
                         dplyr::rename("Sig" = "t5_A.m"), 
                        perm_dat = t5_HiMB_perm,
                           x_col = "testesSpecificity") + coord_cartesian(ylim = c(0,1)) 


# top 10% upregulated vs. others
t10_HiMB_boot <- TwoBoot(High_MB[!is.na(High_MB$testesSpecificity) &
                                  !is.na(High_MB$t10_A.m) &
                                   !is.na(High_MB$exp_SBGE_ase),], 
                        x_col = "testesSpecificity", 
                        groupBy = "t10_A.m")
t10_HiMB_perm <- TwoPerm(High_MB[!is.na(High_MB$testesSpecificity) &
                                  !is.na(High_MB$t10_A.m) &
                                   !is.na(High_MB$exp_SBGE_ase),], 
                        x_col = "testesSpecificity", 
                        groupBy = "t10_A.m", n_perm = 10000)
t10_plot <- pointSEplot(boot_dat = t10_HiMB_boot %>% 
                          dplyr::rename("Sig" = "t10_A.m"), 
                       perm_dat = t10_HiMB_perm,
                       x_col = "testesSpecificity") + coord_cartesian(ylim = c(0,1)) 


# top 25% upregulated vs. others
t25_HiMB_boot <- TwoBoot(High_MB[!is.na(High_MB$testesSpecificity) &
                                  !is.na(High_MB$t25_A.m) &
                                   !is.na(High_MB$exp_SBGE_ase),], 
                         x_col = "testesSpecificity", 
                         groupBy = "t25_A.m")
t25_HiMB_perm <- TwoPerm(High_MB[!is.na(High_MB$testesSpecificity) &
                                  !is.na(High_MB$t25_A.m) &
                                    !is.na(High_MB$exp_SBGE_ase),], 
                         x_col = "testesSpecificity", 
                         groupBy = "t25_A.m")
t25_plot <- pointSEplot(boot_dat = t25_HiMB_boot %>% 
                          dplyr::rename("Sig" = "t25_A.m"), 
                       perm_dat = t25_HiMB_perm,
                       x_col = "testesSpecificity") + coord_cartesian(ylim = c(0,1)) 
  
######



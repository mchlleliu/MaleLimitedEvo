###################################
#
#                             Grieshop et al. 2023
#             DsRed experimental evolution - transcriptomics analysis
#                  TajD and DoS - Signals of Balancing Selection?
# 
# 
###################################

# Do candidate sexually antagonistic genes show signals of balancing selection as 
# predicted by previous theoretical models of sexual antagonism? 
# (but see Connallon and Clark paper about balancing selection)


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

## plotting function
pointSEplot <- function(boot_dat, perm_dat, x_col, SBGE_cat = NA){
  # set y-axis value above each error bar
  if(!is.na(SBGE_cat)){
    y_count_FALSE <- boot_dat[!boot_dat$Sig,] %>% 
      group_by(.[[SBGE_cat]]) %>%
      summarise(max = q95 + 0.05) %>% 
      rename({{SBGE_cat}} := 1)
    y_count_TRUE <- boot_dat[boot_dat$Sig,] %>%
      group_by(.[[SBGE_cat]]) %>%
      summarise(max = q95 + 0.05) %>%
      rename({{SBGE_cat}} := 1)
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
                aes_string(x =SBGE_cat, y = "y_star", label = "sig1"), color = "black", size = 10) +
      geom_text(data = perm_dat,aes_string(x =SBGE_cat, y = "y_count_FALSE", label = "n_FALSE"), 
                color = "black", size = 6, hjust = 1) +
      geom_text(data = perm_dat,aes_string(x =SBGE_cat, y = "y_count_TRUE", label = "n_TRUE"), 
                color = "black", size = 5.5, hjust = -1) +
      # change the labels accordingly if needed.
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


#### DoS candidates vs non-candidates
# Calculations were obtained from Fraisse et al. 2019 and included in the Singh & Agrawal 2023 paper.
# DoS = (Dn/(Dn+Ds)) - (Pn/(Pn+Ps))
# DoS < 0 would suggest more intraspecific polymorphism (balancing selection, relaxed purifying selection)
# DoS > 0 would suggest more interspecific divergence (positive selection)
## Overall DoS > 0
## No significant difference in DoS between candidate vs non-candidate genes
## Candidates that are strong/extremely MB seem to have lower DoS estimates.
#########
# use permutation and bootstrap functions in boot_permute.R
perm_DoS <- TwoPerm(SSAV.geno[!is.na(SSAV.geno$Sig) & 
                                !is.na(SSAV.geno$DoS),], x_col = "DoS", groupBy = "Sig")
boot_DoS <- TwoBoot(SSAV.geno[!is.na(SSAV.geno$Sig) & 
                                !is.na(SSAV.geno$DoS),], x_col = "DoS", groupBy = "Sig")
DoS_all <- pointSEplot(boot_dat = boot_DoS, perm_dat = perm_DoS, x_col = "DoS")


# by SBGE category (according to Osada et al. as used in the Sing & Agrawal paper)
perm_DoS_SBGE <- TwoPerm_SBGE(SSAV.geno[!is.na(SSAV.geno$Sig) & 
                                !is.na(SSAV.geno$DoS) & !is.na(SSAV.geno$SBGEcat.body.Osada),], 
                         x_col = "DoS", 
                         groupBy = "Sig", 
                         SBGE_cat = "SBGEcat.body.Osada")
boot_DoS_SBGE <- TwoBoot_SBGE(SSAV.geno[!is.na(SSAV.geno$Sig) & 
                                !is.na(SSAV.geno$DoS) & !is.na(SSAV.geno$SBGEcat.body.Osada),],
                         x_col = "DoS", groupBy = "Sig", SBGE_cat = "SBGEcat.body.Osada")
DoS_all_SBGE <- pointSEplot(boot_dat = boot_DoS_SBGE, perm_dat = perm_DoS_SBGE, 
                            x_col = "DoS", SBGE_cat = "SBGEcat.body.Osada") + 
  scale_x_discrete(labels = c("Extreme FB","Strong FB","Female-Biased", 
                              "Unbiased", "Male-Biased", "Strong MB", "Extreme MB"))
#########



#### TajD
# D = d/SE(d), where d = pi - watterson's theta
# D < 0 implies directional selection or selective sweeps
# D > 0 implies balancing selection or relaxed purifying selection 
# unresolved problems with numbering here...
## generally negative, but not strongly so.
#########
# use permutation and bootstrap functions in boot_permute.R
perm_TajD <- TwoPerm(SSAV.geno[!is.na(SSAV.geno$Sig) & 
                                !is.na(SSAV.geno$tajD.N),], x_col = "tajD.N", groupBy = "Sig")
boot_TajD <- TwoBoot(SSAV.geno[!is.na(SSAV.geno$Sig) & 
                                !is.na(SSAV.geno$tajD.N),], x_col = "tajD.N", groupBy = "Sig")
# numbers are gone... not rlly sure why.
TajD_all <- pointSEplot(boot_dat = boot_TajD, perm_dat = perm_TajD, x_col = "tajD.N") + coord_cartesian(ylim = c(-0.5, 0.5))


# by SBGE category (according to Osada et al. as used in the Sing & Agrawal paper)
perm_TajD_SBGE <- TwoPerm_SBGE(SSAV.geno[!is.na(SSAV.geno$Sig) & 
                                          !is.na(SSAV.geno$tajD.N) & 
                                           !is.na(SSAV.geno$SBGEcat.body.Osada),], 
                              x_col = "tajD.N", 
                              groupBy = "Sig", 
                              SBGE_cat = "SBGEcat.body.Osada")
boot_TajD_SBGE <- TwoBoot_SBGE(SSAV.geno[!is.na(SSAV.geno$Sig) & 
                                          !is.na(SSAV.geno$tajD.N) & 
                                           !is.na(SSAV.geno$SBGEcat.body.Osada),],
                              x_col = "tajD.N", groupBy = "Sig", SBGE_cat = "SBGEcat.body.Osada")

# numbers are gone... not rlly sure why.
TajD_all_SBGE <- pointSEplot(boot_dat = boot_TajD_SBGE, perm_dat = perm_TajD_SBGE, 
                            x_col = "tajD.N", SBGE_cat = "SBGEcat.body.Osada") + 
  scale_x_discrete(labels = c("Extreme FB","Strong FB","Female-Biased", 
                              "Unbiased", "Male-Biased", "Strong MB", "Extreme MB")) + coord_cartesian(ylim = c(-0.5, 0.5))
#########






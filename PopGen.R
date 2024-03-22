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
library(tidyverse)
library(tidyr)
library(plyr)
library(dplyr)
library(broom)
library(ggplot2)
library(ggblend)
library(ggpubr)
library(RColorBrewer)
########

# load Singh & Agrawal 2023 dataset
######
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
######

# Get Mishra et al.'s data 
##########

ASE <- read.delim(file="~/Desktop/UofT/SSAV_RNA/Data/DifferentialGeneExpression.whole.bodies.tsv", sep="\t", header=TRUE)
# Grab the desired variables from there
ASE <- data.frame(cbind(ASE$log2FoldChange,
                        ASE$lfcSE,
                        ASE$FlyBaseID))
colnames(ASE) <- c("exp_SBGE_ase", "se_SBGE_ase", "FlyBaseID")
# fix formatting
ASE$exp_SBGE_ase <- as.numeric(ASE$exp_SBGE_ase)
ASE$se_SBGE_ase <- as.numeric(ASE$se_SBGE_ase)
str(ASE)


# Define three levels of SBGE categorization 
x1 = 1 # first cut-off (FBG < -1, MBG > 1 , -1 < UBG < 1)
x2 = 5 # second cut-off (extreme FBG < -5, extreme MBG > 5)
# y0 = 1 # tolerance of middle bins ### I DONT GET THIS 
# # xmid1 = (x1 + x2)/2

# Simple (3 levels)
# one level of female-biased gene expression
fbg.keep <- ASE$exp_SBGE_ase < -x1 & (ASE$exp_SBGE_ase + ASE$se_SBGE_ase) < 0
fbg <- ASE[fbg.keep,]
fbg$SBGE_simp <- rep(c("a.fbg"), dim(fbg)[1])
# one unbiased category
ubg.keep <- ASE$exp_SBGE_ase < x1 & ASE$exp_SBGE_ase > -x1 & (ASE$exp_SBGE_ase - ASE$se_SBGE_ase) > -(x1+y0) & (ASE$exp_SBGE_ase + ASE$se_SBGE_ase) < (x1+y0)
ubg <- ASE[ubg.keep,]
ubg$SBGE_simp <- rep(c("b.ubg"), dim(ubg)[1])
# two levels of male-biased gene expression
mbg.keep <- ASE$exp_SBGE_ase > x1 & (ASE$exp_SBGE_ase - ASE$se_SBGE_ase) > 0
mbg <- ASE[mbg.keep,]
mbg$SBGE_simp <- rep(c("c.mbg"), dim(mbg)[1])
# 1 gene is tossed out b/c the uncertainty in its estimate breaches a cutoff boundary 
ASE <- rbind(fbg, mbg, ubg)
str(ASE)

# Complex (5 levels)
more.fbg.keep <- ASE$exp_SBGE_ase < -x2 & (ASE$exp_SBGE_ase + ASE$se_SBGE_ase) < -(x2-y0) # extreme FBG < -5
more.fbg <- ASE[more.fbg.keep,]
more.fbg$SBGE_comp <- rep(c("a.more.fbg"), dim(more.fbg)[1])
#
fbg.keep <- ASE$exp_SBGE_ase < -x1 & ASE$exp_SBGE_ase > -x2 & (ASE$exp_SBGE_ase + ASE$se_SBGE_ase) < 0 & (ASE$exp_SBGE_ase - ASE$se_SBGE_ase) > -(x2+y0)
fbg <- ASE[fbg.keep,]
fbg$SBGE_comp <- rep(c("b.fbg"), dim(fbg)[1])
# one unbiased category
ubg.keep <- ASE$exp_SBGE_ase < x1 & ASE$exp_SBGE_ase > -x1 & (ASE$exp_SBGE_ase - ASE$se_SBGE_ase) > -(x1+y0) & (ASE$exp_SBGE_ase + ASE$se_SBGE_ase) < (x1+y0)
ubg <- ASE[ubg.keep,]
ubg$SBGE_comp <- rep(c("c.ubg"), dim(ubg)[1])
# two levels of male-biased gene expression
mbg.keep <- ASE$exp_SBGE_ase > x1 & ASE$exp_SBGE_ase < x2 & (ASE$exp_SBGE_ase - ASE$se_SBGE_ase) > 0 & (ASE$exp_SBGE_ase + ASE$se_SBGE_ase) < (x2+y0)
mbg <- ASE[mbg.keep,]
mbg$SBGE_comp <- rep(c("d.mbg"), dim(mbg)[1])
#
more.mbg.keep <- ASE$exp_SBGE_ase > x2 & (ASE$exp_SBGE_ase - ASE$se_SBGE_ase) > (x2-y0)
more.mbg <- ASE[more.mbg.keep,]
more.mbg$SBGE_comp <- rep(c("e.more.mbg"), dim(more.mbg)[1])
# 3 genes tossed out  b/c the uncertainty in their estimate breaches a cutoff boundary
ASE <- rbind(more.fbg, fbg, ubg, mbg, more.mbg)
str(ASE)

##########

## prepare dataset
##########
# load results if not loaded in env.
A.f.geno <- read.delim("Results/A.f.geno_candidates.tsv")
A.m.geno <- read.delim("Results/A.m.geno_candidates.tsv")
SSAV.geno <- read.delim("Results/All.geno_candidates.tsv")


jseq.All.geno.tmp <- read.delim("Results/jseq.All.geno.txt")
jseq.All.geno.tmp <- merge(jseq.All.geno.tmp, SDIU, by = "FlyBaseID")
jseq.All.geno.tmp <- jseq.All.geno.tmp[!is.na(jseq.All.geno.tmp$Sig) & !is.na(jseq.All.geno.tmp$tajD.S),]

sum(jseq.All.geno.tmp$Sig)
all.candidates <- SSAV.geno %>% 
  dplyr::mutate(Sig = ifelse(FlyBaseID %in% jseq.All.geno.tmp$FlyBaseID[jseq.All.geno.tmp$Sig], 
                                                    TRUE, Sig))

SSAV.geno <- merge(SSAV.geno, SDIU[,c("FlyBaseID", "Whole.SBGE.Osada", "tajD.S", "DoS")], by = "FlyBaseID", all = TRUE)
SSAV.geno <- merge(SSAV.geno, ASE, by = "FlyBaseID", all = TRUE)

all.candidates <- merge(all.candidates,  SDIU[,c("FlyBaseID", "Whole.SBGE.Osada", "tajD.S", "DoS")], by = "FlyBaseID", all = T)
all.candidates <- merge(all.candidates, ASE, by = "FlyBaseID", all = T)

##########



## plotting function
pointSEplot <- function(boot_dat, perm_dat, x_col, SBGE_cat = NA){
  # set y-axis value above each error bar
  if(!is.na(SBGE_cat)){
    y_count_FALSE <- boot_dat[!boot_dat$Sig,] %>% 
      dplyr::group_by(.[[SBGE_cat]]) %>%
      dplyr::summarise(max = q95 + 0.065) %>% # 0.02, 0.05
      dplyr::rename({{SBGE_cat}} := 1)
    y_count_TRUE <- boot_dat[boot_dat$Sig,] %>%
      dplyr::group_by(.[[SBGE_cat]]) %>%
      dplyr::summarise(max = q95 + 0.065) %>% # 0.02, 0.05
      dplyr::rename({{SBGE_cat}} := 1)
  }
  else{
    y_count_FALSE <- boot_dat[!boot_dat$Sig,] %>% 
      dplyr::summarise(max = q95 + 0.05)
    y_count_TRUE <- boot_dat[boot_dat$Sig,] %>% 
      dplyr:: summarise(max = q95 + 0.05)
    perm_dat <- as.data.frame(t(perm_dat))
    colnames(perm_dat) <- c("obs_diff", "n_TRUE", "n_FALSE", "pval", "Sig")
  }
  # join y-axis coordinate with text dataframes
  perm_dat$y_count_FALSE <- y_count_FALSE$max
  perm_dat$y_count_TRUE <- y_count_TRUE$max
  perm_dat$y_star <- max(perm_dat$y_count_FALSE, perm_dat$y_count_TRUE) + 0.05 # 0.2
  str(perm_dat)
  
  if(!is.na(SBGE_cat)){
    pointSE_plot <- ggplot(boot_dat) + 
      geom_point(aes_string(SBGE_cat, "q50", color = "Sig"), 
                 position = position_jitterdodge(jitter.width = 0.001), size = 6) +
      geom_errorbar(aes_string(SBGE_cat, "q50", color = "Sig"),
                    ymin = boot_dat$q05, ymax = boot_dat$q95, 
                    position = position_jitterdodge(jitter.width = 0.01), width = 0.5) +
      labs(x = "Sex-biased Gene Expression", # "omegaA_MK" = expression(italic("\u03c9A")[MK]); "alpha_MK" = expression(italic("\u03b1")[MK])
           y = print(x_col)) +
      geom_text(data = perm_dat %>% mutate(sig1 = if_else(Sig==1, "*", "")),
                aes_string(x =SBGE_cat, y = "y_star", label = "sig1"), color = "black", size = 16) +
      geom_text(data = perm_dat,aes_string(x =SBGE_cat, y = "y_count_FALSE", label = "n_FALSE"), 
                color = "black", size = 7.5, hjust = 1) +
      geom_text(data = perm_dat,aes_string(x =SBGE_cat, y = "y_count_TRUE", label = "n_TRUE"), 
                color = "black", size = 7.5, hjust = -0.8) +
      # change the labels accordingly if needed.
      scale_x_discrete(labels = c("Highly FB","Female-Biased", "Unbiased", "Male-Biased", "Highly MB")) +
      geom_vline(xintercept = c(seq(1.5, length(perm_dat[[SBGE_cat]]) )), color = "grey")
  }
  else{
    pointSE_plot <- ggplot(boot_dat) +
      geom_point(aes(x = Sig, y = q50, color = Sig), show.legend = FALSE,
                 position = position_jitterdodge(jitter.width = 0.001), size = 6) +
      geom_errorbar(aes(x = Sig, y = q50, color = Sig),
                    ymin = boot_dat$q05, ymax = boot_dat$q95, show.legend = FALSE,
                    position = position_jitterdodge(jitter.width = 0.01), width = 0.5) +
      labs(y = print(x_col)) +
      geom_text(data = perm_dat %>% mutate(sig1 = if_else(Sig==1, "*", "")),
                aes(x = 1.5, y = y_star+0.01, label = sig1), color = "black", size = 10) +
      geom_text(data = perm_dat,aes_string(x = "Sig", y = "y_count_FALSE", label = "n_FALSE"), 
                color = "black", size = 7.5) +
      geom_text(data = perm_dat,aes_string(x = "Sig", y = "y_count_TRUE", label = "n_TRUE"), 
                color = "black", size = 7.5, hjust = -6.5, nudge_x = 0.45) +
      scale_x_discrete(labels = c("Background","Candidates")) +
      labs(x = "\nSex-biased Gene Expression", y = print(x_col)) 
  }
  
  pointSE_plot <- pointSE_plot +
    scale_colour_manual(values = c("orange3", "forestgreen"), # "red3", "steelblue3", "#888888" # "purple3", "chartreuse3", "orange2"
                        labels = c("Background", "Candidates")) + # "Chr-2", "Chr-3", "X-Chr"
    guides(color = guide_legend(override.aes = list(shape = c(16, 16),
                                                    size = c(7.5, 7.5),
                                                    alpha = 1))) +
    # scale_y_continuous(limits = c(0,1)) +
    theme_classic() +
    theme(plot.title.position = c("panel"),
          legend.title = element_blank(),
          legend.position = c("bottom"),
          # legend.box.background = element_rect(),
          legend.text = element_text(size=30, margin = margin(10,0,10,0), color = "black", vjust = -0.2),
          axis.text.x = element_text(size=30, margin = margin(5,0,0,0), color = "black"),
          axis.text.y = element_text(size=20, margin = margin(0,5,0,0), color = "black"),
          axis.title.x = element_text(size=40, margin = margin(10,0,10,0), color = "black", vjust = -0.2),
          axis.title.y = element_text(size=40, margin = margin(0,10,0,0), color = "black"),
          plot.title = element_text(size=40, margin = margin(0,0,0,0), color = "black"),
          plot.margin = margin(6,6,6,6),
          panel.border = element_rect(colour = "black", fill=NA, size=1)
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



#### -------------- TajD ------------------
# D = d/SE(d), where d = pi - watterson's theta
# D < 0 implies directional selection or selective sweeps
# D > 0 implies balancing selection or relaxed purifying selection 
## all negative, suggesting that the population is undergoing adaptive evolution
# Candidate genes seem to have higher TajD.N estimates, suggesting that candidate and bg genes differ
# in some aspect of selection. This could be due to sexually antagonistic selection maybe resulting in 
# balancing selection or weaker purifying selection. Look at next section where we look at pi to test
# whether strength of purifying selection is different between candidates and non-candidates

# For all candidate genes vs non candidate genes
####### 
# use permutation and bootstrap functions in boot_permute.R
# only the AS genes, not significantly different
# only the DE genes also not significantly different
TwoPerm(jseq.All.geno.tmp, x_col = "tajD.S", groupBy = "Sig")

perm_TajD_N
perm_TajD <- TwoPerm(SSAV.geno[!is.na(SSAV.geno$Sig) & 
                                  !is.na(SSAV.geno$tajD.S),], x_col = "tajD.S", groupBy = "Sig")
boot_TajD_N
boot_TajD <- TwoBoot(SSAV.geno[!is.na(SSAV.geno$Sig) & 
                                 !is.na(SSAV.geno$tajD.S),], x_col = "tajD.S", 
                     groupBy = "Sig")
TajD_all <- pointSEplot(boot_dat = boot_TajD, perm_dat = perm_TajD, x_col = "tajD.S") + 
  coord_cartesian(ylim = c(-0.2, 0.1))

ggplot(SSAV.geno[!is.na(SSAV.geno$Sig) & 
                   !is.na(SSAV.geno$tajD.S),]) + 
  geom_point(aes(x = Sig, y = tajD.S, color = Sig), alpha = 0.5)

# by SBGE category 
# (according to Osada et al. as used in the Sing & Agrawal paper)
perm_TajD_SBGE_N 
perm_TajD_SBGE <- TwoPerm_SBGE(SSAV.geno[!is.na(SSAV.geno$Sig) & 
                                           !is.na(SSAV.geno$tajD.S) &
                                           !is.na(SSAV.geno$SBGE_comp),], 
                               x_col = "tajD.S", 
                               groupBy = "Sig", 
                               SBGE_cat = "SBGE_comp")
boot_TajD_SBGE_N
boot_TajD_SBGE <- TwoBoot_SBGE(SSAV.geno[!is.na(SSAV.geno$Sig) & 
                                           !is.na(SSAV.geno$tajD.S) &
                                           !is.na(SSAV.geno$SBGE_comp),],
                               x_col = "tajD.S", groupBy = "Sig", SBGE_cat = "SBGE_comp")
TajD_all_SBGE <- pointSEplot(boot_dat = boot_TajD_SBGE, perm_dat = perm_TajD_SBGE, 
                             x_col = "tajD.S", SBGE_cat = "SBGE_comp") + 
  scale_x_discrete(labels = c("Highly FB","Female-Biased", 
                              "Unbiased", "Male-Biased", "Highly MB")) + 
  coord_cartesian(ylim = c(-0.5, 0.5))


# Combined TajD plots All and separating SBGE 
perm_All_SBGE_TajD <- c("a.all", perm_TajD)
perm_All_SBGE_TajD <- rbind(perm_All_SBGE_TajD, perm_TajD_SBGE_ASE)

boot_All_SBGE_TajD <- boot_TajD
boot_All_SBGE_TajD$SBGE_comp <- "a.all"
boot_All_SBGE_TajD <- rbind(boot_All_SBGE_TajD, boot_TajD_SBGE_ASE)

# Exclude highly sex-biased genes
boot_All_SBGE_TajD <- boot_All_SBGE_TajD[boot_All_SBGE_TajD$SBGE_comp != "a.more.fbg" &
                                           boot_All_SBGE_TajD$SBGE_comp != "e.more.mbg" ,]
perm_All_SBGE_TajD <- perm_All_SBGE_TajD[perm_All_SBGE_TajD$SBGE_comp != "a.more.fbg" &
                                           perm_All_SBGE_TajD$SBGE_comp != "e.more.mbg" ,]


Fig6_main_N 
Fig6_main <- pointSEplot(boot_dat = boot_All_SBGE_TajD, 
                         perm_dat = perm_All_SBGE_TajD, 
                         x_col = "TajD.S",
                         SBGE_cat = "SBGE_comp") +
  scale_x_discrete(labels = c("All","Highly FB","Female-Biased", 
                              "Unbiased", "Male-Biased", "Highly MB")) +
  coord_cartesian(ylim = c(-0.2, 0.1)) +
  labs(y = "Tajima's D") +
  geom_vline(xintercept = 1.5, color = "black", size = 1.5) 

#######

# for male candidates
#######
perm_TajD_A.m <- TwoPerm(SSAV.geno[!is.na(SSAV.geno$A.m.Sig) & 
                                     !is.na(SSAV.geno$tajD.N),], 
                         x_col = "tajD.N", 
                         groupBy = "A.m.Sig")
boot_TajD_A.m <- TwoBoot(SSAV.geno[!is.na(SSAV.geno$A.m.Sig) & 
                                     !is.na(SSAV.geno$tajD.N),], 
                         x_col = "tajD.N", 
                         groupBy = "A.m.Sig")

# by SBGE category (Mishra et al.)
perm_TajD_SBGE_ASE_A.m <- TwoPerm_SBGE(SSAV.geno[!is.na(SSAV.geno$A.m.Sig) & 
                                                   !is.na(SSAV.geno$tajD.N) &
                                                   !is.na(SSAV.geno$SBGE_comp),], 
                                       x_col = "tajD.N", 
                                       groupBy = "A.m.Sig", 
                                       SBGE_cat = "SBGE_comp")
boot_TajD_SBGE_ASE_A.m <- TwoBoot_SBGE(SSAV.geno[!is.na(SSAV.geno$A.m.Sig) & 
                                                   !is.na(SSAV.geno$tajD.N) &
                                                   !is.na(SSAV.geno$SBGE_comp),],
                                       x_col = "tajD.N", 
                                       groupBy = "A.m.Sig", 
                                       SBGE_cat = "SBGE_comp")

## combined all and SBGE-separated plot
perm_All_SBGE_TajD_A.m <- c("a.all", perm_TajD_A.m)
perm_All_SBGE_TajD_A.m <- rbind(perm_All_SBGE_TajD_A.m, perm_TajD_SBGE_ASE_A.m)

boot_All_SBGE_TajD_A.m <- boot_TajD_A.m
boot_All_SBGE_TajD_A.m$SBGE_comp <- "a.all"
boot_All_SBGE_TajD_A.m <- rbind(boot_All_SBGE_TajD_A.m, boot_TajD_SBGE_ASE_A.m)

# Exclude highly sex-biased genes
boot_All_SBGE_TajD_A.m <- boot_All_SBGE_TajD_A.m[boot_All_SBGE_TajD_A.m$SBGE_comp != "a.more.fbg" &
                                                   boot_All_SBGE_TajD_A.m$SBGE_comp != "e.more.mbg" ,]
perm_All_SBGE_TajD_A.m <- perm_All_SBGE_TajD_A.m[perm_All_SBGE_TajD_A.m$SBGE_comp != "a.more.fbg" &
                                                   perm_All_SBGE_TajD_A.m$SBGE_comp != "e.more.mbg" ,]

Fig6A_suppl <- pointSEplot(boot_dat = boot_All_SBGE_TajD_A.m, 
                           perm_dat = perm_All_SBGE_TajD_A.m, 
                           x_col = "TajD.N",
                           SBGE_cat = "SBGE_comp") +
  scale_x_discrete(labels = c("All", "Highly FB","Female-Biased", 
                              "Unbiased", "Male-Biased", "Highly MB")) +
  coord_cartesian(ylim = c(-0.2, 0.1)) +
  labs(y = "Tajima's D") +
  geom_vline(xintercept = 1.5, color = "black", size = 1.5) 

#######

# for female candidates
#######
perm_TajD_A.f <- TwoPerm(SSAV.geno[!is.na(SSAV.geno$A.f.Sig) & 
                                     !is.na(SSAV.geno$tajD.N),], 
                         x_col = "tajD.N", 
                         groupBy = "A.f.Sig")
boot_TajD_A.f <- TwoBoot(SSAV.geno[!is.na(SSAV.geno$A.f.Sig) & 
                                     !is.na(SSAV.geno$tajD.N),], 
                         x_col = "tajD.N", 
                         groupBy = "A.f.Sig")

# by SBGE category (Mishra et al.)
perm_TajD_SBGE_ASE_A.f <- TwoPerm_SBGE(SSAV.geno[!is.na(SSAV.geno$A.f.Sig) & 
                                                   !is.na(SSAV.geno$tajD.N) &
                                                   !is.na(SSAV.geno$SBGE_comp),], 
                                       x_col = "tajD.N", 
                                       groupBy = "A.f.Sig", 
                                       SBGE_cat = "SBGE_comp")
boot_TajD_SBGE_ASE_A.f <- TwoBoot_SBGE(SSAV.geno[!is.na(SSAV.geno$A.f.Sig) & 
                                                   !is.na(SSAV.geno$tajD.N) &
                                                   !is.na(SSAV.geno$SBGE_comp),],
                                       x_col = "tajD.N", 
                                       groupBy = "A.f.Sig", 
                                       SBGE_cat = "SBGE_comp")

## combined all and SBGE-separated plot
perm_All_SBGE_TajD_A.f <- c("a.all", perm_TajD_A.f)
perm_All_SBGE_TajD_A.f <- rbind(perm_All_SBGE_TajD_A.f, perm_TajD_SBGE_ASE_A.f)

boot_All_SBGE_TajD_A.f <- boot_TajD_A.f
boot_All_SBGE_TajD_A.f$SBGE_comp <- "a.all"
boot_All_SBGE_TajD_A.f <- rbind(boot_All_SBGE_TajD_A.f, boot_TajD_SBGE_ASE_A.f)

# Exclude highly sex-biased genes
boot_All_SBGE_TajD_A.f <- boot_All_SBGE_TajD_A.f[boot_All_SBGE_TajD_A.f$SBGE_comp != "a.more.fbg" &
                                                   boot_All_SBGE_TajD_A.f$SBGE_comp != "e.more.mbg" ,]
perm_All_SBGE_TajD_A.f <- perm_All_SBGE_TajD_A.f[perm_All_SBGE_TajD_A.f$SBGE_comp != "a.more.fbg" &
                                                   perm_All_SBGE_TajD_A.f$SBGE_comp != "e.more.mbg" ,]


Fig6B_suppl <- pointSEplot(boot_dat = boot_All_SBGE_TajD_A.f, 
                           perm_dat = perm_All_SBGE_TajD_A.f, 
                           x_col = "TajD.N",
                           SBGE_cat = "SBGE_comp") +
  scale_x_discrete(labels = c("All","Highly FB", "Female-Biased", 
                              "Unbiased", "Male-Biased", "Highly MB")) +
  coord_cartesian(ylim = c(-0.2, 0.1)) +
  labs(y = "Tajima's D") +
  geom_vline(xintercept = 1.5, color = "black", size = 1.5) 

#######



#### -------------- relative piN ------------------
# nucleotide diversity at nonsynonymous sites vs at synonymous sites
# for all candidates
#######
perm_piNfract <- TwoPerm(SSAV.geno[!is.na(SSAV.geno$Sig) & 
                                 !is.na(SSAV.geno$piNfract),], x_col = "piNfract", groupBy = "Sig")
boot_piNfract <- TwoBoot(SSAV.geno[!is.na(SSAV.geno$Sig) & 
                                 !is.na(SSAV.geno$piNfract),], x_col = "piNfract", groupBy = "Sig")
piNfract_all <- pointSEplot(boot_dat = boot_piNfract, perm_dat = perm_piNfract, x_col = "piNfract") + 
  coord_cartesian(ylim = c(-0.1, 0.3))

# by SBGE category (according to Osada et al. as used in the Sing & Agrawal paper)
perm_piNfract_SBGE <- TwoPerm_SBGE(SSAV.geno[!is.na(SSAV.geno$Sig) & 
                                           !is.na(SSAV.geno$piNfract) & 
                                           !is.na(SSAV.geno$SBGEcat.body.Osada),], 
                               x_col = "piNfract", 
                               groupBy = "Sig", 
                               SBGE_cat = "SBGEcat.body.Osada")
boot_piNfract_SBGE <- TwoBoot_SBGE(SSAV.geno[!is.na(SSAV.geno$Sig) & 
                                           !is.na(SSAV.geno$piNfract) & 
                                           !is.na(SSAV.geno$SBGEcat.body.Osada),],
                               x_col = "piNfract", groupBy = "Sig", SBGE_cat = "SBGEcat.body.Osada")
piNfract_all_SBGE <- pointSEplot(boot_dat = boot_piNfract_SBGE, perm_dat = perm_piNfract_SBGE, 
                             x_col = "piNfract", SBGE_cat = "SBGEcat.body.Osada") + 
  scale_x_discrete(labels = c("Extreme FB","Strong FB","Female-Biased", 
                              "Unbiased", "Male-Biased", "Strong MB", "Extreme MB")) + 
  coord_cartesian(ylim = c(-0.2, 0.75))


## according to Mishra et al. SBGE bins
perm_piNfract_SBGE_ASE <- TwoPerm_SBGE(SSAV.geno[!is.na(SSAV.geno$Sig) & 
                                               !is.na(SSAV.geno$piNfract) & 
                                               !is.na(SSAV.geno$SBGE_comp),], 
                                   x_col = "piNfract", 
                                   groupBy = "Sig", 
                                   SBGE_cat = "SBGE_comp")
boot_piNfract_SBGE_ASE <- TwoBoot_SBGE(SSAV.geno[!is.na(SSAV.geno$Sig) & 
                                               !is.na(SSAV.geno$piNfract) & 
                                               !is.na(SSAV.geno$SBGE_comp),],
                                   x_col = "piNfract", groupBy = "Sig", SBGE_cat = "SBGE_comp")
piNfract_all_SBGE_ASE <- pointSEplot(boot_dat = boot_piNfract_SBGE_ASE, perm_dat = perm_piNfract_SBGE_ASE, 
                                 x_col = "piNfract", SBGE_cat = "SBGE_comp") + 
  scale_x_discrete(labels = c("Highly FB","Female-Biased", 
                              "Unbiased", "Male-Biased", "Highly MB")) + 
  coord_cartesian(ylim = c(-0.2, 0.75))



## combined all and SBGE-separated plot
perm_All_SBGE_piNfract <- c("a.all", perm_piNfract)
perm_All_SBGE_piNfract <- rbind(perm_All_SBGE_piNfract, perm_piNfract_SBGE_ASE)

boot_All_SBGE_piNfract <- boot_piNfract
boot_All_SBGE_piNfract$SBGE_comp <- "a.all"
boot_All_SBGE_piNfract <- rbind(boot_All_SBGE_piNfract, boot_piNfract_SBGE_ASE)

# Exclude highly sex-biased genes
boot_All_SBGE_piNfract <- boot_All_SBGE_piNfract[boot_All_SBGE_piNfract$SBGE_comp != "a.more.fbg" &
                                                   boot_All_SBGE_piNfract$SBGE_comp != "e.more.mbg" ,]
perm_All_SBGE_piNfract <- perm_All_SBGE_piNfract[perm_All_SBGE_piNfract$SBGE_comp != "a.more.fbg" &
                                                   perm_All_SBGE_piNfract$SBGE_comp != "e.more.mbg" ,]


Fig6B_main <- pointSEplot(boot_dat = boot_All_SBGE_piNfract, 
                           perm_dat = perm_All_SBGE_piNfract, 
                           x_col = "piNfract",
                           SBGE_cat = "SBGE_comp") +
  scale_x_discrete(labels = c("All", "Highly FB", "Female-Biased", 
                              "Unbiased", "Male-Biased", "Highly MB")) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(y = "piNfract") +
  geom_vline(xintercept = 1.5, color = "black", size = 1.5) 

#######

# for male candidates
#######
perm_piNfract_A.m <- TwoPerm(SSAV.geno[!is.na(SSAV.geno$A.m.Sig) & 
                                     !is.na(SSAV.geno$piNfract),], 
                         x_col = "piNfract", 
                         groupBy = "A.m.Sig")
boot_piNfract_A.m <- TwoBoot(SSAV.geno[!is.na(SSAV.geno$A.m.Sig) & 
                                     !is.na(SSAV.geno$piNfract),], 
                         x_col = "piNfract", 
                         groupBy = "A.m.Sig")

# by SBGE category (Mishra et al.)
perm_piNfract_SBGE_ASE_A.m <- TwoPerm_SBGE(SSAV.geno[!is.na(SSAV.geno$A.m.Sig) & 
                                                   !is.na(SSAV.geno$piNfract) &
                                                   !is.na(SSAV.geno$SBGE_comp),], 
                                       x_col = "piNfract", 
                                       groupBy = "A.m.Sig", 
                                       SBGE_cat = "SBGE_comp")
boot_piNfract_SBGE_ASE_A.m <- TwoBoot_SBGE(SSAV.geno[!is.na(SSAV.geno$A.m.Sig) & 
                                                       !is.na(SSAV.geno$piNfract) &
                                                       !is.na(SSAV.geno$SBGE_comp),],
                                       x_col = "piNfract", 
                                       groupBy = "A.m.Sig", 
                                       SBGE_cat = "SBGE_comp")

## combined all and SBGE-separated plot
perm_All_SBGE_piNfract_A.m <- c("a.all", perm_piNfract_A.m)
perm_All_SBGE_piNfract_A.m <- rbind(perm_All_SBGE_piNfract_A.m, perm_piNfract_SBGE_ASE_A.m)

boot_All_SBGE_piNfract_A.m <- boot_piNfract_A.m
boot_All_SBGE_piNfract_A.m$SBGE_comp <- "a.all"
colnames(boot_All_SBGE_piNfract_A.m)[1] <- "Sig"
boot_All_SBGE_piNfract_A.m <- rbind(boot_All_SBGE_piNfract_A.m, boot_piNfract_SBGE_ASE_A.m)

# Exclude highly sex-biased genes
boot_All_SBGE_piNfract_A.m <- boot_All_SBGE_piNfract_A.m[boot_All_SBGE_piNfract_A.m$SBGE_comp != "a.more.fbg" &
                                                           boot_All_SBGE_piNfract_A.m$SBGE_comp != "e.more.mbg" ,]
perm_All_SBGE_piNfract_A.m <- perm_All_SBGE_piNfract_A.m[perm_All_SBGE_piNfract_A.m$SBGE_comp != "a.more.fbg" &
                                                           perm_All_SBGE_piNfract_A.m$SBGE_comp != "e.more.mbg" ,]


Fig6C_suppl <- pointSEplot(boot_dat = boot_All_SBGE_piNfract_A.m, 
                           perm_dat = perm_All_SBGE_piNfract_A.m, 
                           x_col = "piNfract",
                           SBGE_cat = "SBGE_comp") +
  scale_x_discrete(labels = c("All","Highly FB", "Female-Biased", 
                              "Unbiased", "Male-Biased", "Highly MB")) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(y = "piNfract") +
  geom_vline(xintercept = 1.5, color = "black", size = 1.5) 

#######

# for female candidates
#######
perm_piNfract_A.f <- TwoPerm(SSAV.geno[!is.na(SSAV.geno$A.f.Sig) & 
                                     !is.na(SSAV.geno$piNfract),], 
                         x_col = "piNfract", 
                         groupBy = "A.f.Sig")
boot_piNfract_A.f <- TwoBoot(SSAV.geno[!is.na(SSAV.geno$A.f.Sig) & 
                                     !is.na(SSAV.geno$piNfract),], 
                         x_col = "piNfract", 
                         groupBy = "A.f.Sig")

# by SBGE category (Mishra et al.)
perm_piNfract_SBGE_ASE_A.f <- TwoPerm_SBGE(SSAV.geno[!is.na(SSAV.geno$A.f.Sig) & 
                                                   !is.na(SSAV.geno$piNfract) &
                                                   !is.na(SSAV.geno$SBGE_comp),], 
                                       x_col = "piNfract", 
                                       groupBy = "A.f.Sig", 
                                       SBGE_cat = "SBGE_comp")
boot_piNfract_SBGE_ASE_A.f <- TwoBoot_SBGE(SSAV.geno[!is.na(SSAV.geno$A.f.Sig) & 
                                                       !is.na(SSAV.geno$piNfract) &
                                                       !is.na(SSAV.geno$SBGE_comp),],
                                           x_col = "piNfract", 
                                           groupBy = "A.f.Sig", 
                                           SBGE_cat = "SBGE_comp")

## combined all and SBGE-separated plot
perm_All_SBGE_piNfract_A.f <- c("a.all", perm_piNfract_A.f)
perm_All_SBGE_piNfract_A.f <- rbind(perm_All_SBGE_piNfract_A.f, perm_piNfract_SBGE_ASE_A.f)

boot_All_SBGE_piNfract_A.f <- boot_piNfract_A.f
boot_All_SBGE_piNfract_A.f$SBGE_comp <- "a.all"
colnames(boot_All_SBGE_piNfract_A.f)[1] <- "Sig"
boot_All_SBGE_piNfract_A.f <- rbind(boot_All_SBGE_piNfract_A.f, boot_piNfract_SBGE_ASE_A.f)

# Exclude highly sex-biased genes
boot_All_SBGE_piNfract_A.f <- boot_All_SBGE_piNfract_A.f[boot_All_SBGE_piNfract_A.f$SBGE_comp != "a.more.fbg" &
                                                           boot_All_SBGE_piNfract_A.f$SBGE_comp != "e.more.mbg" ,]
perm_All_SBGE_piNfract_A.f <- perm_All_SBGE_piNfract_A.f[perm_All_SBGE_piNfract_A.f$SBGE_comp != "a.more.fbg" &
                                                           perm_All_SBGE_piNfract_A.f$SBGE_comp != "e.more.mbg" ,]


Fig6D_suppl <- pointSEplot(boot_dat = boot_All_SBGE_piNfract_A.f, 
                           perm_dat = perm_All_SBGE_piNfract_A.f, 
                           x_col = "piNfract",
                           SBGE_cat = "SBGE_comp") +
  scale_x_discrete(labels = c("All", "Highly FB","Female-Biased", 
                              "Unbiased", "Male-Biased", "Highly MB")) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(y = "piNfract") +
  geom_vline(xintercept = 1.5, color = "black", size = 1.5) 

#######


pdf(file = "~/Desktop/UofT/SSAV_RNA/Plots/TajDN_S.pdf",   # The directory you want to save the file in
    width = 18, # 18 The width of the plot in inches
    height = 15) # 15, 25 The height of the plot in inches
ggarrange(Fig6_main_N + theme(axis.text.x = element_blank(), axis.title.x = element_text(size = 5)) + labs(x=""),
          NA, Fig6_main,
          heights = c(1, 0.0005, 1),
          nrow = 3,
          common.legend = TRUE, legend = "bottom",
          font.label = list(size = 30), hjust = -0.01)

# ggarrange(Fig6A_suppl + theme(axis.text.x = element_blank(), axis.title.x = element_text(size = 5)) + labs(x=""),
#           NA,
#           Fig6C_suppl + theme(axis.text.x = element_blank(), axis.title.x = element_text(size = 5)) + labs(x=""),
#           NA,
#           Fig6B_suppl + theme(axis.text.x = element_blank(), axis.title.x = element_text(size = 5)) + labs(x=""),
#           NA,
#           Fig6D_suppl,
#           labels = c("A)", NA, NA, NA, "B)", NA, NA),
#           heights = c(1, 0.0005, 1, 0.0005, 1, 0.0005, 1),
#           nrow = 7,
#           common.legend = TRUE, legend = "bottom",
#           font.label = list(size = 30), hjust = -0.01)
dev.off()

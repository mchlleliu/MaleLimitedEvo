###################################
#
#                             Grieshop et al. 2023
#                             Author: Michelle Liu
#             DsRed experimental evolution - transcriptomics analysis
#                Comparing rMF between candidates vs non-candidates
# 
# 
###################################

rm(list=ls()) # Clears the environment
setwd("~/Desktop/UofT/SSAV_RNA/")

# load bootstrap & permutation functions from "boot_permute.R" first!!

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
y0 = 1 # tolerance of middle bins ### I DONT GET THIS
# xmid1 = (x1 + x2)/2

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

# load rMF calculations by Aneil (also in Singh & Agrawal 2023)
Aneil <- read.delim(file="~/Desktop/UofT/SSAV_RNA/Data/rmf.For.Karl.To.Use.v2.txt", header=TRUE)
colnames(Aneil) <- c("FlyBaseID", "rMF", "Var.F", "pVal.F", "Var.M", "pVal.M")
Aneil <- na.omit(Aneil) # There's an N/A for rMF


SSAV.geno_rMF <- merge(SSAV.geno, Aneil, by = "FlyBaseID", all = TRUE)
SSAV.geno_rMF <- SSAV.geno_rMF[!is.na(SSAV.geno_rMF$Sig) & !is.na(SSAV.geno_rMF$rMF),]

SSAV.geno_rMF_ASE <- merge(SSAV.geno_rMF, ASE, by = "FlyBaseID", all = T)
SSAV.geno_rMF_ASE <- SSAV.geno_rMF_ASE[!is.na(SSAV.geno_rMF_ASE$Sig) & 
                                 !is.na(SSAV.geno_rMF_ASE$rMF) & 
                                 !is.na(SSAV.geno_rMF_ASE$exp_SBGE_ase),]
sum(SSAV.geno$FlyBaseID[SSAV.geno$Sig] %in% Aneil$FlyBaseID)

##########

# split by expression magnitude
#######
SSAV.geno_rMF_ASE_exp <- merge(SSAV.geno_rMF_ASE, SDIU[,c("FlyBaseID", "log.avgExp.AdultLarva")])
SSAV.geno_rMF_ASE_exp <- SSAV.geno_rMF_ASE_exp[!is.na(SSAV.geno_rMF_ASE_exp$log.avgExp.AdultLarva) &
                                                 !is.na(SSAV.geno_rMF_ASE_exp$Sig) & 
                                                 !is.na(SSAV.geno_rMF_ASE_exp$rMF),]

summary(lm(rMF ~ Sig + log.avgExp.AdultLarva + SBGE_comp, data = SSAV.geno_rMF_ASE_exp))

q <- quantile(SSAV.geno_rMF_ASE_exp$log.avgExp.AdultLarva, na.rm = T, probs = c(0, 1/3, 2/3, 1))
SSAV.geno_rMF_ASE_exp_LOW <- SSAV.geno_rMF_ASE_exp[SSAV.geno_rMF_ASE_exp$log.avgExp.AdultLarva <= q[2],]
SSAV.geno_rMF_ASE_exp_MED <- SSAV.geno_rMF_ASE_exp[SSAV.geno_rMF_ASE_exp$log.avgExp.AdultLarva > q[2] &
                                                     SSAV.geno_rMF_ASE_exp$log.avgExp.AdultLarva <= q[3],]
SSAV.geno_rMF_ASE_exp_HI <- SSAV.geno_rMF_ASE_exp[SSAV.geno_rMF_ASE_exp$log.avgExp.AdultLarva > q[3],]
dim(SSAV.geno_rMF_ASE_exp_HI)
#######


## plotting function
pointSEplot <- function(boot_dat, perm_dat, x_col, SBGE_cat = NA){
  
  # set y-axis value above each error bar
  if(!is.na(SBGE_cat)){
    perm_dat$pval <- as.numeric(perm_dat$pval)
    y_count_FALSE <- boot_dat[!boot_dat$Sig,] %>% 
      dplyr::group_by(.[[SBGE_cat]]) %>%
      dplyr::summarise(max = q95 + 0.05) %>% # 0.05
      dplyr::rename({{SBGE_cat}} := 1)
    y_count_TRUE <- boot_dat[boot_dat$Sig,] %>%
      dplyr::group_by(.[[SBGE_cat]]) %>%
      dplyr::summarise(max = q95 + 0.05) %>% # 0.05
      dplyr::rename({{SBGE_cat}} := 1)
  }
  else{
    y_count_FALSE <- boot_dat[!boot_dat$Sig,] %>% 
      dplyr::summarise(max = q95 + 0.05)
    y_count_TRUE <- boot_dat[boot_dat$Sig,] %>% 
      dplyr::summarise(max = q95 + 0.05)
    perm_dat <- as.data.frame(t(perm_dat))
    colnames(perm_dat) <- c("obs_diff", "n_TRUE", "n_FALSE", "pval", "Sig")
  }
  # join y-axis coordinate with text dataframes
  perm_dat$y_count_FALSE <- y_count_FALSE$max
  perm_dat$y_count_TRUE <- y_count_TRUE$max
  perm_dat$y_star <- max(perm_dat$y_count_FALSE, perm_dat$y_count_TRUE) + 0.1
  perm_dat$y_star <- ifelse(perm_dat$y_star > 0.95, 0.85, perm_dat$y_star)
  str(perm_dat)
  
  if(!is.na(SBGE_cat)){
    pointSE_plot <- ggplot(boot_dat) + 
      geom_point(aes_string(SBGE_cat, "q50", color = "Sig"), 
                 position = position_jitterdodge(jitter.width = 0.001), size = 6) +
      geom_errorbar(aes_string(SBGE_cat, "q50", color = "Sig"),
                    ymin = boot_dat$q05, ymax = boot_dat$q95, size = 1, 
                    position = position_jitterdodge(jitter.width = 0.01), width = 0.5) +
      labs(x = "Sex-biased Gene Expression", # "omegaA_MK" = expression(italic("\u03c9A")[MK]); "alpha_MK" = expression(italic("\u03b1")[MK])
           y = print(x_col)) +
      geom_text(data = perm_dat %>% mutate(sig1 = if_else(pval < 0.0001, "***",
                                                          ifelse(pval < 0.01, "**", 
                                                                 ifelse(pval < 0.05, "*", "")))),
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
      geom_text(data = perm_dat %>% mutate(sig1 = if_else(pval < 0.0001, "***",
                                                          ifelse(pval < 0.01, "**", 
                                                                 ifelse(pval < 0.05, "*", "")))),
                aes(x = 1.5, y = y_star+0.01, label = sig1), color = "black", size = 10) +
      geom_text(data = perm_dat,aes_string(x = "Sig", y = "y_count_FALSE", label = "n_FALSE"), 
                color = "black", size = 7.5) +
      geom_text(data = perm_dat,aes_string(x = "Sig", y = "y_count_TRUE", label = "n_TRUE"), 
                color = "black", size = 7.5, hjust = -3) +
      scale_x_discrete(labels = c("BG","DE")) +
      labs(x = "\nSex-biased Gene Expression", y = print(x_col)) 
  }
  
  pointSE_plot <- pointSE_plot +
    scale_colour_manual(values = c("#AA4499", "#009E73"), # "red3", "steelblue3", "#888888" # "purple3", "chartreuse3", "orange2", "#E69F00", "#009E73"
                        labels = c("Background", "DE")) + # "Chr-2", "Chr-3", "X-Chr"
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
          axis.title.y = element_text(size=40, margin = margin(0,12,0,12), color = "black"),
          plot.title = element_text(size=40, margin = margin(0,0,0,0), color = "black"),
          plot.margin = margin(6,6,6,6),
          panel.border = element_rect(colour = "black", fill=NA, size=1)
    )
  return(pointSE_plot)
}


# mean rMF candidates vs non-candidates
########
# use the boot and permute functions in boot_permute.R
# run permutation to test for significant difference between candidates vs non-candidates
perm_rMF <- TwoPerm(SSAV.geno_rMF_ASE, x_col = "rMF", groupBy = "Sig")
# bootstrap 95% confidence interval
boot_rMF <- TwoBoot(SSAV.geno_rMF_ASE, x_col = "rMF", groupBy = "Sig")
# plot the CI and mean.
rMF_all <- pointSEplot(boot_dat = boot_rMF, perm_dat = perm_rMF, x_col = "rMF") + coord_cartesian(ylim = c(0,1))

# by expression magnitude
perm_rMF_LOW <- TwoPerm(SSAV.geno_rMF_ASE_exp_LOW, x_col = "rMF", groupBy = "Sig")
boot_rMF_LOW <- TwoBoot(SSAV.geno_rMF_ASE_exp_LOW, x_col = "rMF", groupBy = "Sig")
rMF_all_LOW <- pointSEplot(boot_dat = boot_rMF_LOW, perm_dat = perm_rMF_LOW, x_col = "rMF") + coord_cartesian(ylim = c(0,1))

perm_rMF_MED <- TwoPerm(SSAV.geno_rMF_ASE_exp_MED, x_col = "rMF", groupBy = "Sig")
boot_rMF_MED <- TwoBoot(SSAV.geno_rMF_ASE_exp_MED, x_col = "rMF", groupBy = "Sig")
rMF_all_MED <- pointSEplot(boot_dat = boot_rMF_MED, perm_dat = perm_rMF_MED, x_col = "rMF") + coord_cartesian(ylim = c(0,1))

perm_rMF_HI <- TwoPerm(SSAV.geno_rMF_ASE_exp_HI, x_col = "rMF", groupBy = "Sig")
boot_rMF_HI <- TwoBoot(SSAV.geno_rMF_ASE_exp_HI, x_col = "rMF", groupBy = "Sig")
rMF_all_HI <- pointSEplot(boot_dat = boot_rMF_HI, perm_dat = perm_rMF_HI, x_col = "rMF") + coord_cartesian(ylim = c(0,1))


exp_levels_rMF <- ggarrange(rMF_all + theme(axis.title.x = element_blank(),
                        axis.title.y = element_blank(), axis.text.x = element_blank()),
          ggarrange(rMF_all_LOW + theme(axis.title.x = element_blank(),
                                        axis.title.y = element_blank(),
                                        axis.text.x = element_blank()), 
                             rMF_all_MED + theme(axis.title.x = element_blank(),
                                                 axis.title.y = element_blank(),
                                                 axis.text.x = element_blank()), 
                             rMF_all_HI + theme(axis.title.x = element_blank(),
                                                axis.title.y = element_blank(),
                                                axis.text.x = element_blank()), 
                    ncol = 3),
          nrow = 2, heights = c(2, 1))


pdf(file = "~/Desktop/UofT/SSAV_RNA/Plots/test_rMF_exp.pdf",   # The directory you want to save the file in
    width = 14, # 12, 24, 20 The width of the plot in inches
    height = 14) # 10, 20, 13 The height of the plot in inches

annotate_figure(exp_levels_rMF + theme(plot.margin = margin(10,10,10,0)), left = text_grob("rMF", 
                                                                                           rot = 90, size = 40))
dev.off()


perm_rMF_simp <- TwoPerm(SSAV.geno_rMF_ASE[SSAV.geno_rMF_ASE$SBGE_comp != "a.more.fbg" & 
                                         SSAV.geno_rMF_ASE$SBGE_comp != "e.more.mbg",], 
                         x_col = "rMF", groupBy = "Sig")
# bootstrap 95% confidence interval
boot_rMF_simp <- TwoBoot(SSAV.geno_rMF_ASE[SSAV.geno_rMF_ASE$SBGE_comp != "a.more.fbg" & 
                                    SSAV.geno_rMF_ASE$SBGE_comp != "e.more.mbg",], x_col = "rMF", groupBy = "Sig")
# plot the CI and mean.
rMF_all_simp <- pointSEplot(boot_dat = boot_rMF_simp, perm_dat = perm_rMF_simp, x_col = "rMF")


# Separately for male and female candidates
perm_rMF_A.m <- TwoPerm(SSAV.geno_rMF_ASE[!is.na(SSAV.geno_rMF_ASE$A.m.Sig),], "rMF", "A.m.Sig")
boot_rMF_A.m <- TwoBoot(SSAV.geno_rMF_ASE[!is.na(SSAV.geno_rMF_ASE$A.m.Sig),], "rMF", "A.m.Sig")

perm_rMF_A.f <- TwoPerm(SSAV.geno_rMF_ASE[!is.na(SSAV.geno_rMF_ASE$A.f.Sig),], "rMF", "A.f.Sig")
boot_rMF_A.f <- TwoBoot(SSAV.geno_rMF_ASE[!is.na(SSAV.geno_rMF_ASE$A.f.Sig),], "rMF", "A.f.Sig")

########



## ------------- rMF by SBGE categories ----------------
#########
## for all candidates (combined male and female candidates)
# run permutation
perm_rMF_SBGE <- TwoPerm_SBGE(SSAV.geno_rMF_ASE, x_col = "rMF", groupBy = "Sig", SBGE_cat = "SBGE_comp")
# bootstrap 95% conf. interval
boot_rMF_SBGE <- TwoBoot_SBGE(SSAV.geno_rMF_ASE, x_col = "rMF", groupBy = "Sig", SBGE_cat = "SBGE_comp")
# plot
rMF_SBGE <- pointSEplot(boot_dat = boot_rMF_SBGE[boot_rMF_SBGE$SBGE_comp != "a.more.fbg" & 
                                                   boot_rMF_SBGE$SBGE_comp != "e.more.mbg",], 
                        perm_dat = perm_rMF_SBGE[perm_rMF_SBGE$SBGE_comp != "a.more.fbg" & 
                                                   perm_rMF_SBGE$SBGE_comp != "e.more.mbg",],
                        x_col = "rMF", SBGE_cat = "SBGE_comp") + 
  scale_x_discrete(labels = c("Female-Biased", "Unbiased", "Male-Biased"))


## only female candidates
# run permutation
perm_rMF_SBGE_A.f <- TwoPerm_SBGE(SSAV.geno_rMF_ASE[!is.na(SSAV.geno_rMF_ASE$A.f.Sig),], 
                                  x_col = "rMF", groupBy = "A.f.Sig", SBGE_cat = "SBGE_comp")
# bootstrap 95% conf. interval
boot_rMF_SBGE_A.f <- TwoBoot_SBGE(SSAV.geno_rMF_ASE[!is.na(SSAV.geno_rMF_ASE$A.f.Sig),],
                                  x_col = "rMF", groupBy = "A.f.Sig", SBGE_cat = "SBGE_comp")
# plot
rMF_SBGE_A.f <- pointSEplot(boot_dat = boot_rMF_SBGE_A.f[boot_rMF_SBGE_A.f$SBGE_comp != "a.more.fbg" & 
                                                           boot_rMF_SBGE_A.f$SBGE_comp != "e.more.mbg",], 
                            perm_dat = perm_rMF_SBGE_A.f[perm_rMF_SBGE_A.f$SBGE_comp != "a.more.fbg" & 
                                                           perm_rMF_SBGE_A.f$SBGE_comp != "e.more.mbg",], 
                        x_col = "rMF", SBGE_cat = "SBGE_comp") + # might just cut off the Highly SB genes 
  scale_x_discrete(labels = c("Female-Biased", "Unbiased", "Male-Biased")) 


## only male candidates
# run permutation
perm_rMF_SBGE_A.m <- TwoPerm_SBGE(SSAV.geno_rMF_ASE[!is.na(SSAV.geno_rMF_ASE$A.m.Sig),], 
                                  x_col = "rMF", groupBy = "A.m.Sig", SBGE_cat = "SBGE_comp")
# bootstrap 95% conf. interval
boot_rMF_SBGE_A.m <- TwoBoot_SBGE(SSAV.geno_rMF_ASE[!is.na(SSAV.geno_rMF_ASE$A.m.Sig),],
                                  x_col = "rMF", groupBy = "A.m.Sig", SBGE_cat = "SBGE_comp")
# plot
rMF_SBGE_A.m <- pointSEplot(boot_dat = boot_rMF_SBGE_A.m[boot_rMF_SBGE_A.m$SBGE_comp != "a.more.fbg" & 
                                                           boot_rMF_SBGE_A.m$SBGE_comp != "e.more.mbg",], 
                            perm_dat = perm_rMF_SBGE_A.m[perm_rMF_SBGE_A.m$SBGE_comp != "a.more.fbg" & 
                                                           perm_rMF_SBGE_A.m$SBGE_comp != "e.more.mbg",],
                            x_col = "rMF", SBGE_cat = "SBGE_comp") + # might just cut off the Highly SB genes
  theme(axis.text.x = element_blank())

#########



## ------------- plot All and SBGE together -----------------
# All genes
###### 
# prepare the combined bootstrap and permutation dataset
perm_All_SBGE_rMF <- c("a.all", perm_rMF)
perm_All_SBGE_rMF <- rbind(perm_All_SBGE_rMF, perm_rMF_SBGE)

boot_All_SBGE_rMF <- boot_rMF
boot_All_SBGE_rMF$SBGE_comp <- "a.all"
boot_All_SBGE_rMF <- rbind(boot_All_SBGE_rMF, boot_rMF_SBGE)

boot_All_SBGE_rMF_simp <- boot_All_SBGE_rMF[boot_All_SBGE_rMF$SBGE_comp != "a.more.fbg" &
                                         boot_All_SBGE_rMF$SBGE_comp != "e.more.mbg" ,]
perm_All_SBGE_rMF_simp <- perm_All_SBGE_rMF[perm_All_SBGE_rMF$SBGE_comp != "a.more.fbg" &
                                              perm_All_SBGE_rMF$SBGE_comp != "e.more.mbg" ,]


Fig4_main <- pointSEplot(boot_dat = boot_All_SBGE_rMF,
                         perm_dat = perm_All_SBGE_rMF, 
                         x_col = "rMF", SBGE_cat = "SBGE_comp") +
  coord_cartesian(ylim = c(-0.2,1)) +
  scale_x_discrete(labels = c("All", "Highly FB","Female-Biased", 
                              "Unbiased", "Male-Biased", "Highly MB")) +
  geom_vline(xintercept = 1.5, color = "black", size = 1.5) +
  ylab(expression(italic("r"["MF"])))

#######


# Male candidates
######
perm_All_SBGE_rMF_A.m <- c("a.all", perm_rMF_A.m)
perm_All_SBGE_rMF_A.m <- rbind(perm_All_SBGE_rMF_A.m, perm_rMF_SBGE_A.m)

boot_All_SBGE_rMF_A.m <- boot_rMF_A.m
colnames(boot_All_SBGE_rMF_A.m)[1] <- "Sig" 
boot_All_SBGE_rMF_A.m$SBGE_comp <- "a.all"
boot_All_SBGE_rMF_A.m <- rbind(boot_All_SBGE_rMF_A.m, boot_rMF_SBGE_A.m)
# boot_All_SBGE_rMF_A.m <- boot_All_SBGE_rMF_A.m[boot_All_SBGE_rMF_A.m$SBGE_comp != "a.more.fbg" &
#                                               boot_All_SBGE_rMF_A.m$SBGE_comp != "e.more.mbg" ,]
# perm_All_SBGE_rMF_A.m <- perm_All_SBGE_rMF_A.m[perm_All_SBGE_rMF_A.m$SBGE_comp != "a.more.fbg" &
#                                                  perm_All_SBGE_rMF_A.m$SBGE_comp != "e.more.mbg" ,]

Fig4A_suppl <- pointSEplot(boot_dat = boot_All_SBGE_rMF_A.m,
                           perm_dat = perm_All_SBGE_rMF_A.m, 
                           x_col = "rMF", SBGE_cat = "SBGE_comp") +
  scale_x_discrete(labels = c("All", "Highly FB", "Female-Biased", 
                              "Unbiased", "Male-Biased", "Highly MB")) +
  coord_cartesian(ylim = c(-0.2, 1)) +
  geom_vline(xintercept = 1.5, color = "black", size = 1.5) +
  ylab(expression(italic("r"["MF"])))
#######


# Female candidates
######
perm_All_SBGE_rMF_A.f <- c("a.all", perm_rMF_A.f)
perm_All_SBGE_rMF_A.f <- rbind(perm_All_SBGE_rMF_A.f, perm_rMF_SBGE_A.f)

boot_All_SBGE_rMF_A.f <- boot_rMF_A.f
colnames(boot_All_SBGE_rMF_A.f)[1] <- "Sig" 
boot_All_SBGE_rMF_A.f$SBGE_comp <- "a.all"
boot_All_SBGE_rMF_A.f <- rbind(boot_All_SBGE_rMF_A.f, boot_rMF_SBGE_A.f)

# exclude highly sex-biased genes
# boot_All_SBGE_rMF_A.f <- boot_All_SBGE_rMF_A.f[boot_All_SBGE_rMF_A.f$SBGE_comp != "a.more.fbg" &
#                                                  boot_All_SBGE_rMF_A.f$SBGE_comp != "e.more.mbg" ,]
# perm_All_SBGE_rMF_A.f <- perm_All_SBGE_rMF_A.f[perm_All_SBGE_rMF_A.f$SBGE_comp != "a.more.fbg" &
#                                                  perm_All_SBGE_rMF_A.f$SBGE_comp != "e.more.mbg" ,]

Fig4B_suppl <- pointSEplot(boot_dat = boot_All_SBGE_rMF_A.f,
                           perm_dat = perm_All_SBGE_rMF_A.f, 
                           x_col = "rMF", SBGE_cat = "SBGE_comp") +
  scale_x_discrete(labels = c("All", "Highly FB", "Female-Biased", 
                              "Unbiased", "Male-Biased", "Highly MB")) +
  coord_cartesian(ylim = c(-0.5,1)) +
  geom_vline(xintercept = 1.5, color = "black", size = 1.5) +
  ylab(expression(italic("r"["MF"])))
#######





pdf(file = "~/Desktop/UofT/SSAV_RNA/Plots/final_2/Fig3_main.pdf",   # The directory you want to save the file in
    width = 17, # The width of the plot in inches
    height = 9) # 9 18 The height of the plot in inches
# ggarrange(Fig4A_suppl + theme(axis.text.x = element_blank()) + labs(x= ""),
#           Fig4B_suppl + theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 25)),
#           labels = c("A)", "B)"),
#           nrow = 2, heights = c(1, 0.9),
#           common.legend = TRUE, legend = "bottom",
#           font.label = list(size = 30), hjust = -0.01)
Fig4_main + theme(axis.title.x = element_blank(), legend.text = element_text(size = 30, vjust = 1),
                  axis.text.x = element_text(size = 25))
dev.off()

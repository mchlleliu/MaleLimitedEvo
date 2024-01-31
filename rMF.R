###################################
#
#                             Grieshop et al. 2023
#                             Author: Michelle Liu
#             DsRed experimental evolution - transcriptomics analysis
#                Comparing rMF between candidates vs non-candidates
# 
# 
###################################

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

##########


## plotting function
pointSEplot <- function(boot_dat, perm_dat, x_col, SBGE_cat = NA){
  # set y-axis value above each error bar
  if(!is.na(SBGE_cat)){
    y_count_FALSE <- boot_dat[!boot_dat$Sig,] %>% 
      group_by(.[[SBGE_cat]]) %>%
      summarise(max = q95 + 0.05) %>% # 0.05
      rename({{SBGE_cat}} := 1)
    y_count_TRUE <- boot_dat[boot_dat$Sig,] %>%
      group_by(.[[SBGE_cat]]) %>%
      summarise(max = q95 + 0.05) %>% # 0.05
      rename({{SBGE_cat}} := 1)
  }
  else{
    y_count_FALSE <- boot_dat[!boot_dat$Sig,] %>% 
      summarise(max = q95 + 0.05)
    y_count_TRUE <- boot_dat[boot_dat$Sig,] %>% 
      summarise(max = q95 + 0.05)
    perm_dat <- as.data.frame(t(perm_dat))
    colnames(perm_dat) <- c("obs_diff", "n_TRUE", "n_FALSE", "pval", "Sig")
  }
  # join y-axis coordinate with text dataframes
  perm_dat$y_count_FALSE <- y_count_FALSE$max
  perm_dat$y_count_TRUE <- y_count_TRUE$max
  perm_dat$y_star <- max(perm_dat$y_count_FALSE, perm_dat$y_count_TRUE) + 0.2
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
perm_rMF <- TwoPerm(SSAV.geno_rMF, x_col = "rMF", groupBy = "Sig")
# bootstrap 95% confidence interval
boot_rMF <- TwoBoot(SSAV.geno_rMF, x_col = "rMF", groupBy = "Sig")
# plot the CI and mean.
rMF_all <- pointSEplot(boot_dat = boot_rMF, perm_dat = perm_rMF, x_col = "rMF")


# Separately for male and female candidates
perm_rMF_A.m <- TwoPerm(SSAV.geno_rMF[!is.na(SSAV.geno_rMF$A.m.Sig),], "rMF", "A.m.Sig")
boot_rMF_A.m <- TwoBoot(SSAV.geno_rMF[!is.na(SSAV.geno_rMF$A.m.Sig),], "rMF", "A.m.Sig")

perm_rMF_A.f <- TwoPerm(SSAV.geno_rMF[!is.na(SSAV.geno_rMF$A.f.Sig),], "rMF", "A.f.Sig")
boot_rMF_A.f <- TwoBoot(SSAV.geno_rMF[!is.na(SSAV.geno_rMF$A.f.Sig),], "rMF", "A.f.Sig")

########



## ------------- rMF by SBGE categories ----------------
#########
# prepare dataset
str(ASE) # load ASE data from Mishra et al. (look at External_data.R)
SSAV.geno_ASE <- merge(SSAV.geno_rMF, ASE, by = "FlyBaseID", all = T)
SSAV.geno_ASE <- SSAV.geno_ASE[!is.na(SSAV.geno_ASE$Sig) & 
                                 !is.na(SSAV.geno_ASE$rMF) & 
                                 !is.na(SSAV.geno_ASE$exp_SBGE_ase),]

## for all candidates (combined male and female candidates)
# run permutation
perm_rMF_SBGE <- TwoPerm_SBGE(SSAV.geno_ASE, x_col = "rMF", groupBy = "Sig", SBGE_cat = "SBGE_comp")
# bootstrap 95% conf. interval
boot_rMF_SBGE <- TwoBoot_SBGE(SSAV.geno_ASE, x_col = "rMF", groupBy = "Sig", SBGE_cat = "SBGE_comp")
# plot
rMF_SBGE <- pointSEplot(boot_dat = boot_rMF_SBGE[boot_rMF_SBGE$SBGE_comp != "a.more.fbg" & 
                                                   boot_rMF_SBGE$SBGE_comp != "e.more.mbg",], 
                        perm_dat = perm_rMF_SBGE[perm_rMF_SBGE$SBGE_comp != "a.more.fbg" & 
                                                   perm_rMF_SBGE$SBGE_comp != "e.more.mbg",],
                        x_col = "rMF", SBGE_cat = "SBGE_comp") + 
  scale_x_discrete(labels = c("Female-Biased", "Unbiased", "Male-Biased"))


## only female candidates
# run permutation
perm_rMF_SBGE_A.f <- TwoPerm_SBGE(SSAV.geno_ASE[!is.na(SSAV.geno_ASE$A.f.Sig),], 
                                  x_col = "rMF", groupBy = "A.f.Sig", SBGE_cat = "SBGE_comp")
# bootstrap 95% conf. interval
boot_rMF_SBGE_A.f <- TwoBoot_SBGE(SSAV.geno_ASE[!is.na(SSAV.geno_ASE$A.f.Sig),],
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
perm_rMF_SBGE_A.m <- TwoPerm_SBGE(SSAV.geno_ASE[!is.na(SSAV.geno_ASE$A.m.Sig),], 
                                  x_col = "rMF", groupBy = "A.m.Sig", SBGE_cat = "SBGE_comp")
# bootstrap 95% conf. interval
boot_rMF_SBGE_A.m <- TwoBoot_SBGE(SSAV.geno_ASE[!is.na(SSAV.geno_ASE$A.m.Sig),],
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


Fig4_main <- pointSEplot(boot_dat = boot_All_SBGE_rMF_simp,
                         perm_dat = perm_All_SBGE_rMF_simp, 
                         x_col = "rMF", SBGE_cat = "SBGE_comp") +
  coord_cartesian(ylim = c(0,1)) +
  scale_x_discrete(labels = c("All", "Female-Biased", 
                              "Unbiased", "Male-Biased")) +
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
boot_All_SBGE_rMF_A.m <- boot_All_SBGE_rMF_A.m[boot_All_SBGE_rMF_A.m$SBGE_comp != "a.more.fbg" &
                                              boot_All_SBGE_rMF_A.m$SBGE_comp != "e.more.mbg" ,]
perm_All_SBGE_rMF_A.m <- perm_All_SBGE_rMF_A.m[perm_All_SBGE_rMF_A.m$SBGE_comp != "a.more.fbg" &
                                                 perm_All_SBGE_rMF_A.m$SBGE_comp != "e.more.mbg" ,]

Fig4A_suppl <- pointSEplot(boot_dat = boot_All_SBGE_rMF_A.m,
                           perm_dat = perm_All_SBGE_rMF_A.m, 
                           x_col = "rMF", SBGE_cat = "SBGE_comp") +
  scale_x_discrete(labels = c("All", "Female-Biased", 
                              "Unbiased", "Male-Biased")) +
  coord_cartesian(ylim = c(0,1)) +
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
boot_All_SBGE_rMF_A.f <- boot_All_SBGE_rMF_A.f[boot_All_SBGE_rMF_A.f$SBGE_comp != "a.more.fbg" &
                                                 boot_All_SBGE_rMF_A.f$SBGE_comp != "e.more.mbg" ,]
perm_All_SBGE_rMF_A.f <- perm_All_SBGE_rMF_A.f[perm_All_SBGE_rMF_A.f$SBGE_comp != "a.more.fbg" &
                                                 perm_All_SBGE_rMF_A.f$SBGE_comp != "e.more.mbg" ,]

Fig4B_suppl <- pointSEplot(boot_dat = boot_All_SBGE_rMF_A.f,
                           perm_dat = perm_All_SBGE_rMF_A.f, 
                           x_col = "rMF", SBGE_cat = "SBGE_comp") +
  scale_x_discrete(labels = c("All", "Female-Biased", 
                              "Unbiased", "Male-Biased")) +
  coord_cartesian(ylim = c(0,1)) +
  geom_vline(xintercept = 1.5, color = "black", size = 1.5) +
  ylab(expression(italic("r"["MF"])))
#######





pdf(file = "~/Desktop/UofT/SSAV_RNA/Plots/finals/Fig4_suppl.pdf",   # The directory you want to save the file in
    width = 16, # 16 The width of the plot in inches
    height = 20) # 10 The height of the plot in inches
ggarrange(Fig4A_suppl + theme(axis.text.x = element_blank()) + labs(x= ""),
          Fig4B_suppl,
          labels = c("A)", "B)"),
          nrow = 2,
          common.legend = TRUE, legend = "bottom",
          font.label = list(size = 30), hjust = -0.01)
# Fig4_main + theme(axis.title.x = element_blank(), legend.text = element_text(size = 22.5))
dev.off()


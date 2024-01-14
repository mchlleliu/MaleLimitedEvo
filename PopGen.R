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
SSAV.geno <- merge(SSAV.geno, ASE, by = "FlyBaseID", all = TRUE)
##########

## plotting function
pointSEplot <- function(boot_dat, perm_dat, x_col, SBGE_cat = NA){
  # set y-axis value above each error bar
  if(!is.na(SBGE_cat)){
    y_count_FALSE <- boot_dat[!boot_dat$Sig,] %>% 
      group_by(.[[SBGE_cat]]) %>%
      summarise(max = q95 + 0.065) %>% # 0.02, 0.05
      rename({{SBGE_cat}} := 1)
    y_count_TRUE <- boot_dat[boot_dat$Sig,] %>%
      group_by(.[[SBGE_cat]]) %>%
      summarise(max = q95 + 0.065) %>% # 0.02, 0.05
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
perm_TajD <- TwoPerm(SSAV.geno[!is.na(SSAV.geno$Sig) & 
                                !is.na(SSAV.geno$tajD.N),], x_col = "tajD.N", groupBy = "Sig")
boot_TajD <- TwoBoot(SSAV.geno[!is.na(SSAV.geno$Sig) & 
                                !is.na(SSAV.geno$tajD.N),], x_col = "tajD.N", groupBy = "Sig")
TajD_all <- pointSEplot(boot_dat = boot_TajD, perm_dat = perm_TajD, x_col = "tajD.N") + 
  coord_cartesian(ylim = c(-0.2, 0.1))

# by SBGE category 
# (according to Osada et al. as used in the Sing & Agrawal paper)
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
TajD_all_SBGE <- pointSEplot(boot_dat = boot_TajD_SBGE, perm_dat = perm_TajD_SBGE, 
                             x_col = "tajD.N", SBGE_cat = "SBGEcat.body.Osada") + 
  scale_x_discrete(labels = c("Extreme FB","Strong FB","Female-Biased", 
                              "Unbiased", "Male-Biased", "Strong MB", "Extreme MB")) + coord_cartesian(ylim = c(-0.5, 0.5))


# by SBGE category according to Mishra et al. population
# merge datasets
SSAV.geno <- merge(SSAV.geno, ASE, by = "FlyBaseID", all = T)
perm_TajD_SBGE_ASE <- TwoPerm_SBGE(SSAV.geno[!is.na(SSAV.geno$Sig) & 
                                               !is.na(SSAV.geno$tajD.N) & 
                                               !is.na(SSAV.geno$SBGE_comp),], 
                                   x_col = "tajD.N", 
                                   groupBy = "Sig", 
                                   SBGE_cat = "SBGE_comp")
boot_TajD_SBGE_ASE <- TwoBoot_SBGE(SSAV.geno[!is.na(SSAV.geno$Sig) & 
                                               !is.na(SSAV.geno$tajD.N) & 
                                               !is.na(SSAV.geno$SBGE_comp),],
                                   x_col = "tajD.N", groupBy = "Sig", SBGE_cat = "SBGE_comp")
TajD_all_SBGE_ASE <- pointSEplot(boot_dat = boot_TajD_SBGE_ASE, perm_dat = perm_TajD_SBGE_ASE, 
                                 x_col = "tajD.N", SBGE_cat = "SBGE_comp") + 
  scale_x_discrete(labels = c("Extreme FB","Strong FB","Female-Biased", 
                              "Unbiased", "Male-Biased", "Strong MB", "Extreme MB")) + coord_cartesian(ylim = c(-0.5, 0.5))

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


Fig6_main <- pointSEplot(boot_dat = boot_All_SBGE_TajD, 
                         perm_dat = perm_All_SBGE_TajD, 
                         x_col = "TajD.N",
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


pdf(file = "~/Desktop/UofT/SSAV_RNA/Plots/finals/Fig6_main.pdf",   # The directory you want to save the file in
    width = 18, # 18 The width of the plot in inches
    height = 15) # 15, 25 The height of the plot in inches
ggarrange(Fig6_main + theme(axis.text.x = element_blank(), axis.title.x = element_text(size = 5)) + labs(x=""),
          NA, Fig6B_main,
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

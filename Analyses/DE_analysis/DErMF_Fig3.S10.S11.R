###################################
#
#                             Grieshop et al. 2025
#                             Author: Michelle Liu
#             DsRed experimental evolution - transcriptomics analysis
#                Comparing rMF between candidates vs non-candidates
# 
# 
###################################

rm(list=ls()) # Clears the environment
setwd("~/Desktop/UofT/SSAV_RNA/")

# load bootstrap & permutation functions
source("boot_permute.R")

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
source("Mishra_et.al_SBGE.R")


## prepare dataset
##########
source("DE_Plotting_data.R")

# load rMF calculations in Singh & Agrawal 2023
# (estimated from the Huang et al. 2015 data set) 
rMF <- read.delim(file="~/Desktop/UofT/SSAV_RNA/Data/rmf.For.Karl.To.Use.v2.txt", header=TRUE)
colnames(rMF) <- c("FlyBaseID", "rMF", "Var.F", "pVal.F", "Var.M", "pVal.M")
rMF <- na.omit(rMF) # There's an N/A for rMF

# combine rMF info with Exp. population dataset
SSAV.geno_rMF <- merge(SSAV.geno, rMF, by = "FlyBaseID", all = TRUE)
SSAV.geno_rMF <- SSAV.geno_rMF[!is.na(SSAV.geno_rMF$Sig) & !is.na(SSAV.geno_rMF$rMF),]

# how many genes in our dataset has an estimated rMF value? 
sum(SSAV.geno$FlyBaseID[SSAV.geno$Sig] %in% rMF$FlyBaseID)
rm(rMF) # clear this object since it has already been incorporated elsewhere
##########



## plotting function
pointSEplot <- function(boot_dat, perm_dat, x_col, SBGE_cat = NA){
  
  # set y-axis value above each error bar
  if(!is.na(SBGE_cat)){ # stratified by SBGE category
    perm_dat$pval <- as.numeric(perm_dat$pval)
    y_count_FALSE <- boot_dat[!boot_dat$Sig,] %>% 
      dplyr::group_by(.[[SBGE_cat]]) %>%
      dplyr::summarise(y_count_FALSE = q95 + 0.05) %>% # 0.05
      dplyr::rename({{SBGE_cat}} := 1)
    
    y_count_TRUE <- boot_dat[boot_dat$Sig,] %>%
      dplyr::group_by(.[[SBGE_cat]]) %>%
      dplyr::summarise(y_count_TRUE = q95 + 0.05) %>% # 0.05
      dplyr::rename({{SBGE_cat}} := 1)
    
    
    # join y-axis coordinate with text data frames (comment in if plotting the stratified by exp plot)
    perm_dat <- merge(perm_dat, y_count_FALSE, by = SBGE_cat)
    perm_dat <- merge(perm_dat, y_count_TRUE, by = SBGE_cat)
    
  }
  else{ # not stratified by SBGE category
    y_count_FALSE <- boot_dat[!boot_dat$Sig,] %>% 
      dplyr::summarise(y_count_FALSE = q95 + 0.05)
    y_count_TRUE <- boot_dat[boot_dat$Sig,] %>% 
      dplyr::summarise(y_count_TRUE = q95 + 0.05)
    perm_dat <- as.data.frame(t(perm_dat))
    colnames(perm_dat) <- c("obs_diff", "n_TRUE", "n_FALSE", "pval", "Sig")
    
    perm_dat <- merge(perm_dat, y_count_FALSE)
    perm_dat <- merge(perm_dat, y_count_TRUE)
    
  }
  
  perm_dat$y_star <- max(perm_dat$y_count_FALSE, perm_dat$y_count_TRUE) + 0.1
  perm_dat$y_star <- ifelse(perm_dat$y_star > 0.95, 0.85, perm_dat$y_star)
  str(perm_dat)
  
  # plot with grouping
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
  
  # regular plot (without grouping)
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
  
  
  # general theme settings
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
# run permutation to test for significant difference between candidates vs non-candidates
perm_rMF <- TwoPerm(SSAV.geno_rMF, x_col = "rMF", groupBy = "Sig")
# bootstrap 95% confidence interval
boot_rMF <- TwoBoot(SSAV.geno_rMF, x_col = "rMF", groupBy = "Sig")
# plot the CI and mean.
rMF_all <- pointSEplot(boot_dat = boot_rMF, 
                       perm_dat = perm_rMF, 
                       x_col = "rMF") + coord_cartesian(ylim = c(0,1))


# exclude highly sex-biased genes
perm_rMF_simp <- TwoPerm(SSAV.geno_rMF[SSAV.geno_rMF$SBGE_comp != "a.more.fbg" & 
                                         SSAV.geno_rMF$SBGE_comp != "e.more.mbg",], 
                         x_col = "rMF", groupBy = "Sig")
# bootstrap 95% confidence interval
boot_rMF_simp <- TwoBoot(SSAV.geno_rMF[SSAV.geno_rMF$SBGE_comp != "a.more.fbg" & 
                                        SSAV.geno_rMF$SBGE_comp != "e.more.mbg",], x_col = "rMF", groupBy = "Sig")
# plot the CI and mean.
rMF_all_simp <- pointSEplot(boot_dat = boot_rMF_simp, perm_dat = perm_rMF_simp, x_col = "rMF")


# Separately for male and female candidates
perm_rMF_A.m <- TwoPerm(SSAV.geno_rMF[!is.na(SSAV.geno_rMF$A.m.Sig),], "rMF", "A.m.Sig")
boot_rMF_A.m <- TwoBoot(SSAV.geno_rMF[!is.na(SSAV.geno_rMF$A.m.Sig),], "rMF", "A.m.Sig")

perm_rMF_A.f <- TwoPerm(SSAV.geno_rMF[!is.na(SSAV.geno_rMF$A.f.Sig),], "rMF", "A.f.Sig")
boot_rMF_A.f <- TwoBoot(SSAV.geno_rMF[!is.na(SSAV.geno_rMF$A.f.Sig),], "rMF", "A.f.Sig")

########


## ------------- rMF by SBGE categories ----------------
## for all candidates (combined male and female candidates)
# run permutation
perm_rMF_SBGE <- TwoPerm_SBGE(SSAV.geno_rMF, x_col = "rMF", groupBy = "Sig", SBGE_cat = "SBGE_comp")
# bootstrap 95% conf. interval
boot_rMF_SBGE <- TwoBoot_SBGE(SSAV.geno_rMF, x_col = "rMF", groupBy = "Sig", SBGE_cat = "SBGE_comp")
# plot
rMF_SBGE <- pointSEplot(boot_dat = boot_rMF_SBGE, 
                        perm_dat = perm_rMF_SBGE,
                        x_col = "rMF", SBGE_cat = "SBGE_comp") + 
  coord_cartesian(ylim = c(-0.25, 1))


## only female candidates
# run permutation
perm_rMF_SBGE_A.f <- TwoPerm_SBGE(SSAV.geno_rMF[!is.na(SSAV.geno_rMF$A.f.Sig),], 
                                  x_col = "rMF", groupBy = "A.f.Sig", SBGE_cat = "SBGE_comp")
# bootstrap 95% conf. interval
boot_rMF_SBGE_A.f <- TwoBoot_SBGE(SSAV.geno_rMF[!is.na(SSAV.geno_rMF$A.f.Sig),],
                                  x_col = "rMF", groupBy = "A.f.Sig", SBGE_cat = "SBGE_comp")
# plot
rMF_SBGE_A.f <- pointSEplot(boot_dat = boot_rMF_SBGE_A.f, 
                            perm_dat = perm_rMF_SBGE_A.f, 
                        x_col = "rMF", SBGE_cat = "SBGE_comp")


## only male candidates
# run permutation
perm_rMF_SBGE_A.m <- TwoPerm_SBGE(SSAV.geno_rMF[!is.na(SSAV.geno_rMF$A.m.Sig),], 
                                  x_col = "rMF", groupBy = "A.m.Sig", SBGE_cat = "SBGE_comp")
# bootstrap 95% conf. interval
boot_rMF_SBGE_A.m <- TwoBoot_SBGE(SSAV.geno_rMF[!is.na(SSAV.geno_rMF$A.m.Sig),],
                                  x_col = "rMF", groupBy = "A.m.Sig", SBGE_cat = "SBGE_comp")
# plot
rMF_SBGE_A.m <- pointSEplot(boot_dat = boot_rMF_SBGE_A.m, 
                            perm_dat = perm_rMF_SBGE_A.m,
                            x_col = "rMF", SBGE_cat = "SBGE_comp") +
  theme(axis.text.x = element_blank()) # remove x-axis labels for final plot

#########


## ------------- plot All and SBGE together -----------------
# All genes
###### 
# prepare the combined bootstrap and permutation dataset
# perm_All_SBGE_rMF is a data frame with:
#   Sig: DE status
#   q05: lower confidence interval from 95% bootstrap
#   q50: point estimate of the mean
#   q95: upper confidence interval from 95% bootstrap
#   SBGE_comp: analysis group (all, or by SBGE categories)
perm_All_SBGE_rMF <- c("a.all", perm_rMF)
perm_All_SBGE_rMF <- rbind(perm_All_SBGE_rMF, perm_rMF_SBGE)

boot_All_SBGE_rMF <- boot_rMF
boot_All_SBGE_rMF$SBGE_comp <- "a.all"
boot_All_SBGE_rMF <- rbind(boot_All_SBGE_rMF, boot_rMF_SBGE)


Fig3_main <- pointSEplot(boot_dat = boot_All_SBGE_rMF,
                         perm_dat = perm_All_SBGE_rMF, 
                         x_col = "rMF", SBGE_cat = "SBGE_comp") +
  coord_cartesian(ylim = c(-0.2,1)) +
  scale_x_discrete(labels = c("All", "Highly FB","Female-Biased", 
                              "Unbiased", "Male-Biased", "Highly MB")) +
  geom_vline(xintercept = 1.5, color = "black", size = 1.5) +
  ylab(expression(italic("r"["MF"])))

#######


# Considering DE genes from male and female comparisons separately (Fig. S10):

# Male candidates
######
perm_All_SBGE_rMF_A.m <- c("a.all", perm_rMF_A.m)
perm_All_SBGE_rMF_A.m <- rbind(perm_All_SBGE_rMF_A.m, perm_rMF_SBGE_A.m)

boot_All_SBGE_rMF_A.m <- boot_rMF_A.m
colnames(boot_All_SBGE_rMF_A.m)[1] <- "Sig" 
boot_All_SBGE_rMF_A.m$SBGE_comp <- "a.all"
boot_All_SBGE_rMF_A.m <- rbind(boot_All_SBGE_rMF_A.m, boot_rMF_SBGE_A.m)

Fig5A_suppl <- pointSEplot(boot_dat = boot_All_SBGE_rMF_A.m,
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

Fig5B_suppl <- pointSEplot(boot_dat = boot_All_SBGE_rMF_A.f,
                           perm_dat = perm_All_SBGE_rMF_A.f, 
                           x_col = "rMF", SBGE_cat = "SBGE_comp") +
  scale_x_discrete(labels = c("All", "Highly FB", "Female-Biased", 
                              "Unbiased", "Male-Biased", "Highly MB")) +
  coord_cartesian(ylim = c(-0.5,1)) +
  geom_vline(xintercept = 1.5, color = "black", size = 1.5) +
  ylab(expression(italic("r"["MF"])))
#######




# Save figures
png(file = "~/Desktop/UofT/SSAV_RNA/Plots/revision/Fig_S10_rMF3Mb.png",   # The directory you want to save the file in
    width = 15, # 17 The width of the plot in inches
    height = 16, # 9 18 The height of the plot in inches
    units = "in", res = 300)

# Suppl. Figure S10
ggarrange(Fig5A_suppl + theme(axis.text.x = element_blank()) + labs(x= ""),
          Fig5B_suppl + theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 25)),
          labels = c("A)", "B)"),
          nrow = 2, heights = c(1, 0.9),
          common.legend = TRUE, legend = "bottom",
          font.label = list(size = 30), hjust = -0.01)

# Fig. 3
# Fig3_main + theme(axis.title.x = element_blank(), legend.text = element_text(size = 30, vjust = 1),
#                   axis.text.x = element_text(size = 25))


dev.off()



# (Figure S11)
# ------------- rMF by exp. magnitude ----------------
# expression magnitude was calculated in Singh & Agrawal 2023, using data from SEBIDA
SinghAgrawal <- read.csv(file="~/Desktop/UofT/SSAV_RNA/Data/SBGEandSSSdataForMBE.csv", sep=",", header=TRUE)
colnames(SinghAgrawal)[2] <- "FlyBaseID"

# save this in another object for now...
SSAV.geno_rMF_exp <- merge(SSAV.geno_rMF, SinghAgrawal[,c("FlyBaseID", "log.avgExp.AdultLarva")], by = "FlyBaseID")
SSAV.geno_rMF_exp <- SSAV.geno_rMF_exp[!is.na(SSAV.geno_rMF_exp$log.avgExp.AdultLarva) &
                                         !is.na(SSAV.geno_rMF_exp$Sig) & 
                                         !is.na(SSAV.geno_rMF_exp$rMF),]

rm(SinghAgrawal) # remove this data set from env

# 1.
# Linear model testing whether SA candidate status is still a significant predictor of rMF,
# controlling for rMF variation caused by expression level
summary(lm(rMF ~ Sig + log.avgExp.AdultLarva + SBGE_comp, data = SSAV.geno_rMF_exp))


# 2.
# Another way to test for this is to split the analysis based on exp levels
# define expression level category (Low/Medium/High)
q <- quantile(SSAV.geno_rMF_exp$log.avgExp.AdultLarva, na.rm = T, probs = c(0, 1/3, 2/3, 1))
SSAV.geno_rMF_exp <- SSAV.geno_rMF_exp %>%
  mutate(exp = ifelse(log.avgExp.AdultLarva <= q[2], "b.Low",
                      ifelse(log.avgExp.AdultLarva > q[2] & log.avgExp.AdultLarva <= q[3], "c.Medium", 
                             "d.High")))


# remember to call the bootstrapping and permutation functions in "boot_permute.R"
# NOTE: I made these functions originally to stratify based on SBGE categories (hence, "SBGE_cat"), 
#       but they can also be used for grouping other categories, as long as you specify the
#       column name which contains the categorical levels
# permute the two groups (DE/Non-DE)
perm_rMF_exp <- TwoPerm_SBGE(SSAV.geno_rMF_exp, x_col = "rMF", groupBy = "Sig", SBGE_cat = "exp")
# bootstrap 95% conf. interval
boot_rMF_exp <- TwoBoot_SBGE(SSAV.geno_rMF_exp, x_col = "rMF", groupBy = "Sig", SBGE_cat = "exp")

# combine with unstratified analysis result
perm_All_exp_rMF <- c("a.all", perm_rMF)
perm_All_exp_rMF <- rbind(perm_All_exp_rMF, perm_rMF_exp)

boot_All_exp_rMF <- boot_rMF
boot_All_exp_rMF$exp <- "a.all"
boot_All_exp_rMF <- rbind(boot_All_exp_rMF, boot_rMF_exp)


# plot Figure S11.
Fig_S11 <- pointSEplot(boot_dat = boot_All_exp_rMF,
                      perm_dat = perm_All_exp_rMF, 
                      x_col = "rMF", SBGE_cat = "exp") +
  coord_cartesian(ylim = c(0,1)) +
  scale_x_discrete(labels = c("All", "Low", "Medium", "High")) +
  geom_vline(xintercept = 1.5, color = "black", size = 1.5) +
  ylab(expression(italic("r"["MF"]))) +
  theme(axis.title.x = element_blank(), 
        legend.text = element_text(size = 30, vjust = 1),
        axis.text.x = element_text(size = 25))


# save plot
png(file = "~/Desktop/UofT/SSAV_RNA/Plots/final_2/rMF_Fig_S11.png",   # The directory you want to save the file in
    width = 12, # 12, 24, 20 The width of the plot in inches
    height = 7, # 10, 20, 13 The height of the plot in inches
    units = "in", res = 300)
    
Fig_S6

dev.off()
########

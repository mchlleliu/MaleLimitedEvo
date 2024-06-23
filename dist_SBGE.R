###################################
#
#                             Grieshop et al. 2023
#                             Author: Michelle Liu
#             DsRed experimental evolution - transcriptomics analysis
#                  Distribution of SA candidates by sex-bias
# 
# 
###################################

# This uses Mishra et al.'s data, 
# Or can be used with any data which have log2FC value of SBGE
# Bins for SBGE categories should be defined first. 
#     -> See External_data.R under "Get Mishra et al.'s data" for how to prepare the external dataset.
rm(list=ls()) # Clears the environment
setwd("~/Desktop/UofT/SSAV_RNA/")

# packages required
#########
library(tidyr)
library(plyr)
library(dplyr)
library(broom)
library(ggplot2)
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

# Prepare plotting dataset 
# (Run this first before anything else!!, make sure you have the correct files/path to them)
#########
# load results if not loaded in env.
A.f.geno <- read.delim("Results/A.f.geno_candidates.tsv")
A.m.geno <- read.delim("Results/A.m.geno_candidates.tsv")
SSAV.geno <- read.delim("Results/All.geno_candidates.tsv")


# include SBGE categories (using Mishra et al. dataset. Look at External_data.R)
A.m.geno_ASE <- merge(A.m.geno, ASE, by = "FlyBaseID", all = TRUE)
A.m.geno_ASE <- A.m.geno_ASE[!is.na(A.m.geno_ASE$Sig) & !is.na(A.m.geno_ASE$exp_SBGE_ase),]
A.m.geno_ASE$SBGE_comp <- as.factor(A.m.geno_ASE$SBGE_comp)
A.m.geno_ASE$SBGE_simp <- as.factor(A.m.geno_ASE$SBGE_simp)
str(A.m.geno_ASE)

A.f.geno_ASE <- merge(A.f.geno, ASE, by = "FlyBaseID", all = TRUE)
A.f.geno_ASE <- A.f.geno_ASE[!is.na(A.f.geno_ASE$Sig) & !is.na(A.f.geno_ASE$exp_SBGE_ase),]
A.f.geno_ASE$SBGE_comp <- as.factor(A.f.geno_ASE$SBGE_comp)
A.f.geno_ASE$SBGE_simp <- as.factor(A.f.geno_ASE$SBGE_simp)
str(A.f.geno_ASE)

# Genes present in both SSAV males and SSAV females data
SSAV.geno_ASE <- merge(SSAV.geno, ASE, by = "FlyBaseID", all = TRUE)
SSAV.geno_ASE <- SSAV.geno_ASE[!is.na(SSAV.geno_ASE$Sig) & !is.na(SSAV.geno_ASE$exp_SBGE_ase),]
SSAV.geno_ASE$SBGE_comp <- as.factor(SSAV.geno_ASE$SBGE_comp)
SSAV.geno_ASE$SBGE_simp <- as.factor(SSAV.geno_ASE$SBGE_simp)
str(SSAV.geno_ASE)


SSAV.geno_ASE_exp <- merge(SSAV.geno_ASE, SDIU[,c("FlyBaseID", "log.avgExp.AdultLarva")])
SSAV.geno_ASE_exp <- SSAV.geno_ASE_exp[!is.na(SSAV.geno_ASE_exp$log.avgExp.AdultLarva) &
                                   !is.na(SSAV.geno_ASE_exp$Sig) & 
                                    !is.na(SSAV.geno_ASE_exp$exp_SBGE_ase),]

q <- quantile(SSAV.geno_ASE_exp$log.avgExp.AdultLarva, na.rm = T, probs = c(0, 1/3, 2/3, 1))
SSAV.geno_ASE_exp_LOW <- SSAV.geno_ASE_exp[SSAV.geno_ASE_exp$log.avgExp.AdultLarva <= q[2],]
SSAV.geno_ASE_exp_MED <- SSAV.geno_ASE_exp[SSAV.geno_ASE_exp$log.avgExp.AdultLarva > q[2] &
                                             SSAV.geno_ASE_exp$log.avgExp.AdultLarva <= q[3],]
SSAV.geno_ASE_exp_HI <- SSAV.geno_ASE_exp[SSAV.geno_ASE_exp$log.avgExp.AdultLarva > q[3],]
dim(SSAV.geno_ASE_exp_MED)
#########



# Density plot functions
##########
# TwoBootDens:
# Essentially the same as TwoBoot (boot_permute.R), 
# only with more adjustments for the density function.
# args: boot_dat = dataset to bootstrap
#       x_col = the column containing the data to bootstrap
#       groupBy = logical column which separates group 1 and group 2
#       boot_n = number of bootstrap replicates. Default is 1,000
# returns dataframe with bootstrapped densities for each point in x and 95% conf. intervals
TwoBootDens <- function(boot_dat, x_col, groupBy, 
                        boot_n = 1000){
  densities.within <- as_tibble(data_frame(bs = 1:boot_n) %>% 
                                 dplyr:: group_by(bs) %>% # make 1000 bootstrap replicates
                                  # sample randomly for each bootstrap
                                  dplyr::mutate(data = list(boot_dat %>% 
                                                       dplyr::group_by(.[[groupBy]]) %>% 
                                                       dplyr::sample_frac(size = 1, replace = T)))) %>% 
    unnest(c(bs, data)) %>% # create separate rows for each bootstrap replicate in the list
    dplyr::group_by(bs, .[[groupBy]]) %>% # group the data by bootstrap replicate and sig/non-sig
    # for each bootstrap replicate, 
    # calculate the density using 512 points (default) on the x-axis of each bootstrap dataset
    dplyr::do(tidy(density(.[[x_col]], 
                    from = min(boot_dat[[x_col]]), 
                    to = max(boot_dat[[x_col]]), 
                    n = 512))) %>%
    dplyr::rename(group := 2)
  
  # calculate the quantiles for each point from the bootstrapped estimates
  densities.qtiles <- densities.within %>% 
    dplyr::rename(!!x_col := x, dens = y) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(group, .[[x_col]]) %>% 
    dplyr::summarise(q05 = quantile(dens, 0.025),
              q50 = quantile(dens, 0.5),
              q95 = quantile(dens, 0.975)) %>% 
    dplyr::rename(!!x_col := 2)
  
  return(densities.qtiles) # return density extimates and quantiles
}

# DiffDens: 
# calculates the difference between bootstrapped densities of two groups
# args: bs_dat = dataframe from TwoBootDens 
#       x_col = column containing the bootstrapped data
# returns dataframe with calculated differences of bootstrapped densities 
#   and conf intervals (q05 - q95) between two groups 
DiffDens <- function(bs_dat, x_col){
  Sig <- bs_dat[bs_dat$group == TRUE, ]  # define group 1
  NS <- bs_dat[bs_dat$group == FALSE, ]  # define group 2
  # (if you want to calculate the opposite way (i.e. y-x, just switch the order you merge))
  diff <- merge(Sig, NS, by = {{x_col}}) # merge dataframe to one row 
  
  dens_diff <- data.frame(diff[[x_col]]) # make dataframe for the diff
  dens_diff <- dens_diff %>%
    dplyr::mutate(q05 = diff$q05.x - diff$q95.y, # calculate q05
           q50 = diff$q50.x - diff$q50.y, # calculate q05
           q95 = diff$q95.x - diff$q05.y) # calculate q95
  
  dens_diff$col <- "color" # pseudo column for plotting purposes (see plot_densities function)
  colnames(dens_diff)[1] <- print(x_col) # change column name to the right variable
  
  return(dens_diff) # return dataframe of diff
}


# density plotting function for downstream analyses
# args: dat = dataframe with log2FC SBGE data
#       bs_dat = dataframe from TwoBootDens
#       diff_dat = dataframe from DiffDens
#       x_col = bootstrapped data?
#       groupBy = logical column which separates group 1 and group 2
#       x_lab = label for x axis
plotDens <- function(dat, bs_dat, diff_dat, x_col, groupBy, x_lab){
  dens_plot <- ggplot(bs_dat, aes_string(x=x_col, y=bs_dat$q50)) + 
    # shade SE around density line
    geom_ribbon(data=bs_dat, 
                aes(ymin=q05, ymax=q95, color = group),
                size = 0.25, alpha = 0.7, fill = "grey") +
    # density lines for each group
    geom_density(data = dat, aes_string(x=x_col, y="..density..", color=groupBy), size = 1) +
    # shade SE around difference in densities line
    geom_ribbon(data=diff_dat, 
                aes_string(x=x_col, y="q50", ymin="q05", ymax="q95", colour = "col"), 
                size = 0.1, alpha = 0.7, fill = "grey") +
    # difference in densities line
    geom_line(data = diff_dat, aes_string(x=x_col, y = "q50", colour = "col"), size = 1.5) +
    # adjust colout
    scale_color_manual(name=NULL, labels=c("diff","background","candidates"),
                       values=c("black","orange2", "forestgreen")) +
    # adjust x and y axes labels
    ylab("Density") + xlab(print(x_lab)) +
    # horizontal line at 0
    geom_hline(yintercept = 0) +
    # theme settings below:
    theme_classic() +
    theme(legend.title = element_blank(),
          legend.position = c(0.85,0.85),
          legend.text = element_text(size = 20, color = "black"),
          axis.text.x = element_text(size=20, margin = margin(5,0,0,0), color = "black"),
          axis.text.y = element_text(size=20, margin = margin(0,5,0,0), color = "black"),
          axis.title.x = element_text(size=40, margin = margin(10,0,0,0), color = "black"),
          axis.title.y = element_text(size=40, margin = margin(0,10,0,0), color = "black"),
          plot.title = element_text(size=40, margin = margin(0,0,0,0), color = "black"),
          plot.margin = margin(6,6,6,6))
  return(dens_plot)
}

##########


# Density plots
##########
# Candidates vs non-candidate genes in SSAV males
boot_ALL <- TwoBootDens(SSAV.geno_ASE[!is.na(SSAV.geno_ASE$exp_SBGE_ase) & 
                                       !is.na(SSAV.geno_ASE$Sig),], 
                        "exp_SBGE_ase", "Sig")
diff_ALL <- DiffDens(boot_ALL, x_col = "exp_SBGE_ase")
ALL <- plotDens(dat = SSAV.geno_ASE[!is.na(SSAV.geno_ASE$exp_SBGE_ase) &
                                          !is.na(SSAV.geno_ASE$Sig),], 
                     bs_dat = boot_ALL, 
                     diff_dat = diff_ALL, 
                     x_col = "exp_SBGE_ase",
                     groupBy = "Sig", 
                     x_lab = "log2FC SBGE (ASE)") +
  scale_x_continuous(limits = c(-20, 20), breaks = seq(-25,25, by = 5)) +
  scale_y_continuous(limits = c(-0.175, 0.25))

boot_male_ALL <- TwoBootDens(A.m.geno_ASE[!is.na(A.m.geno_ASE$exp_SBGE_ase) & 
                                            !is.na(A.m.geno_ASE$Sig),], 
                             "exp_SBGE_ase", "Sig")
diff_male_ALL <- DiffDens(boot_male_ALL, x_col = "exp_SBGE_ase")
male_All <- plotDens(dat = A.m.geno_ASE[!is.na(A.m.geno_ASE$exp_SBGE_ase) &
                                    !is.na(A.m.geno_ASE$Sig),], 
                     bs_dat = boot_male_ALL, 
                     diff_dat = diff_male_ALL, 
                     x_col = "exp_SBGE_ase",
                     groupBy = "Sig", 
                     x_lab = "log2FC SBGE (ASE)") +
  scale_x_continuous(limits = c(-20, 20), breaks = seq(-25,25, by = 5)) +
  scale_y_continuous(limits = c(-0.175, 0.25))


# Candidates vs non-candidates in SSAV females
boot_fem_ALL <- TwoBootDens(A.f.geno_ASE[!is.na(A.f.geno_ASE$exp_SBGE_ase) & 
                                         !is.na(A.f.geno_ASE$Sig),], 
                          "exp_SBGE_ase", "Sig")
diff_fem_ALL <- DiffDens(boot_fem_ALL, "exp_SBGE_ase")
fem_All <- plotDens(A.f.geno_ASE[!is.na(A.f.geno_ASE$exp_SBGE_ase) & 
                                   !is.na(A.f.geno_ASE$Sig),], 
                    boot_fem_ALL, 
                    diff_fem_ALL, 
                    "exp_SBGE_ase", 
                    groupBy = "Sig", 
                    x_lab = "log2FC SBGE (ASE)") +
  scale_x_continuous(limits = c(-20, 20), breaks = seq(-25,25, by = 5)) +
  scale_y_continuous(limits = c(-0.175, 0.25))

##########


# Dot & whisker plot for proportion of candidates vs non-candidates per SBGE bin
##########

propSBGE <- function(dat, SBGE_cat){
  total_All <- dim(dat)[1]
  total_Sig <- dim(dat[dat$Sig,])[1]
  total_NS <- total_All - total_Sig
  
  total_SBGE <- dat %>% 
    dplyr::group_by(.[[SBGE_cat]]) %>%
    dplyr::summarise(SBGE_count = n()) %>%
    dplyr::rename({{SBGE_cat}} := 1)
  
  fract <- dat %>% 
    dplyr::group_by(.[[SBGE_cat]], Sig) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::mutate(frac_Sig = ifelse(Sig, count/total_Sig, count/total_NS)) %>%
    dplyr::rename({{SBGE_cat}} := 1) 
  
  fract <- merge(fract, total_SBGE, by = SBGE_cat)
  prop_SIG <- sum(fract[fract$Sig,]$count)/sum(fract[fract$Sig,]$SBGE_count)
  prop_NS <- sum(fract[!fract$Sig,]$count)/sum(fract[!fract$Sig,]$SBGE_count)
  
  fract <- fract %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(frac_by_SBGE_bin = count/SBGE_count,
           lower = list(binom.test(x=count, n=SBGE_count, p=ifelse(Sig, prop_SIG, prop_NS))), 
           upper = lower$conf.int[2], 
           lower = lower$conf.int[1])
  
  return(fract)
}


plotSBGEprop <- function(dat, SBGE_cat, xlab){
  plot_dat <- propSBGE(dat, SBGE_cat)
  
  prop_all <- data.frame(dim(dat[dat$Sig,])[1]/dim(dat)[1])
  prop_all$lower <- binom.test(dim(dat[dat$Sig,])[1], dim(dat)[1])$conf.int[1]
  prop_all$upper <- binom.test(dim(dat[dat$Sig,])[1], dim(dat)[1])$conf.int[2]
  
  plot_vec <- ggplot(plot_dat[plot_dat$Sig,], aes_string(x = SBGE_cat, y = "frac_by_SBGE_bin")) + 
    geom_hline(yintercept = prop_all[,1],
               linetype = "dashed", show.legend = TRUE, size = 0.75) +
    annotate('ribbon', x = c(-Inf, Inf), ymin = prop_all$lower, ymax = prop_all$upper, 
             alpha = 0.20, fill = 'grey30') +
    geom_errorbar(ymin = plot_dat[plot_dat$Sig,]$lower, ymax = plot_dat[plot_dat$Sig,]$upper,
                  width = 0.5, size = 1, color = "#009E73") +
    geom_point(fill = "#009E73", color = "#009E73", size = 7) + 
    labs(x = "Sex-biased Gene Expression", # "omegaA_MK" = expression(italic("\u03c9A")[MK]); "alpha_MK" = expression(italic("\u03b1")[MK])
         y = "Proportion of DE Genes") +
    scale_x_discrete(labels = c("H.FB", "FB", "UB", "MB", "H.MB")) +
    scale_y_continuous(limits = c(0, 0.20)) +
    geom_text(data = plot_dat[plot_dat$Sig,], aes_string(x = SBGE_cat, y = "upper" , label = "count" ), 
              size = 7.5, vjust = -0.7) +
    # geom_vline(xintercept = c(seq(1.5, length(levels(dat[[SBGE_cat]])) )), color = "grey") +
    theme_classic() +
    theme(plot.title.position = c("panel"),
          legend.title = element_blank(),
          legend.position = c("bottom"),
          #legend.justification = c("right", "bottom"),
          #legend.box.just = "left",
          #legend.box.background = element_rect(),
          legend.box.background = element_rect(),
          #legend.box.margin = margin(4, 6, 6, 6),
          legend.text = element_text(size = 20, color = "black"),
          axis.text.x = element_text(size=20, margin = margin(5,0,5,0), color = "black"),
          axis.text.y = element_text(size=20, margin = margin(0,5,0,0), color = "black"),
          axis.title.x = element_text(size=30, margin = margin(10,0,10,0), color = "black"),
          axis.title.y = element_text(size=30, margin = margin(0,10,0,0), color = "black"),
          plot.title = element_text(size=40, margin = margin(0,0,0,0), color = "black"),
          plot.margin = margin(6,6,6,6))
  return(plot_vec)
}



propSBGE(SSAV.geno_ASE_exp_HI, "SBGE_comp")
propSBGE(A.f.geno_ASE, "SBGE_comp")
propSBGE(A.m.geno_ASE, "SBGE_comp")

bin_All <- plotSBGEprop(SSAV.geno_ASE, "SBGE_comp", "SBGE (ASE)")
bin_All_LOW <- plotSBGEprop(SSAV.geno_ASE_exp_LOW, "SBGE_comp", "SBGE (ASE)")
bin_All_MED <- plotSBGEprop(SSAV.geno_ASE_exp_MED, "SBGE_comp", "SBGE (ASE)") +   
  scale_x_discrete(labels = c("FB", "UB", "MB", "H.MB"))
bin_All_HI <- plotSBGEprop(SSAV.geno_ASE_exp_HI, "SBGE_comp", "SBGE (ASE)") + 
  coord_cartesian(ylim = c(0, 0.4))


exp_levels_SBGE <- ggarrange(bin_All + theme(axis.title.x = element_blank(),
                                            axis.title.y = element_blank()),
                            ggarrange(bin_All_LOW + theme(axis.title.x = element_blank(),
                                                          axis.title.y = element_blank()), 
                                      bin_All_MED + theme(axis.title.x = element_blank(),
                                                          axis.title.y = element_blank()), 
                                      bin_All_HI + theme(axis.title.x = element_blank(),
                                                         axis.title.y = element_blank()), 
                                      ncol = 3),
                            nrow = 2, heights = c(2, 1))


pdf(file = "~/Desktop/UofT/SSAV_RNA/Plots/test_SBGE_exp.pdf",   # The directory you want to save the file in
    width = 14, # 12, 24, 20 The width of the plot in inches
    height = 14) # 10, 20, 13 The height of the plot in inches

annotate_figure(exp_levels_SBGE + theme(plot.margin = margin(10,10,10,0)), 
                left = text_grob("freq",rot = 90, size = 40),
                bottom = text_grob("Sex-biased gene expression", size = 40))
dev.off()


bin_A.f <- plotSBGEprop(A.f.geno_ASE, "SBGE_comp", "SBGE (ASE)")
bin_A.m <- plotSBGEprop(A.m.geno_ASE, "SBGE_comp", "SBGE (ASE)") + coord_cartesian(ylim = c(0, 0.1))

##########




pdf(file = "~/Desktop/UofT/SSAV_RNA/Plots/final_2/Fig2_main_filtered.pdf",   # The directory you want to save the file in
    width = 14, # 14 24 The width of the plot in inches
    height = 10) # 10 20 The height of the plot in inches
# ggarrange(bin_A.f, NA, bin_A.m, NA, NA, NA, fem_All, NA, male_All,
#           labels = c("A)", NA, "B)", NA, NA, NA, "C)", NA, "D)"),
#           widths = c(1, 0.05, 1),
#           heights = c(1, 0.05, 1),
#           ncol = 3, nrow = 3,
#           font.label = list(size = 30), hjust = -0.01)

bin_All + coord_cartesian(ylim = c(0,0.15))
# 
# ggarrange(bin_A.m, NA, bin_A.f,
#           labels = c("A)", NA, "B)"),
#           widths = c(1, 0.05, 1),
#           ncol = 3,
#           font.label = list(size = 30), hjust = -0.01)
dev.off()

##########


# calculate average log2FC male female in our data and all 
# figure 6 "Experimental vs. Control Difference: Red Males"
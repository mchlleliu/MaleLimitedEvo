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


# packages required
#########
library(tidyr)
library(plyr)
library(dplyr)
library(broom)
library(ggplot2)
########



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
                                  group_by(bs) %>% # make 1000 bootstrap replicates
                                  # sample randomly for each bootstrap
                                  mutate(data = list(boot_dat %>% 
                                                       group_by(.[[groupBy]]) %>% 
                                                       sample_frac(size = 1, replace = T)))) %>% 
    unnest(c(bs, data)) %>% # create separate rows for each bootstrap replicate in the list
    group_by(bs, .[[groupBy]]) %>% # group the data by bootstrap replicate and sig/non-sig
    # for each bootstrap replicate, 
    # calculate the density using 512 points (default) on the x-axis of each bootstrap dataset
    do(tidy(density(.[[x_col]], 
                    from = min(boot_dat[[x_col]]), 
                    to = max(boot_dat[[x_col]]), 
                    n = 512))) %>%
    rename(group := 2)
  
  # calculate the quantiles for each point from the bootstrapped estimates
  densities.qtiles <- densities.within %>% 
    rename(!!x_col := x, dens = y) %>%
    ungroup() %>%
    group_by(group, .[[x_col]]) %>% 
    summarise(q05 = quantile(dens, 0.025),
              q50 = quantile(dens, 0.5),
              q95 = quantile(dens, 0.975)) %>% 
    rename(!!x_col := 2)
  
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
    mutate(q05 = diff$q05.x - diff$q95.y, # calculate q05
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
  
  total_SBGE <- dat %>% group_by(.[[SBGE_cat]]) %>%
    summarise(SBGE_count = n()) %>%
    rename({{SBGE_cat}} := 1)
  
  fract <- dat %>% group_by(.[[SBGE_cat]], Sig) %>%
    summarise(count = n()) %>%
    mutate(frac_Sig = ifelse(Sig, count/total_Sig, count/total_NS)) %>%
    rename({{SBGE_cat}} := 1) 
  
  fract <- merge(fract, total_SBGE, by = SBGE_cat)
  prop_SIG <- sum(fract[fract$Sig,]$count)/sum(fract[fract$Sig,]$SBGE_count)
  prop_NS <- sum(fract[!fract$Sig,]$count)/sum(fract[!fract$Sig,]$SBGE_count)
  
  fract <- fract %>% 
    rowwise() %>% 
    mutate(frac_by_SBGE_bin = count/SBGE_count,
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
                  width = 0.5, size = 0.75) +
    geom_point(fill = "forestgreen", color = "forestgreen", size = 7) + 
    labs(x = "Sex-biased Gene Expression", # "omegaA_MK" = expression(italic("\u03c9A")[MK]); "alpha_MK" = expression(italic("\u03b1")[MK])
         y = "Proportion of Candidate Genes") +
    scale_x_discrete(labels = c("Highly FB", "Female-Biased", "Unbiased", "Male-Biased", "Highly MB")) +
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



propSBGE(SSAV.geno_ASE, "SBGE_comp")
propSBGE(A.f.geno_ASE, "SBGE_comp")
propSBGE(A.m.geno_ASE, "SBGE_comp")

bin_All <- plotSBGEprop(SSAV.geno_ASE, "SBGE_comp", "SBGE (ASE)")
bin_A.f <- plotSBGEprop(A.f.geno_ASE, "SBGE_comp", "SBGE (ASE)")
bin_A.m <- plotSBGEprop(A.m.geno_ASE, "SBGE_comp", "SBGE (ASE)") + coord_cartesian(ylim = c(0, 0.05))

##########




pdf(file = "~/Desktop/UofT/SSAV_RNA/Plots/finals/Fig2_main.pdf",   # The directory you want to save the file in
    width = 14, # 14 24 The width of the plot in inches
    height = 10) # 10 20 The height of the plot in inches
# ggarrange(bin_A.f, NA, bin_A.m, NA, NA, NA, fem_All, NA, male_All,
#           labels = c("A)", NA, "B)", NA, NA, NA, "C)", NA, "D)"),
#           widths = c(1, 0.05, 1),
#           heights = c(1, 0.05, 1),
#           ncol = 3, nrow = 3,
#           font.label = list(size = 30), hjust = -0.01)

bin_All_DsRed

# ggarrange(bin_A.m_DsRed, NA, bin_A.f_DsRed,
#           labels = c("A)", NA, "B)"),
#           widths = c(1, 0.05, 1),
#           ncol = 3,
#           font.label = list(size = 30), hjust = -0.01)
dev.off()

##########
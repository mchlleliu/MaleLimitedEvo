###################################
#
#                             Grieshop et al. 2023
#                              Author: Michelle Liu
#             DsRed experimental evolution - transcriptomics analysis
#     Correlation of Red/NR changes in Experimental samples vs Control samples
#                                   Figure S7 
#
# 
###################################

rm(list=ls())
setwd("~/Desktop/UofT/SSAV_RNA/")

# load packages
#########
library(broom)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggblend)
library(ggpubr)
#########


# the two functions below are used in the corr plot function below

# quad_count: counts the number of point in each of the 8 quadrants
# variables: dat = data.frame object containing the values to be plotted
#            x = values for group 1, plotted on the x-axis
#            y = values for group 2, plotted on the y-axis
#            lim = for plotting code. specificy the x-lim and y-lim of the plot
# return: a data frame with number of points, 
#         proportion relative to total points, 
#         and plotting coordinates
quad_count <- function(dat, x, y, lim){
  count <- dat[!is.na(dat[[x]]) & !is.na(dat[[y]]),] %>%
    # Count how many with each combination of X and Y being positive
    dplyr::count(right = .[[x]] > 0, # on the right side of plot?
                 top = .[[y]] > 0) %>%  # on the top side of plot?
    
    dplyr::mutate(perc = n/sum(n)) %>% # calculate percentage of points relative to total number of points
    
    # specificy coordinates for texts on plot
    dplyr::mutate(!!x := lim/2*(2.55*(right - 0.5)), 
           !!y := lim/2*(2.55*(top - 0.5)))
  
  return(count)
}


# colour_quadrant: add a column to specify colour for each point based on where they are located,
#                 for points in the concordant quadrants (I and III)
# variables: dat = data.frame object containing the values to be plotted
#            x = values for group 1, plotted on the x-axis
#            y = values for group 2, plotted on the y-axis
#            colNonCon = colour for quadrants II and IV
# return: a data frame with color specified for each data point
colour_quadrant <-  function(dat, x, y, colCon, colNonCon){
  col <- dat %>%
    # logical columns to define where the point is locates
    dplyr::mutate(right = .[[x]] > 0, # on the right part of plot?
           top = .[[y]] > 0) %>% # on the top part of plot?
           
    # add the colour
    dplyr::mutate(quadrant = ifelse(right & top | (!right & !top), colCon, colNonCon))
  return(col)
}


# Bootstrapping correlation (modified from boot_permute.R)
##########
# Bootstrapped 95% confidence intervals for correlation between two groups
BootCorr <- function(boot_dat, x_col, y_col, boot_n = 1000){
  boot_tabs <- as_tibble(data_frame(bs = 1:boot_n) %>% # make boot_n bootstrap replicates
                           dplyr::group_by(bs) %>% # for each bootstrap,
                           # sample randomly x boot_n times from each group
                           dplyr::mutate(data = list(boot_dat %>% 
                                                     dplyr::sample_frac(size = 1, replace = T)))) %>% 
    unnest(c(bs, data)) %>% # create separate rows for each bootstrap replicate in the list
    dplyr::group_by(bs) %>% # group the data by bootstrap replicate
    # for each bootstrap replicate, calculate the correlation
    dplyr::do(tidy(cor.test(.[[x_col]], .[[y_col]])$estimate))
  
  print(head(boot_tabs))
  
  # summarise bootstrap replicates
  boot_SE <- boot_tabs %>%
    dplyr::ungroup() %>%
    dplyr::summarise(q05 = quantile(x, 0.025, na.rm = TRUE),
                     q50 = quantile(x, 0.5, na.rm = TRUE),
                     q95 = quantile(x, 0.975, na.rm = TRUE))
  return(boot_SE)
}


# plot_corr: add a column to specify colour for each point based on where they are located,
#                 for points in the concordant quadrants (I and III)
# variables: dat = data.frame object containing the values to be plotted
#            x = values for group 1, plotted on the x-axis
#            y = values for group 2, plotted on the y-axis
#            colx = colour for x > y
#            coly = colour for x < y
#            colNonCon = colour for quadrants II and IV
#            xlab = label for x-axis
#            ylab = label for y-axis
#            lim = x and y axes limit
#            title = of graph
plot_corr <- function(dat, x, y, colCon, colNonCon, xlab, ylab, lim, title = ""){
  # do correlation test
  pear_cor <- cor.test(dat[[x]], dat[[y]])$estimate
  print(pear_cor)
  
  # bootstrap the confidence interval for the correlation
  boot_CI <- BootCorr(boot_dat = dat, x_col = x, y_col = y, boot_n = 1000)
  print(boot_CI)

  # count the percentages
  quad_n <- quad_count(dat, x, y, lim)
  # manage the colour of points
  quad_col <- colour_quadrant(dat, x, y, colCon, colNonCon)
  
  # plot
  corr <- ggplot(dat, aes_string(x = x, y = y)) +
    geom_point(size = 3, shape = 16, alpha = 0.75, color = quad_col$quadrant) +
    geom_abline(intercept = 0, slope = 0,  size = 0.5, linetype="solid", color = "black") +
    geom_hline(yintercept = 0,  size = 0.5, linetype="solid", color = "black") +
    geom_vline(xintercept = 0,  size = 0.5, linetype="solid", color = "black") +
    # geom_abline(intercept = 0, slope = 1,  size = 0.5, linetype="dashed", color = "black") +
    # geom_abline(intercept = 0, slope = -1,  size = 0.5, linetype="dashed", color = "black") +
    coord_cartesian(xlim=c(-lim, lim), ylim = c(-lim,lim)) +
    # add correlation
    geom_label(aes(x = 1.5, y = -1.5,
                   label = c(paste("r:", round(pear_cor, digits = 3),
                                   "\n[", round(boot_CI$q05, digits = 3),
                                   ",", round(boot_CI$q95, digits = 3),
                                   "]", sep = " "))), data = quad_n, size = 10) +


    # add percentages
    # geom_text(aes(label = paste(round(perc*100,digits=1),"%",sep="")), data = quad_n, size = 10) +
    # coord_cartesian(xlim=c(-lim, lim), ylim = c(-lim,lim)) + # change to (0, lim) or (-lim, lim) depending on the figure
    # labs(x = print(xlab),
    #      y = print(ylab) ,
    #      title = print(title)) +
    # guides(color = guide_legend(override.aes = list(shape = c(NA, NA), # c(16, 16)
    #                                                 size = c(4, 4),
    #                                                 alpha = 1))) +

    labs(x = print(xlab),
         y = print(ylab) ,
         title = print(title)) +
    guides(color = guide_legend(override.aes = list(shape = c(NA, NA), # c(16, 16)
                                                    size = c(4, 4),
                                                    alpha = 1))) +

    # some theme settings...
    theme_classic() +
    theme(legend.title = element_blank(),
          legend.position = c("None"),
          #legend.justification = c("right", "bottom"),
          #legend.box.just = "left",
          legend.box.background = element_rect(),
          #legend.box.margin = margin(4, 6, 6, 6),
          legend.text = element_text(size = 20, color = "black"),
          plot.tag = element_text(size = 20, color = "black"),
          axis.text.x = element_text(size=20, margin = margin(5,0,0,0), color = "black"),
          axis.text.y = element_text(size=20, margin = margin(0,5,0,0), color = "black"),
          axis.title.x = element_text(size=30, margin = margin(10,0,0,0), color = "black"),
          axis.title.y = element_text(size=30, margin = margin(0,10,0,0), color = "black"),
          plot.title = element_text(size=40, margin = margin(0,0,0,0), color = "black"),
          plot.margin = margin(6,6,6,6),
          panel.border = element_rect(colour = "black", fill=NA, size=1)
    )

  return(corr) # return plot
}



# set up data frames for correlation plot
########
# load datasets if not loaded yet
A.f.geno <- read.delim("Results/A.f.geno_DE.candidates.tsv")
A.m.geno <- read.delim("Results/A.m.geno_DE.candidates.tsv")
C.m.geno <- read.delim("Results/C.m.geno_DE.candidates.tsv")

tmp.males <- A.m.geno
colnames(tmp.males) <- c("m.exp_geno", "m.se_geno", "m.padj", "FlyBaseID", "m.TopSig", "m.Sig")
head(tmp.males)

tmp.females <- A.f.geno
colnames(tmp.females) <- c("f.exp_geno", "f.se_geno", "f.padj", "FlyBaseID", "f.Sig")
head(tmp.females)

tmp.C.males <- C.m.geno
colnames(tmp.C.males) <- c("C.m.exp_geno", "C.m.se_geno", "C.m.padj", "FlyBaseID", "C.m.Sig")

########


# plotting only female candidate SA genes
########
# Exp females vs Control males
corr.plot <- merge(tmp.females, tmp.C.males, by = "FlyBaseID", all = T)
CmAf_fem_cand <- plot_corr(corr.plot[corr.plot$FlyBaseID %in% A.f.geno[A.f.geno$Sig,]$FlyBaseID, ],
                   x = "f.exp_geno", y = "C.m.exp_geno", 
                   "black", "grey65",
                   "Red/NR in Exp females", "Red/NR in Control males", 
                   2.5, "") +
  labs(x = expression(Log["2"]*"FC ("~italic(Red)~"/"~italic(NR)~") in Exp. Females"),
       y = expression(Log["2"]*"FC ("~italic(Red)~"/"~italic(NR)~") in Ctrl. Males"))

# Exp females vs Exp males
corr.plot <- merge(tmp.females, tmp.males, by = "FlyBaseID", all = T)
AmAf_fem_cand <- plot_corr(corr.plot[corr.plot$FlyBaseID %in% A.f.geno[A.f.geno$Sig,]$FlyBaseID, ], 
                   "f.exp_geno", "m.exp_geno", 
                   "black", "grey65",
                   "Red/NR in Exp females", "Red/NR in Exp males", 
                   2.5, "") +
  labs(x = expression(Log["2"]*"FC ("~italic(Red)~"/"~italic(NR)~") in Exp. Females"),
       y = expression(Log["2"]*"FC ("~italic(Red)~"/"~italic(NR)~") in Exp. Males"))

########


# plotting only male candidate SA genes
########

# SSAV males vs Control males
corr.plot <- merge(tmp.males, tmp.C.males, by = "FlyBaseID", all = T)
CmAm_male_cand <- plot_corr(corr.plot[corr.plot$FlyBaseID %in% A.m.geno[A.m.geno$Sig,]$FlyBaseID, ], 
                            "m.exp_geno", "C.m.exp_geno",  
                            "black", "grey65",
                            "Red/NR in Exp males", "Red/NR in Control males", 
                            2.5, "") +
  labs(x = expression(Log["2"]*"FC ("~italic(Red)~"/"~italic(NR)~") in Exp. Males"),
       y = expression(Log["2"]*"FC ("~italic(Red)~"/"~italic(NR)~") in Ctrl. Males"))

# SSAV males vs SSAV females
corr.plot <- merge(tmp.males, tmp.females, by = "FlyBaseID", all = T)
AfAm_male_cand <- plot_corr(corr.plot[corr.plot$FlyBaseID %in% A.m.geno[A.m.geno$Sig,]$FlyBaseID, ], 
                            "m.exp_geno", "f.exp_geno", 
                            "black", "grey65",
                            "Red/NR in Exp. Males", "Red/NR in Exp. Females", 
                            2.5, "") +
  labs(x = expression(Log["2"]*"FC ("~italic(Red)~"/"~italic(NR)~") in Exp. Males"),
       y = expression(Log["2"]*"FC ("~italic(Red)~"/"~italic(NR)~") in Exp. Females"))

########



# Save the plots
png(file = "~/Desktop/UofT/SSAV_RNA/Plots/revision/Fig_S7.png",   # The directory you want to save the file in
    width = 17, # The width of the plot in inches
    height = 17, # The height of the plot in inches
    units = "in", res = 200) 
ggarrange(AfAm_male_cand, NA, AmAf_fem_cand, NA,
          NA, NA, NA, NA,
          CmAm_male_cand, NA, CmAf_fem_cand, NA,
          labels = c("  A)", NA, "  B)",
                     NA,NA,NA,
                     "  C)", NA, "  D)"),
          ncol = 4, nrow = 3,
          widths = c(1, 0.15, 1, 0.05),
          heights = c(1, 0.1, 1),
          font.label = list(size = 30)) 
dev.off()


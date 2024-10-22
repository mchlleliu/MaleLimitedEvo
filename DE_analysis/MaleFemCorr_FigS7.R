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
  count <- dat %>%
    # Count how many with each combination of X and Y being positive
    dplyr::count(right = .[[x]] > 0, # on the right side of plot?
          top = .[[y]] > 0, # on the top side of plot?
          # for each quadrant, divide into two: 
          # (this is a bit weird, but needed so the numbers can be plotted at the right coordinates)
          UP = !top & abs(.[[x]]) > abs(.[[y]]) | # quadrant III up-left or quadrant IV up-right
            top & abs(.[[x]]) < abs(.[[y]])) %>%  # quadrant I up-right or quadrant II up-left
    
    dplyr::mutate(perc = n/sum(n)) %>% # calculate percentage of points relative to total number of points
    
    # this is another strange one for setting up the coordinates
    dplyr::mutate(conc = right & top | (!right & !top), # quadrant I and quadrant III (the concordant changes)
           dir_UP = conc & UP, # concordant changes where x > y
           dir_DOWN = conc & !UP) %>% # concordant changes where x < y
    
    # TRUE = 1, FALSE = 0
    # specificy coordinates for texts on plot
    dplyr::mutate(!!x := lim/2*(2*(right - 0.5)+(UP - 0.5)+((conc-0.001)*0.5)-(dir_UP*1.5)+(dir_DOWN*0.5)), 
           !!y := lim/2*(2*(top - 0.5)+(UP - 0.5)))
  
  return(count)
}


# colour_quadrant: add a column to specify colour for each point based on where they are located,
#                 for points in the concordant quadrants (I and III)
# variables: dat = data.frame object containing the values to be plotted
#            x = values for group 1, plotted on the x-axis
#            y = values for group 2, plotted on the y-axis
#            colx = colour for x > y
#            coly = colour for x < y
#            colNonCon = colour for quadrants II and IV
# return: a data frame with color specified for each data point
colour_quadrant <-  function(dat, x, y, colx, coly, colNonCon){
  col <- dat %>%
    # logical columns to define where the point is locates
    dplyr::mutate(right = .[[x]] > 0, # on the right part of plot?
           top = .[[y]] > 0, # on the top part of plot?
           # for the concordant quadrants (I & III)... 
           DOWN = !top & abs(.[[x]]) > abs(.[[y]]) | # quadrant III where x > y
             top & abs(.[[x]]) > abs(.[[y]])) %>% # quadrant I where x > y
    # add the colour
    dplyr::mutate(quadrant = ifelse(right & top | (!right & !top), 
                             ifelse(DOWN, colx, coly), 
                             colNonCon))
  return(col)
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
plot_corr <- function(dat, x, y, colx, coly, colNonCon, xlab, ylab, lim, title){
  # do correlation test
  pear_cor <- cor.test(dat[[x]], dat[[y]],method = "pearson")
  # count the percentages
  quad_n <- quad_count(dat, x, y, lim)
  # manage the colour of points
  quad_col <- colour_quadrant(dat, x, y, colx, coly, colNonCon)
  
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
                   label = c(paste("r:", round(pear_cor$estimate, digits = 3), 
                                   "\n[", round(pear_cor$conf.int[1], digits = 3),
                                   ",", round(pear_cor$conf.int[2], digits = 3), 
                                   "]", sep = " "))), data = quad_n, size = 10) +
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
          #legend.box.background = element_rect(),
          legend.box.background = element_rect(),
          #legend.box.margin = margin(4, 6, 6, 6),
          legend.text = element_text(size = 20, color = "black"),
          plot.tag = element_text(size = 20, color = "black"),
          axis.text.x = element_text(size=20, margin = margin(5,0,0,0), color = "black"),
          axis.text.y = element_text(size=20, margin = margin(0,5,0,0), color = "black"), 
          axis.title.x = element_text(size=30, margin = margin(10,0,0,0), color = "black"),
          axis.title.y = element_text(size=30, margin = margin(0,10,0,0), color = "black"),
          plot.title = element_text(size=40, margin = margin(0,0,0,0), color = "black"),
          plot.margin = margin(6,6,6,6)
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

# SSAV females vs Control males
corr.plot <- merge(tmp.females, tmp.C.males, by = "FlyBaseID", all = T)
CmAf_fem_cand <- plot_corr(corr.plot[corr.plot$FlyBaseID %in% A.f.geno[A.f.geno$Sig,]$FlyBaseID, ],
                   x = "f.exp_geno", y = "C.m.exp_geno", 
                   "black", "black", "grey65",
                   "Red/NR in SSAV females", "Red/NR in Control males", 
                   2.5, "") +
  labs(x = expression(Log["2"]*"FC ("~italic(Red)~"/"~italic(NR)~") in Exp. Females"),
       y = expression(Log["2"]*"FC ("~italic(Red)~"/"~italic(NR)~") in Ctrl. Males"))

# SSAV females vs SSAV males
corr.plot <- merge(tmp.females, tmp.males, by = "FlyBaseID", all = T)
AmAf_fem_cand <- plot_corr(corr.plot[corr.plot$FlyBaseID %in% A.f.geno[A.f.geno$Sig,]$FlyBaseID, ], 
                   "f.exp_geno", "m.exp_geno", 
                   "black", "black", "grey65",
                   "Red/NR in SSAV females", "Red/NR in SSAV males", 
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
                            "black", "black", "grey65",
                            "Red/NR in SSAV males", "Red/NR in Control males", 
                            2.5, "") +
  labs(x = expression(Log["2"]*"FC ("~italic(Red)~"/"~italic(NR)~") in Exp. Males"),
       y = expression(Log["2"]*"FC ("~italic(Red)~"/"~italic(NR)~") in Ctrl. Males"))

# SSAV males vs SSAV females
corr.plot <- merge(tmp.males, tmp.females, by = "FlyBaseID", all = T)
AfAm_male_cand <- plot_corr(corr.plot[corr.plot$FlyBaseID %in% A.m.geno[A.m.geno$Sig,]$FlyBaseID, ], 
                            "m.exp_geno", "f.exp_geno", 
                            "black", "black", "grey65",
                            "Red/NR in Exp. Males", "Red/NR in Exp. Females", 
                            2.5, "") +
  labs(x = expression(Log["2"]*"FC ("~italic(Red)~"/"~italic(NR)~") in Exp. Males"),
       y = expression(Log["2"]*"FC ("~italic(Red)~"/"~italic(NR)~") in Exp. Females"))

########


# Save the plots
png(file = "~/Desktop/UofT/SSAV_RNA/Plots/final_2/png_version/Fig_S7.png",   # The directory you want to save the file in
    width = 17, # The width of the plot in inches
    height = 17, # The height of the plot in inches
    units = "in", res = 200) 
ggarrange(AfAm_male_cand, NA, AmAf_fem_cand,
          NA, NA, NA,
          CmAm_male_cand, NA, CmAf_fem_cand,
          labels = c("A)", NA, "B)",
                     NA,NA,NA,
                     "C)", NA, "D)"),
          ncol = 3, nrow = 3,
          widths = c(1, 0.05, 1),
          heights = c(1, 0.05, 1),
          font.label = list(size = 30))
dev.off()


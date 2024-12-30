###################################
#
#                             Grieshop et al. 2023
#                             Author: Michelle Liu  
#             DsRed experimental evolution - transcriptomics analysis
#          Correlation of Evolutionary changes in Red vs NonRed samples 
#                             relative to Controls
#                               Figure 6A & S8
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
library(ggrepel)
library(ggblend)
library(grid)
library(gridtext)
#########


# Get Mishra et al.2022 data for SBGE estimates
source("Mishra_et.al_SBGE.R")

# Load & set up plotting datasets
########
# Load data sets if not yet in env
tmp.Red <- read.delim("Results/Red.m.trt_raw.tsv")
colnames(tmp.Red) <- c("Red.exp_trt", "Red.se_trt", "Red.padj", "FlyBaseID")

tmp.NR <- read.delim("Results/NR.m.trt_raw.tsv")
colnames(tmp.NR) <- c("NR.exp_trt", "NR.se_trt", "NR.padj", "FlyBaseID")

corr.plot <- merge(tmp.Red, tmp.NR, by = "FlyBaseID")
# combine Red and NonRed data to create a plotting dataset


corr.plot <- merge(corr.plot, Mishra, by = "FlyBaseID")
# add in SBGE values from Mishra et al. 2022 data
########



# the two functions below are used in the corr plot function below
# quad_count: counts the number of point in each of the 8 quadrants
# variables: dat = data.frame object containing the values to be plotted
#            x = values for group 1, plotted on the x-axis
#            y = values for group 2, plotted on the y-axis
#            lim = for plotting code. specificy the x-lim and y-lim of the plot
# return: a data frame with number of points, 
#         proportion relative to total points, 
#         and plotting coordinates
quad_count <- function(dat, x, y, lim = 5){
  
  # dat <- dat[dat[[x]] != dat[[y]],]
  # 
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
  
  print(count)
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
    mutate(right = .[[x]] > 0, # on the right part of plot?
           top = .[[y]] > 0, # on the top part of plot?
           # for the concordant quadrants (I & III)... 
           DOWN = !top & abs(.[[x]]) > abs(.[[y]]) | # quadrant III where x > y
             top & abs(.[[x]]) > abs(.[[y]])) %>% # quadrant I where x > y
    # add the colour
    mutate(quadrant = ifelse(right & top | (!right & !top), 
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
  # count the percentages
  quad_n <- quad_count(dat, x, y, lim)
  # manage the colour of points
  quad_col <- colour_quadrant(dat, x, y, colx, coly, colNonCon)
  # plot
  corr <- ggplot(dat, aes_string(x = x, y = y)) +
    geom_point(size = 2, shape = 16, alpha = 0.5, color = quad_col$quadrant) +  
    
    # add lines to separate quadrants
    geom_hline(yintercept = 0,  size = 0.5, linetype="solid", color = "black") + # comment out for Fig 6A
    geom_vline(xintercept = 0,  size = 0.5, linetype="solid", color = "black") + # comment out for Fig 6A
    geom_abline(intercept = 0, slope = 1,  size = 0.5, linetype="dashed", color = "black") +
    geom_abline(intercept = 0, slope = -1,  size = 0.5, linetype="dashed", color = "black") + # comment out for Fig. 6A
    
    # add percentages
    geom_text(aes(label = paste(round(perc*100,digits=0),"%",sep="")), data = quad_n, size = 10) +
    coord_cartesian(xlim=c(-lim, lim), ylim = c(-lim,lim)) + # change to (0, lim) or (-lim, lim) depending on the figure
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
          legend.box.background = element_rect(),
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



# plotting Figure S13
######
# remember to edit the plotting function (plot_corr) accordingly
Figure_S13 <- plot_corr(corr.plot, 
                       "Red.exp_trt", "NR.exp_trt", 
                       "red3", "black", "darkgrey",
                       "Log2FC Exp/Ctrl Red males", 
                       "Log2FC Exp/Ctrl NonRed males", 
                       2.5, "") +
  labs(x = expression(atop(paste(log["2"]*"FC (Exp:Ctrl)"), paste("in ",italic("Red")," Males"))),
       y = expression(atop(paste(log["2"]*"FC (Exp:Ctrl)"), paste("in ",italic("NonRed")," Males")))) +
  theme(title = element_blank())


# comment in to save plot
png(file = "~/Desktop/UofT/SSAV_RNA/Plots/final_2/png_version/Fig_S13.png",   # The directory you want to save the file in
    width = 10,
    height = 10,
    units = "in", res = 300)

Figure_S13

dev.off()

######



# plotting Fig. 6A
######
source("Fig6_effectFuns.R")
# remember to edit the plotting function (plot_corr) accordingly
# get only genes where Red and NonRed samples evolved in the same regulatory direction 
# (i.e., both Red and NR UPregulated relative to control or both DOWNregulated)
concordant <- corr.plot[!is.na(corr.plot$Red.exp_trt) & !is.na(corr.plot$NR.exp_trt) &
                          ((corr.plot$Red.exp_trt >= 0 & corr.plot$NR.exp_trt >= 0) |
                             (corr.plot$Red.exp_trt <= 0 & corr.plot$NR.exp_trt <= 0)),]
hist(concordant$Red.exp_trt - concordant$NR.exp_trt)

t.test(concordant$Red.exp_trt, concordant$NR.exp_trt, paired = T)
# correlation between evolutionary changes in Red and NonRed lines is lowly (but significantly) negative

# plot only the concordant changes.
Figure_6A <- plot_corr(concordant, 
                       "Red.exp_trt", "NR.exp_trt", 
                       "red3", "grey9", "darkgrey",
                       "SSAV/Control in Red males", 
                       "SSAV/Control in NonRed males", 
                       1.78, "") +
  labs(x = expression(atop(paste("Absolute "*log["2"]*"FC"), paste("in ",italic("Red")," Males"))),
       y = expression(atop(paste("Absolute "*log["2"]*"FC"), paste("in ",italic("NonRed")," Males")))) +
  theme(title = element_blank())  +
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2)) +
  scale_y_continuous(breaks = c(0, 0.5, 1, 1.5, 2)) +
  coord_cartesian(xlim = c(0, 1.75), ylim = c(0, 1.75))


# put together figure above and figure from the corresponding analysis on splicing changes
# (see: "SplicingPhi_Fig5.6B_TableS6.R" to generate Figure 6B)
Figure6 <- ggarrange(NA,NA,NA,NA,NA, 
                     NA, Figure_6A + ggtitle(expression(bold("Differential Expression"))) +
                       theme(axis.title.x = element_blank(), 
                             axis.title.y = element_blank(),
                             plot.title = element_text(hjust = 0.5, size = 30, vjust = 3)),
                     NA, Figure_6B + ggtitle(expression(bold("Differential Splicing"))) +
                       theme(axis.title.x = element_blank(), 
                             axis.title.y = element_blank(),
                             plot.title = element_text(hjust = 0.5, size = 30, vjust = 3)), NA,
                     ncol = 5,nrow = 3, widths = c(0.1, 1, 0.075, 1, 0.1), heights = c(0.1,1, 0.05),
                     labels = c(NA,NA,NA,NA,NA,
                                NA, "A)", NA, "B)", NA), hjust=0.75, vjust = -0.25,
                     font.label = list(size = 40))



# comment in to save figure 6
png(file = "~/Desktop/UofT/SSAV_RNA/Plots/final_2/png_version/Fig6_main.png",   # The directory you want to save the file in
    width = 20, # 20; The width of the plot in inches
    height = 10, # 30; The height of the plot in inches
    units = "in", res = 300)
    
annotate_figure(Figure6, left = richtext_grob("<span style='font-size:30pt; color:black'>Experimental vs. Control Difference:<br>*NonRed* Males</span>",
                                          rot = 90, vjust = 1),
                bottom = richtext_grob("<span style='font-size:30pt; color:black'>Experimental vs. Control Difference:<br>*Red* Males</span>",
                                       r = unit(50, "pt")))

dev.off()
######
t.test(abs(corr.plot$Red.exp_trt), abs(corr.plot$NR.exp_trt), paired = TRUE, alternative = "less")
t.test(abs(concordant$Red.exp_trt), abs(concordant$NR.exp_trt), paired = TRUE, alternative = "less")

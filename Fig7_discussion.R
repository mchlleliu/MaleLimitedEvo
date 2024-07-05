###################################
#
#                             Grieshop et al. 2023
#                             Author: Michelle Liu
#             DsRed experimental evolution - transcriptomics analysis
#                     Figure 7 - Theoretical discussion figures
# 
# 
###################################


# load packages
#########
library(broom)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggblend)
library(ggpubr)
library(cowplot)
#########



# plot the distributions 
trait_dist <-  ggplot(data.frame(x = c(-950, 950)), aes(x)) + 
  stat_function(fun = dnorm, args = list(mean = 300, sd = 200), 
                geom = 'area', fill = "#D55E00", alpha = 0.5) +
  stat_function(fun = dnorm, args = list(mean = -300, sd =100),
                geom = 'area', fill = "steelblue3", alpha = 0.5) +
  coord_cartesian(ylim = c(0, 0.0045)) +
  theme_classic() +
  theme(plot.title.position = c("panel"),
        legend.title = element_blank(),
        legend.position = c("bottom"),
        legend.text = element_text(size = 25, color = "black"),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size=40, margin = margin(0,0,0,0), color = "black"),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3),
        axis.ticks = element_blank())


# Prediction 1: No conflict
A <- trait_dist +   
  geom_vline(xintercept = -300, color = "steelblue3", linetype = "dashed", size = 1.2) +
  geom_vline(xintercept = 300, color = "#D55E00", linetype = "dashed", size = 1.2) +
  geom_text(aes(x = -940, y = 0.0044), label = expression(bold("A")), size = 10)

A
  

# Prediction 2: there is conflict for further dimorphism
B <- trait_dist +
  geom_vline(xintercept = -450, color = "steelblue3", linetype = "dashed", size = 1.2) +
  geom_vline(xintercept = 700, color = "#D55E00", linetype = "dashed", size = 1.2) +
  geom_segment(aes(x = -300, y = 0.0012, xend = -425, yend = 0.0012),
               lineend = "round", linejoin = "round", color = "steelblue",
               arrow = arrow(length = unit(0.5, "cm")), size = 3.5) +
  geom_segment(aes(x = 300, y = 0.0012, xend = 675, yend = 0.0012),
               lineend = "round", linejoin = "round", color = "#D55E00",
               arrow = arrow(length = unit(0.5, "cm")), size = 3.5) +
  # annotate(geom = "text", label = expression(italic("s")[italic("m")]), y = 0.0015, x = -200, size = 7.5) +
  geom_text(aes(x = -740, y = 0.0044), label = expression(bold("B")), size = 10) +
  coord_cartesian(ylim = c(0, 0.0045), xlim = c(-750, 1000)) 

B  


# Prediction 3: There is conflict for less dimorphism
C <- trait_dist +
  geom_vline(xintercept = -300, color = "steelblue3", linetype = "dashed", size = 1.2) +
  geom_vline(xintercept = 50, color = "#D55E00", linetype = "dashed", size = 1.2) +
  geom_segment(aes(x = 300, y = 0.0012, xend = 75, yend = 0.0012),
               lineend = "round", linejoin = "round", color = "#D55E00",
               arrow = arrow(length = unit(0.5, "cm")), size = 3.5) +
  # geom_segment(aes(x = -500, y = 0.0012, xend = -350, yend = 0.0012),
  #              lineend = "round", linejoin = "round", color = "steelblue",
  #              arrow = arrow(length = unit(0.5, "cm")), size = 3.5) +
  # geom_segment(aes(x = -100, y = 0.0012, xend = -250, yend = 0.0012),
  #              lineend = "round", linejoin = "round", color = "steelblue",
  #              arrow = arrow(length = unit(0.5, "cm")), size = 3.5) +
  geom_text(aes(x = -940, y = 0.0044), label = expression(bold("C")), size = 10)

C  


Fig_7 <- ggarrange(NA, A, NA, 
                   NA, NA, NA,
                   NA, B, NA,
                   NA, NA, NA,
                   NA, C, NA,
                   ncol = 3, nrow = 5,
                   widths = c(0.05, 1, 0.05), 
                   heights = c(1, 0.01, 1, 0.01, 1))
Fig_7



pdf(file = "~/Desktop/UofT/SSAV_RNA/Plots/final_2/Fig7_main.pdf",   # The directory you want to save the file in
    width = 7.5, # 12, 24, 20 The width of the plot in inches
    height = 14) # 10, 20, 13 The height of the plot in inches

annotate_figure(Fig_7 + theme(plot.margin = margin(10,10,10,0)), left = text_grob("Frequency", 
                                            rot = 90, size = 40),
                bottom = text_grob("Trait Value", 
                                   size = 40))
dev.off()

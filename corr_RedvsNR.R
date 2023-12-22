###################################
#
#                             Grieshop et al. 2023
#             DsRed experimental evolution - transcriptomics analysis
#          Correlation of SSAV/Control changes in Red vs NonRed samples
# 
# 
###################################


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
    count(right = .[[x]] > 0, # on the right side of plot?
          top = .[[y]] > 0, # on the top side of plot?
          # for each quadrant, divide into two: 
          # (this is a bit weird, but needed so the numbers can be plotted at the right coordinates)
          UP = !top & abs(.[[x]]) > abs(.[[y]]) | # quadrant III up-left or quadrant IV up-right
            top & abs(.[[x]]) < abs(.[[y]])) %>%  # quadrant I up-right or quadrant II up-left
    mutate(perc = n/sum(n)) %>% # calculate percentage of points relative to total number of points
    
    # this is another strange one for setting up the coordinates
    mutate(conc = right & top | (!right & !top), # quadrant I and quadrant III (the concordant changes)
           dir_UP = conc & UP, # concordant changes where x > y
           dir_DOWN = conc & !UP) %>% # concordant changes where x < y
    
    # TRUE = 1, FALSE = 0
    # specificy coordinates for texts on plot
    mutate(!!x := lim/2*(2*(right - 0.5)+(UP - 0.5)+((conc-0.001)*0.5)-(dir_UP*1.5)+(dir_DOWN*0.5)), 
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
    geom_abline(intercept = 0, slope = 0,  size = 0.5, linetype="solid", color = "black") +
    geom_hline(yintercept = 0,  size = 0.5, linetype="solid", color = "black") +
    geom_vline(xintercept = 0,  size = 0.5, linetype="solid", color = "black") +
    geom_abline(intercept = 0, slope = 1,  size = 0.5, linetype="dashed", color = "black") +
    geom_abline(intercept = 0, slope = -1,  size = 0.5, linetype="dashed", color = "black") +
    
    # add percentages
    geom_text(aes(label = paste(round(perc*100,digits= 2),"%",sep="")), data = quad_n, size = 10) +
    coord_cartesian(xlim=c(-lim, lim), ylim = c(-lim,lim)) +
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

# Set up data frames for plotting
########
tmp.Red <- Red.m.trt
colnames(tmp.Red) <- c("Red.exp_trt", "Red.se_trt", "Red.padj", "FlyBaseID", "Red.Sig")

tmp.NR <- NR.m.trt
colnames(tmp.NR) <- c("NR.exp_trt", "NR.se_trt", "NR.padj", "FlyBaseID", "NR.Sig")

corr.plot <- merge(tmp.Red, tmp.NR, by = "FlyBaseID")
corr.plot <- merge(corr.plot, ASE, by = "FlyBaseID")
########

# plotting based on SBGE categories
########
# you can subset & plot just the candidate genes e.g. = 
# corr.plot[corr.plot$SBGE_comp == "a.more.fbg" & 
#           corr.plot$FlyBaseID %in% A.f.geno[A.f.geno$Sig,]$FlyBaseID

more_fbg <- plot_corr(corr.plot[corr.plot$SBGE_comp == "a.more.fbg",], 
                      "Red.exp_trt", "NR.exp_trt", 
                      "red3", "black", "darkgrey",
                      "SSAV/Control in Red males", 
                      "SSAV/Control in NonRed males", 
                      2.5, "Highly FB (n = 159)")


fbg <- plot_corr(corr.plot[corr.plot$SBGE_comp == "b.fbg",], 
                 "Red.exp_trt", "NR.exp_trt", 
                 "red3", "black", "darkgrey",
                 "SSAV/Control in Red males", 
                 "SSAV/Control in NonRed males", 
                 2.5, "Femaled-biased (n = 3,281)")


ubg <- plot_corr(corr.plot[corr.plot$SBGE_comp == "c.ubg",], 
                 "Red.exp_trt", "NR.exp_trt", 
                 "red3", "black", "darkgrey",
                 "SSAV/Control in Red males", 
                 "SSAV/Control in NonRed males", 
                 2.5, "Unbiased (n = 3,223)")


mbg <- plot_corr(corr.plot[corr.plot$SBGE_comp == "d.mbg",], 
                 "Red.exp_trt", "NR.exp_trt", 
                 "red3", "black", "darkgrey",
                 "SSAV/Control in Red males", 
                 "SSAV/Control in NonRed males", 
                 2.5, "Male-biased (n = 3,110)")


more_mbg <- plot_corr(corr.plot[corr.plot$SBGE_comp == "e.more.mbg",], 
                      "Red.exp_trt", "NR.exp_trt", 
                      "red3", "black", "darkgrey",
                      "SSAV/Control in Red males", 
                      "SSAV/Control in NonRed males", 
                      2.5, "Highly MB (n = 2,372)")
########


pdf(file = "~/Desktop/UofT/SSAV_RNA/Plots/Corr_plots/RedvsNR_bySBGEcat.pdf",   # The directory you want to save the file in
    width = 20, # The width of the plot in inches
    height = 30) # The height of the plot in inches
ggarrange(more_fbg, NA, fbg,
          NA, NA, NA,
          ubg, NA, mbg,
          NA, NA, NA,
          more_mbg, NA, NA,
          labels = c("A)", NA, "B)",
                     NA,NA,NA,
                     "C)", NA, "D)",
                     NA,NA,NA,
                     "E)", NA, NA),
          ncol = 3, nrow = 5,
          widths = c(1, 0.05, 1),
          heights = c(1, 0.05, 1, 0.05, 1),
          font.label = list(size = 30))
dev.off()



####### Binned dot-plot for SBGE (A/C)

# Prepare dataset
##########
Red.m.trt$geno <- "Red"
NR.m.trt$geno <- "NR"
binned.plot <- rbind(Red.m.trt, NR.m.trt)
binned.plot <- merge(binned.plot, ASE, by = "FlyBaseID", all = T)
binned.plot <- merge(binned.plot, Chrs, by = "FlyBaseID")
binned.plot <- binned.plot[!is.na(binned.plot$exp_trt) &
                             !is.na(binned.plot$exp_SBGE_ase) &
                             !is.na(binned.plot$Chr),]

## One sample permutation test for differences n.e to 0 (see boot_permute.R for function)
## by SBGE category
permed_Red.trt <- OnePerm_SBGE(perm_dat = binned.plot[binned.plot$geno == "Red",],
                                x_col = "exp_trt", SBGE_cat = "SBGE_comp")
permed_NR.trt <- OnePerm_SBGE(perm_dat = binned.plot[binned.plot$geno == "NR",],
                               x_col = "exp_trt", SBGE_cat = "SBGE_comp")
permed_Red.trt$geno <- "Red"
permed_NR.trt$geno <- "NR"
permed_RedNR.trt <- rbind(permed_Red.trt, permed_NR.trt)
permed_RedNR.trt <- permed_RedNR.trt %>%
  mutate(holm_padj = p.adjust(permed_RedNR.trt$pval, method = "bonferroni")) %>% 
  mutate(holm_Sig = ifelse(holm_padj < 0.0015, TRUE, FALSE))
rm(permed_Red.trt, permed_NR.trt) # remove clutter

# Permutation to test difference between Red and NonRed
# use temporary dataset for permutation
tmp_RedNR <- binned.plot
tmp_RedNR$IsRed <- tmp_RedNR$geno == "Red"
permDiff_RedNR <- TwoPerm_SBGE(perm_dat = tmp_RedNR, 
                               x_col = "exp_trt", 
                               groupBy = "IsRed", 
                               SBGE_cat = "SBGE_comp")


#########


# Plotting codes
#########
# For 2 groups dot plot with regression line
p <- ggplot(binned.plot, aes(exp_SBGE_ase, exp_trt, color = geno)) +
  geom_point(size = 0.4, shape = 16, alpha = 0.5) +
  geom_smooth(method = "lm",
              formula = y ~ x,
              show.legend = FALSE) +
  labs(x = "SBGE (ASE)", # "omegaA_MK" = expression(italic("\u03c9A")[MK]); "alpha_MK" = expression(italic("\u03b1")[MK])
       y = "Exp. diff. (A/C)") +
  scale_colour_manual(values = c("#888888", "red3"), # "purple3", "chartreuse3", "orange2"
                      labels = c("NR-males", "Red-males")) +
  guides(color = guide_legend(override.aes = list(shape = c(16, 16),
                                                  size = c(4, 4),
                                                  alpha = 1))) +
  geom_abline(intercept = 0, slope = 0,  size = 0.5, linetype="solid", color = "black") +
  scale_x_continuous(limits = c(-5, 5), breaks = c(-1, 0, 1)) +
  scale_y_continuous(limits = c(-1.5, 1.5), breaks = c(-1, -0.5, 0, 0.5, 1)) +
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
        axis.text.x = element_text(size=20, margin = margin(5,0,0,0), color = "black"),
        axis.text.y = element_text(size=20, margin = margin(0,5,0,0), color = "black"), 
        axis.title.x = element_text(size=40, margin = margin(10,0,0,0), color = "black"),
        axis.title.y = element_text(size=40, margin = margin(0,10,0,0), color = "black"),
        plot.title = element_text(size=40, margin = margin(0,0,0,0), color = "black"),
        plot.margin = margin(6,6,6,6)
  )

p



# Binned dot-plot for A-versus-C in males
b <- ggplot(binned.plot, aes(SBGE_comp, exp_trt, color = geno)) +
  geom_point(aes(color = geno), size = 1, shape = 16, alpha = 0.3, 
             position = position_jitterdodge(jitter.width = 0.2)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) + 
  labs(x = "SBGE (ASE)", # "omegaA_MK" = expression(italic("\u03c9A")[MK]); "alpha_MK" = expression(italic("\u03b1")[MK])
       y = "Exp. diff. (A/C)") + 
  scale_colour_manual(values = c("#888888", "red3"),
                      labels = c("NR-males", "Red-males")) +
  guides(color = guide_legend(override.aes = list(shape = c(16, 16),
                                                  size = c(4, 4),
                                                  alpha = 1))) +
  # geom_text(data = permed_RedNR.trt %>% mutate(sig1 = if_else(holm_Sig == 1, "*", "")),
            # aes(x =SBGE_comp, y = 0.95, label = sig1), size = 9, position = position_dodge(width = 0.8), vjust = -4) +
  geom_abline(intercept = 0, slope = 0,  size = 0.5, linetype="solid", color = "black") +
  scale_x_discrete(labels = c("Highly FB", "Female-biased", "Unbiased", "Male-biased", "Highly MB")) + # "Highly FB", "Highly MB", Strongly ant.", "Antagonistic", "Neutral", "Concordant", "Strongly con."
  scale_y_continuous(limits = c(-1.5, 1.5), breaks = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5)) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = c("bottom"),
        #legend.justification = c("right", "bottom"),
        #legend.box.just = "left",
        #legend.box.background = element_rect(),
        legend.box.background = element_rect(),
        #legend.box.margin = margin(4, 6, 6, 6),
        legend.text = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size=20, margin = margin(5,0,0,0), color = "black"),
        axis.text.y = element_text(size=20, margin = margin(0,5,0,0), color = "black"), 
        axis.title.x = element_text(size=40, margin = margin(10,0,0,0), color = "black"),
        axis.title.y = element_text(size=40, margin = margin(0,10,0,0), color = "black"),
        plot.margin = margin(6,6,6,6)
  )



pdf(file = "~/Desktop/UofT/SSAV_RNA/Plots/Corr_plots/RedvsNR_binned.pdf",   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 8) # The height of the plot in inches
b
dev.off()

##########


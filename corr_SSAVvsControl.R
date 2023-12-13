###################################
#
#                             Grieshop et al. 2023
#             DsRed experimental evolution - transcriptomics analysis
#            Correlation of Red/NR changes in SSAV vs Control samples
# 
# 
###################################

# the two functions below are used in the corr plot function below

# quad_count: counts the number of point in each of the 8 quadrants
# variables: dat = data.frame object containing the values to be plotted
#            x = values for group 1, plotted on the x-axis
#            y = values for group 2, plotted on the y-axis
#            lim = for plotting code. specificy the x-lim and y-lim of the plot
quad_count <- function(dat, x, y, lim){
  count <- dat %>%
    # Count how many with each combination of X and Y being positive
    count(right = .[[x]] > 0, # on the right side of plot?
          top = .[[y]] > 0, # on the top side of plot?
          # for each quadrant, divide into two: 
          # (this is a bit weird, but needed so the numbers can be plotted at the right coordinates)
          UP = !top & abs(.[[x]]) > abs(.[[y]]) | # quadrant III up-left or quadrant IV up-right
            top & abs(.[[x]]) < abs(.[[y]])) %>%  # quadrant I up-right or quadrant II up-left
    mutate(perc = n/sum(n)) %>% # calculate percentage of 
    mutate(conc = right & top | (!right & !top), dir_UP = conc & UP, dir_DOWN = conc & !UP) %>%
    # TRUE = 1, FALSE = 0
    # specificy coordinates for texts on plot
    mutate(!!x := lim/2*(2*(right - 0.5)+(UP - 0.5)+((conc-0.001)*0.5)-(dir_UP*1.5)+(dir_DOWN*0.5)), 
           !!y := lim/2*(2*(top - 0.5)+(UP - 0.5)))
  return(count)
}


colour_quadrant <-  function(dat, x, y, colx, coly){
  col <- dat %>%
    # Count how many with each combination of X and Y being positive
    mutate(right = .[[x]] > 0, 
           top = .[[y]] > 0,
           DOWN = !top & abs(.[[x]]) > abs(.[[y]]) |
             top & abs(.[[x]]) > abs(.[[y]])) %>%
    mutate(quadrant = ifelse(right & top | (!right & !top), 
                             ifelse(DOWN, colx, coly), 
                             "red3"))
  return(col)
}

plot_corr <- function(dat, x, y, colx, coly, xlab, ylab, lim, title){
  pear_cor <- cor.test(dat[[x]], dat[[y]],method = "pearson")
  quad_n <- quad_count(dat, x, y, lim)
  quad_col <- colour_quadrant(dat, x, y, colx, coly)
  corr <- ggplot(dat, aes_string(x = x, y = y)) +
    geom_point(size = 2, shape = 16, alpha = 0.5, color = quad_col$quadrant) +  
    geom_abline(intercept = 0, slope = 0,  size = 0.5, linetype="solid", color = "black") +
    geom_hline(yintercept = 0,  size = 0.5, linetype="solid", color = "black") +
    geom_vline(xintercept = 0,  size = 0.5, linetype="solid", color = "black") +
    # geom_abline(intercept = 0, slope = 1,  size = 0.5, linetype="dashed", color = "black") +
    # geom_abline(intercept = 0, slope = -1,  size = 0.5, linetype="dashed", color = "black") +
    coord_cartesian(xlim=c(-lim, lim), ylim = c(-lim,lim)) +
    # geom_text(aes(label = paste(round(perc*100,digits= 2),"%",sep="")), data = quad_n, size = 10) +
    geom_label(aes(x = 1.5, y = -1.5, label = c(paste("r:", round(pear_cor$estimate, digits = 3), "\n[", round(pear_cor$conf.int[1], digits = 3),
                                                      ",", round(pear_cor$conf.int[2], digits = 3), "]", sep = " "))), data = quad_n, size = 10) +
    labs(x = print(xlab), # "omegaA_MK" = expression(italic("\u03c9A")[MK]); "alpha_MK" = expression(italic("\u03b1")[MK])
         y = print(ylab) ,
         title = print(title)) +
    # scale_colour_manual(values = c("#888888", "red3"), # "purple3", "chartreuse3", "orange2" # "red3", "#888888", "steelblue3"
    #                     labels = c("Up-reg NR", "Up-reg Red")) + # "FBG", "MBG", "UBG" # "Chr-2", "Chr-3", "X-Chr"
    guides(color = guide_legend(override.aes = list(shape = c(NA, NA), # c(16, 16)
                                                    size = c(4, 4),
                                                    alpha = 1))) +
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
  return(corr)
}


AmCmf <- plot_corr(corr.plot[corr.plot$FlyBaseID %in% A.f.geno[A.f.geno$Sig,]$FlyBaseID, ], "m.exp_geno", "C.m.exp_geno", 
                   "black", "black", "Red/NR in SSAV males", "Red/NR in Control males", 2.5, "")
AmCmm <- plot_corr(corr.plot[corr.plot$FlyBaseID %in% tmp.males[tmp.males$m.Sig,]$FlyBaseID, ], 
                   "m.exp_geno", "C.m.exp_geno", "black", "black", "Red/NR in SSAV males", 
                   "Red/NR in Control males", 2.5, "")

AmAff <- plot_corr(corr.plot[corr.plot$FlyBaseID %in% A.f.geno[A.f.geno$Sig,]$FlyBaseID, ], "f.exp_geno", "m.exp_geno", 
                   "black", "black", "Red/NR in SSAV females", "Red/NR in SSAV males", 2.5, "")
AmAfm <- plot_corr(corr.plot[corr.plot$FlyBaseID %in% A.m.geno[A.m.geno$Sig,]$FlyBaseID, ], "m.exp_geno", "f.exp_geno", 
                   "black", "black", "Red/NR in SSAV males", "Red/NR in SSAV females", 2.5, "")

AfCmf <- plot_corr(corr.plot[corr.plot$FlyBaseID %in% A.f.geno[A.f.geno$Sig,]$FlyBaseID, ], "f.exp_geno", "C.m.exp_geno", 
                   "black", "black", "Red/NR in SSAV females", "Red/NR in Control males", 2.5, "")
AfCmm <- plot_corr(corr.plot[corr.plot$FlyBaseID %in% A.m.geno[A.m.geno$Sig,]$FlyBaseID, ], "f.exp_geno", "C.m.exp_geno", 
                   "black", "black", "Red/NR in SSAV females", "Red/NR in Control males", 2.5, "")

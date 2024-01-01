###################################
#
#                             Grieshop et al. 2023
#                             Author: Michelle Liu
#             DsRed experimental evolution - transcriptomics analysis
#                  Distribution of SA candidates by chromosomes
# 
# 
###################################

# Prepare plotting dataset
# load results if not loaded in env.
A.f.geno <- read.delim("Results/A.f.geno_candidates.tsv")
A.m.geno <- read.delim("Results/A.m.geno_candidates.tsv")

# include SBGE categories (using Mishra et al. dataset. Look at External_data.R)
A.m.geno_Chr <- merge(A.m.geno, Chrs_All, by = "FlyBaseID", all = TRUE)
A.f.geno_Chr <- merge(A.f.geno, Chrs_All, by = "FlyBaseID", all = TRUE)



## binned proportion of candidate genes
propChr <- function(dat){
  total_All <- dim(dat)[1]
  total_Sig <- dim(dat[dat$Sig,])[1]
  total_NS <- total_All - total_Sig
  
  total_Chr <- dat %>% group_by(Chr) %>%
    summarise(Chr_count = n()) 
  
  fract <- dat %>% group_by(Chr, Sig) %>%
    summarise(count = n()) %>%
    mutate(frac_Sig = ifelse(Sig, count/total_Sig, count/total_NS))
  
  fract <- merge(fract, total_Chr, by = "Chr") %>%
    rowwise() %>%
    mutate(frac_by_Chr = count/Chr_count,
           lower = list(binom.test(count, Chr_count)),
           upper = lower$conf.int[2],
           lower = lower$conf.int[1])

  return(fract)
}
plotChrprop <- function(dat, xlab){
  plot_dat <- propChr(dat)
  
  prop_all <- data.frame(dim(dat[dat$Sig,])[1]/dim(dat)[1])
  prop_all$lower <- prop.test(dim(dat[dat$Sig,])[1], dim(dat)[1])$conf.int[1]
  prop_all$upper <- prop.test(dim(dat[dat$Sig,])[1], dim(dat)[1])$conf.int[2]
  
  plot_vec <- ggplot(plot_dat[plot_dat$Sig,], aes_string(x = "Chr", y = "frac_by_Chr")) + 
    geom_errorbar(ymin = plot_dat[plot_dat$Sig,]$lower, ymax = plot_dat[plot_dat$Sig,]$upper,
                  width = 0.5, size = 0.75) +
    geom_point(fill = "forestgreen", color = "forestgreen", size = 7) + 
    labs(x = print(xlab), # "omegaA_MK" = expression(italic("\u03c9A")[MK]); "alpha_MK" = expression(italic("\u03b1")[MK])
         y = "fraction of candidate genes") +
    scale_y_continuous(limits = c(0, 0.20)) +
    geom_hline(yintercept = prop_all[,1],
               linetype = "dashed", show.legend = TRUE, size = 0.6) +
    annotate('ribbon', x = c(-Inf, Inf), ymin = prop_all$lower, ymax = prop_all$upper, 
             alpha = 0.15, fill = 'forestgreen') +
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
          axis.title.x = element_text(size=30, margin = margin(10,0,0,0), color = "black"),
          axis.title.y = element_text(size=30, margin = margin(0,10,0,0), color = "black"),
          plot.title = element_text(size=40, margin = margin(0,0,0,0), color = "black"),
          plot.margin = margin(6,6,6,6))
  return(plot_vec)
}

A.f_Sig <- plotChrprop(A.f.geno_Chr[!is.na(A.f.geno_Chr$Sig) &
                                       !is.na(A.f.geno_Chr$Chr),], "SBGE (ASE)")
A.m_Sig <- plotChrprop(A.m.geno_Chr[!is.na(A.m.geno_Chr$Sig) &
                           !is.na(A.m.geno_Chr$Chr) & A.m.geno_Chr$Chr != "Y",], 
            "SBGE (ASE)") + coord_cartesian(ylim = c(0, 0.05))



pdf(file = "~/Desktop/UofT/SSAV_RNA/Plots/Chromosomal_dist.pdf",   # The directory you want to save the file in
    width = 24, # 12, 24 The width of the plot in inches
    height = 10) # 10, 20 The height of the plot in inches
ggarrange(A.f_Sig, NA, A.m_Sig,
          labels = c("A)", NA, "B)"),
          widths = c(1, 0.05, 1),
          ncol = 3,
          font.label = list(size = 30))
dev.off()





###################################
#
#                             Grieshop et al. 2023
#             DsRed experimental evolution - transcriptomics analysis
#                  Direction of Red/NR changes by SBGE catergory
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
#########

# Prepare plotting dataset containing all samples
#########
# merge dataset with out population log2FC estimates. 
## Using Mishra et al.'s data. See "External_data.R" for code to bin SBGE categories
A.m.geno_ASE <- merge(A.m.geno, ASE, by = "FlyBaseID", all = T)
A.f.geno_ASE <- merge(A.f.geno, ASE, by = "FlyBaseID", all = T)
C.m.geno_ASE <- merge(C.m.geno, ASE, by = "FlyBaseID", all = T)
A.m.geno_ASE$trt2 = "Am"
A.f.geno_ASE$trt2 = "Af"
C.m.geno_ASE$trt2 = "Cm"

# only keep genes with data available in both datasets
A.m.geno_ASE <- A.m.geno_ASE[!is.na(A.m.geno_ASE$exp_geno) &
                               !is.na(A.m.geno_ASE$exp_SBGE_ase),]
A.f.geno_ASE <- A.f.geno_ASE[!is.na(A.f.geno_ASE$exp_geno) &
                               !is.na(A.f.geno_ASE$exp_SBGE_ase),]
C.m.geno_ASE <- C.m.geno_ASE[!is.na(C.m.geno_ASE$exp_geno) &
                               !is.na(C.m.geno_ASE$exp_SBGE_ase),]
dim(C.m.geno_ASE) # check how many cut off

# merge all to one dataframe object
All.geno <- rbind(A.m.geno_ASE[-5], A.f.geno_ASE, C.m.geno_ASE)
#########


# Permute datasets
########
## One sample permutation test for differences n.e to 0 (see boot_permute.R for function)
## by SBGE category
permed_A.f.geno <- OnePerm_SBGE(perm_dat = All.geno[All.geno$trt2 == "Af",],
                                x_col = "exp_geno", SBGE_cat = "SBGE_comp")
permed_A.m.geno <- OnePerm_SBGE(perm_dat = All.geno[All.geno$trt2 == "Am",],
                                x_col = "exp_geno", SBGE_cat = "SBGE_comp")
permed_C.m.geno <- OnePerm_SBGE(perm_dat = All.geno[All.geno$trt2 == "Cm",],
                                x_col = "exp_geno", SBGE_cat = "SBGE_comp")
permed_A.f.geno$trt2 <- "Af"
permed_A.m.geno$trt2 <- "Am"
permed_C.m.geno$trt2 <- "Cm"
permed_All.geno <- rbind(permed_A.f.geno, permed_A.m.geno, permed_C.m.geno)
permed_All.geno <- permed_All.geno %>%
  mutate(holm_padj = p.adjust(permed_All.geno$pval, method = "holm")) %>% 
  mutate(holm_Sig = ifelse(holm_padj < 0.005, TRUE, FALSE))
rm(permed_A.f.geno, permed_A.m.geno, permed_C.m.geno) # remove clutter


# only candidate genes found in SSAV females
permed_A.f.geno_AfSig <- OnePerm_SBGE(perm_dat = All.geno[All.geno$trt2 == "Af" &
                                                      All.geno$FlyBaseID %in% A.f.geno[A.f.geno$Sig,]$FlyBaseID,],
                                x_col = "exp_geno", SBGE_cat = "SBGE_comp")
permed_A.m.geno_AfSig <- OnePerm_SBGE(perm_dat = All.geno[All.geno$trt2 == "Am" &
                                                      All.geno$FlyBaseID %in% A.f.geno[A.f.geno$Sig,]$FlyBaseID,],
                                x_col = "exp_geno", SBGE_cat = "SBGE_comp")
permed_C.m.geno_AfSig <- OnePerm_SBGE(perm_dat = All.geno[All.geno$trt2 == "Cm" &
                                                      All.geno$FlyBaseID %in% A.f.geno[A.f.geno$Sig,]$FlyBaseID,],
                                x_col = "exp_geno", SBGE_cat = "SBGE_comp")
permed_A.f.geno_AfSig$trt2 <- "Af"
permed_A.m.geno_AfSig$trt2 <- "Am"
permed_C.m.geno_AfSig$trt2 <- "Cm"
permed_All.geno_AfSig <- rbind(permed_A.f.geno_AfSig, permed_A.m.geno_AfSig, permed_C.m.geno_AfSig)
permed_All.geno_AfSig <- permed_All.geno_AfSig %>%
  mutate(holm_padj = p.adjust(permed_All.geno_AfSig$pval, method = "holm")) %>% 
  mutate(holm_Sig = ifelse((n > 30 & holm_padj < 0.05), TRUE, FALSE))
rm(permed_A.f.geno_AfSig, permed_A.m.geno_AfSig, permed_C.m.geno_AfSig) # remove clutter



permed_A.f.geno_AmSig <- OnePerm_SBGE(perm_dat = All.geno[All.geno$trt2 == "Af" &
                                                            All.geno$FlyBaseID %in% A.m.geno[A.m.geno$Sig,]$FlyBaseID,],
                                      x_col = "exp_geno", SBGE_cat = "SBGE_comp")
permed_A.m.geno_AmSig <- OnePerm_SBGE(perm_dat = All.geno[All.geno$trt2 == "Am" &
                                                            All.geno$FlyBaseID %in% A.m.geno[A.m.geno$Sig,]$FlyBaseID,],
                                      x_col = "exp_geno", SBGE_cat = "SBGE_comp")
permed_C.m.geno_AmSig <- OnePerm_SBGE(perm_dat = All.geno[All.geno$trt2 == "Cm" &
                                                            All.geno$FlyBaseID %in% A.m.geno[A.m.geno$Sig,]$FlyBaseID,],
                                      x_col = "exp_geno", SBGE_cat = "SBGE_comp")
permed_A.f.geno_AmSig$trt2 <- "Af"
permed_A.m.geno_AmSig$trt2 <- "Am"
permed_C.m.geno_AmSig$trt2 <- "Cm"
permed_All.geno_AmSig <- rbind(permed_A.f.geno_AmSig, permed_A.m.geno_AmSig, permed_C.m.geno_AmSig)
permed_All.geno_AmSig <- permed_All.geno_AmSig %>%
  mutate(holm_padj = p.adjust(permed_All.geno_AmSig$pval, method = "holm")) %>% 
  mutate(holm_Sig = ifelse((n > 30 & holm_padj < 0.05), TRUE, FALSE))
rm(permed_A.f.geno_AmSig, permed_A.m.geno_AmSig, permed_C.m.geno_AmSig) # remove clutter

########




## plotting codes:
# store in object: All.exp_geno, A.f.sig.exp_geno, A.m.sig.exp_geno
# Binned dot-plot
########
binPlot_RedNR <- function(dat, perm_dat){
  ggplot(dat, aes(SBGE_comp, exp_geno, color = trt2)) +
  geom_point(aes(color = trt2), size = 1, shape = 16, 
             alpha = 0.3, position = position_jitterdodge(jitter.width = 0.17)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) + 
  labs(x = "SBGE (ASE)", # "omegaA_MK" = expression(italic("\u03c9A")[MK]); "alpha_MK" = expression(italic("\u03b1")[MK])
       y = "Exp. diff. (Red/NR)") +
  # title = "C-males") +
  scale_colour_manual(values = c("red3", "steelblue3", "#888888"), # "red3", "steelblue3", "#888888" # "purple3", "chartreuse3", "orange2"
                      labels = c("SSAV females", "SSAV males", "Control males")) + # "Chr-2", "Chr-3", "X-Chr"
  guides(color = guide_legend(override.aes = list(shape = c(16, 16, 16),
                                                  size = c(4, 4, 4),
                                                  alpha = 1))) +
  # add star do signify significant difference from 0
  geom_text(data = perm_dat %>% mutate(sig1 = if_else(holm_Sig , "*", "")),
            aes(x =SBGE_comp, y = 0.95, label = sig1), 
            size = 10, position = position_dodge(width = 0.8), vjust = -3.5, show.legend = FALSE) +
  # add number of genes per category
  geom_text(data = perm_dat, aes(label = n, y = Inf, group = trt2), 
            position = position_dodge(width = 0.8), vjust = 6, size =4.5, show.legend = FALSE) +
  # line at y = 0
  geom_abline(intercept = 0, slope = 0,  size = 0.5, linetype= "solid", color = "black") +
  scale_x_discrete(labels = c("Highly FB", "Female-Biased", "Unbiased", "Male-Biased", "Highly MB")) + # "Highly ant.", "Antagonistic", "Uncorrelated", "Concordant", "Highly con." .... "Strong pur.", "Purifying sel.", "Neutral", "Positive sel.", "Strong pos."
  scale_y_continuous(limits = c(-2, 2), breaks = c(-2, -1, 0, 1, 2)) +
  # default theme settings:
  theme_classic() +
  theme(plot.title.position = c("panel"),
        legend.title = element_blank(),
        legend.position = c("bottom"),
        legend.text = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size=20, margin = margin(5,0,0,0), color = "black"),
        axis.text.y = element_text(size=20, margin = margin(0,5,0,0), color = "black"),
        axis.title.x = element_text(size=30, margin = margin(10,0,0,0), color = "black"),
        axis.title.y = element_text(size=30, margin = margin(0,10,0,0), color = "black"),
        plot.title = element_text(size=40, margin = margin(0,0,0,0), color = "black"),
        plot.margin = margin(6,6,6,6)
  )
}

All.exp_geno <- binPlot_RedNR(All.geno, permed_All.geno)
A.f.sig.exp_geno <- binPlot_RedNR(All.geno[All.geno$FlyBaseID %in% A.f.geno[A.f.geno$Sig,]$FlyBaseID,],
                                  permed_All.geno_AfSig)
A.m.sig.exp_geno <- binPlot_RedNR(All.geno[All.geno$FlyBaseID %in% A.m.geno[A.m.geno$Sig,]$FlyBaseID,],
                                  permed_All.geno_AmSig)


pdf(file = "~/Desktop/UofT/SSAV_RNA/Plots/RedvsNR_All.pdf",  # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 17) # The height of the plot in inches
ggarrange(All.exp_geno + theme(axis.title.x = element_blank(), legend.position = c("none")),
          NA, 
          A.f.sig.exp_geno + theme(axis.title.x = element_blank(), legend.position = c("none")), 
          NA, 
          A.m.sig.exp_geno,         
          labels = c("A)", NA, "B)", NA, "C)"),
          heights = c(1, 0.05, 1, 0.05, 1.25), ncol =1, nrow = 5, 
          font.label = list(size = 30)) 
dev.off()
########


# LOESS plots
##########
plotLoess <- function(dat_loess, sig_genes = NA){
  if(!is.null(dim(sig_genes))){
    plot_dat <- dat_loess[dat_loess$FlyBaseID %in% sig_genes$FlyBaseID,]
  }
  else{
    plot_dat <- dat_loess
  }
  min_x <- 5 * round((min(plot_dat$exp_SBGE_ase) - 1) / 5)
  max_x <-  5 * round((max(plot_dat$exp_SBGE_ase) + 1) / 5)
  loess <- ggplot(data = plot_dat, aes(exp_SBGE_ase, exp_geno, color=trt2)) + 
    geom_point(size = 2, shape = 16, alpha = 0.7, colour = "grey") +
    geom_smooth(method = "loess", span = 0.5) +
    scale_color_manual(name=NULL, values=c("red3","blue3", "darkgreen")) +
    # scale_x_continuous(limits = c(min_x, max_x), breaks = seq(min_x, max_x, by = 2.5)) +
    labs(y="exp. diff (Red/NR)", x="SBGE(ASE)") +
    geom_hline(yintercept = 0, size = 0.5, linetype= "dashed", color = "black", alpha = 0.7) +
    geom_vline(xintercept = -5, size = 0.5, linetype= "dashed", color = "black", alpha = 0.7) +
    geom_vline(xintercept = 0, size = 0.5, linetype= "dashed", color = "black", alpha = 0.7) +
    geom_vline(xintercept = 5, size = 0.5, linetype= "dashed", color = "black", alpha = 0.7) +
    theme_classic() +
    theme( # legend.title = element_blank(),
      #       legend.position = c("None"),
      #       #legend.justification = c("right", "bottom"),
      #       #legend.box.just = "left",
      #       #legend.box.background = element_rect(),
      #       legend.box.background = element_rect(),
      #       #legend.box.margin = margin(4, 6, 6, 6),
      legend.text = element_text(size = 20, color = "black"),
      axis.text.x = element_text(size=20, margin = margin(5,0,0,0), color = "black"),
      axis.text.y = element_text(size=20, margin = margin(0,5,0,0), color = "black"), 
      axis.title.x = element_text(size=40, margin = margin(10,0,0,0), color = "black"),
      axis.title.y = element_text(size=40, margin = margin(0,10,0,0), color = "black"),
      plot.title = element_text(size=40, margin = margin(0,0,0,0), color = "black"),
      plot.margin = margin(6,6,6,6)) + coord_cartesian(xlim = c(-10, 10))
  return(loess)
}


loess_All.geno_Af_sig <- plotLoess(All.geno, A.f.geno[A.f.geno$Sig==TRUE,])
loess_All.geno_Am_sig <- plotLoess(All.geno, A.m.geno[A.m.geno$Sig==TRUE,])
loess_All.geno <- plotLoess(All.geno)

##########


# Looking at standard error distribution
######
dens_se_plot <- ggplot() + 
  geom_density(data=All.geno, aes(se_geno, fill= trt2, color = trt2), alpha = 0.2) +
  labs(x="SE of log2FC") +
  scale_colour_manual(values = c("darkred", "darkblue", "black"), # "red3", "steelblue3", "#888888" # "purple3", "chartreuse3", "orange2"
                      labels = c("SSAV females", "SSAV males", "Control males")) +  # "Chr-2", "Chr-3", "X-Chr"
  scale_fill_manual(values = c("red3", "steelblue3", "#888888"), # "red3", "steelblue3", "#888888" # "purple3", "chartreuse3", "orange2"
                    labels = c("SSAV females", "SSAV males", "Control males")) +  # "Chr-2", "Chr-3", "X-Chr"
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = c(0.75, 0.85),
        #       #legend.justification = c("right", "bottom"),
        #       #legend.box.just = "left",
        #       #legend.box.background = element_rect(),
        #       legend.box.background = element_rect(),
        #       #legend.box.margin = margin(4, 6, 6, 6),
        legend.text = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size=20, margin = margin(5,0,0,0), color = "black"),
        axis.text.y = element_text(size=20, margin = margin(0,5,0,0), color = "black"), 
        axis.title.x = element_text(size=40, margin = margin(10,0,0,0), color = "black"),
        axis.title.y = element_text(size=40, margin = margin(0,10,0,0), color = "black"),
        plot.title = element_text(size=40, margin = margin(0,0,0,0), color = "black"),
        plot.margin = margin(6,6,6,6)
  )
######


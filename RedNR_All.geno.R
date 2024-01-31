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
library(cowplot)
#########

# Prepare plotting dataset containing all samples
#########
# load datasets
A.f.geno <- read.delim("Results/A.f.geno_candidates.tsv")
A.m.geno <- read.delim("Results/A.m.geno_candidates.tsv")
C.m.geno <- read.delim("Results/C.m.geno_candidates.tsv")

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
All.geno$SBGE_simp <- as.factor(All.geno$SBGE_simp)
All.geno$SBGE_comp <- as.factor(All.geno$SBGE_comp)
All.geno$trt2 <- as.factor(All.geno$trt2)
str(All.geno)
#########


# Permute against 0
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
write.table(permed_All.geno, file = "~/Desktop/UofT/SSAV_RNA/Results/permed_All.geno.tsv", sep = "\t", # Fix file name accordingly
            row.names = FALSE, col.names = TRUE)



# all candidate genes
permed_A.f.geno_Sig <- OnePerm_SBGE(perm_dat = All.geno[All.geno$trt2 == "Af" &
                                                          All.geno$FlyBaseID %in% SSAV.geno[SSAV.geno$Sig,]$FlyBaseID,],
                                x_col = "exp_geno", SBGE_cat = "SBGE_comp")
permed_A.m.geno_Sig <- OnePerm_SBGE(perm_dat = All.geno[All.geno$trt2 == "Am" &
                                                          All.geno$FlyBaseID %in% SSAV.geno[SSAV.geno$Sig,]$FlyBaseID,],
                                x_col = "exp_geno", SBGE_cat = "SBGE_comp")
permed_C.m.geno_Sig <- OnePerm_SBGE(perm_dat = All.geno[All.geno$trt2 == "Cm" &
                                                          All.geno$FlyBaseID %in% SSAV.geno[SSAV.geno$Sig,]$FlyBaseID,],
                                x_col = "exp_geno", SBGE_cat = "SBGE_comp")
permed_A.f.geno_Sig$trt2 <- "Af"
permed_A.m.geno_Sig$trt2 <- "Am"
permed_C.m.geno_Sig$trt2 <- "Cm"
permed_All.geno_Sig <- rbind(permed_A.f.geno_Sig, permed_A.m.geno_Sig, permed_C.m.geno_Sig)
permed_All.geno_Sig <- permed_All.geno_Sig %>%
  mutate(holm_padj = p.adjust(permed_All.geno_Sig$pval, method = "holm")) %>% 
  mutate(holm_Sig = ifelse(holm_padj < 0.005, TRUE, FALSE))
rm(permed_A.f.geno_Sig, permed_A.m.geno_Sig, permed_C.m.geno_Sig) # remove clutter

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


# Permute SSAV males against Control males
########
perm_A.m_C.m <- TwoPerm_SBGE(perm_dat = All.geno[All.geno$trt2 != "Af",] %>%
                               mutate(Am = ifelse(trt2 == "Am", TRUE, FALSE)), 
                             x_col = "exp_geno", groupBy = "Am", SBGE_cat = "SBGE_comp")

perm_A.m_C.m_Sig <- TwoPerm_SBGE(perm_dat = All.geno[All.geno$trt2 != "Af" &
                                                       All.geno$FlyBaseID %in% SSAV.geno[SSAV.geno$Sig,]$FlyBaseID,] %>%
                               mutate(Am = ifelse(trt2 == "Am", TRUE, FALSE)), 
                             x_col = "exp_geno", groupBy = "Am", SBGE_cat = "SBGE_comp")

perm_A.m.C.m_FemSig <- TwoPerm_SBGE(perm_dat = All.geno[All.geno$trt2 != "Af" &
                                                   All.geno$FlyBaseID %in% A.f.geno[A.f.geno$Sig,]$FlyBaseID,] %>%
                               mutate(Am = ifelse(trt2 == "Am", TRUE, FALSE)), 
                             x_col = "exp_geno", groupBy = "Am", SBGE_cat = "SBGE_comp")
perm_A.m.C.m_MaleSig <- TwoPerm_SBGE(perm_dat = All.geno[All.geno$trt2 != "Af" &
                                                          All.geno$FlyBaseID %in% A.m.geno[A.m.geno$Sig,]$FlyBaseID,] %>%
                                      mutate(Am = ifelse(trt2 == "Am", TRUE, FALSE)), 
                                    x_col = "exp_geno", groupBy = "Am", SBGE_cat = "SBGE_comp")
########


## plotting codes:
# store in object: All.exp_geno, A.f.sig.exp_geno, A.m.sig.exp_geno
# Binned dot-plot
########
binPlot_RedNR <- function(dat, perm_dat){
  ggplot(dat, aes(SBGE_comp, exp_geno, color = trt2)) +
  geom_point(aes(color = trt2), size = 1, shape = 16, 
             alpha = 0.3, position = position_jitterdodge(jitter.width = 0.35)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) + 
  labs(x = "Sex-biased Gene Expression", # "omegaA_MK" = expression(italic("\u03c9A")[MK]); "alpha_MK" = expression(italic("\u03b1")[MK])
       y = expression(Log["2"]*"FC (Red/NR)")) +
  # title = "C-males") +
  scale_colour_manual(values = c("red3", "steelblue3", "#888888"), # "red3", "steelblue3", "#888888" # "purple3", "chartreuse3", "orange2"
                      labels = c("SM females", "SM males", "C males")) + # "Chr-2", "Chr-3", "X-Chr"
  guides(color = guide_legend(override.aes = list(shape = c(18, 18, 18),
                                                  size = c(5, 5, 5),
                                                  alpha = 1))) +
  # add star do signify significant difference from 0
  geom_text(data = perm_dat %>% mutate(sig1 = if_else(holm_Sig , "*", "ns")),
            aes(x =SBGE_comp, y = 2, label = sig1), size = 7.5,
            position = position_dodge(width = 0.8), show.legend = FALSE) +
  # add number of genes per category
  geom_text(data = perm_dat, aes(label = n, y = Inf, group = trt2), 
            position = position_dodge(width = 0.8), vjust = 5, size = 6, show.legend = FALSE) +
  # line at y = 0
  geom_abline(intercept = 0, slope = 0,  size = 0.5, linetype= "solid", color = "black") +
    geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5), color = "grey") +
  scale_x_discrete(labels = c("Highly FB", "Female-Biased", "Unbiased", "Male-Biased", "Highly MB")) + # "Highly ant.", "Antagonistic", "Uncorrelated", "Concordant", "Highly con." .... "Strong pur.", "Purifying sel.", "Neutral", "Positive sel.", "Strong pos."
  scale_y_continuous(limits = c(-2, 2), breaks = c(-2, -1, 0, 1, 2)) +
  # default theme settings:
  theme_classic() +
  theme(plot.title.position = c("panel"),
        legend.title = element_blank(),
        legend.position = c("bottom"),
        legend.text = element_text(size = 25, color = "black"),
        axis.text.x = element_text(size=20, margin = margin(5,0,0,0), color = "black"),
        axis.text.y = element_text(size=20, margin = margin(0,5,0,0), color = "black"),
        axis.title.x = element_text(size=30, margin = margin(10,0,0,0), color = "black"),
        axis.title.y = element_text(size=30, margin = margin(0,10,0,0), color = "black"),
        plot.title = element_text(size=40, margin = margin(0,0,0,0), color = "black"),
        plot.margin = margin(6,6,6,6),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
  ) 
}

All.exp_geno <- binPlot_RedNR(All.geno, permed_All.geno) + 
  geom_signif(y_position = c(-1.2, -1.2, -1.3, -1.7, -1.35), xmin = c(1, 2, 3, 4, 5), 
              xmax = c(1.3, 2.3, 3.3, 4.3, 5.3),
  annotation = c("*", "*", "*", "*", "*"), tip_length = -0.01, 
  textsize = 10, size = 0.75, vjust = 1.85, color = "darkblue")

All.sig.exp_geno <- binPlot_RedNR(All.geno[All.geno$FlyBaseID %in% SSAV.geno[SSAV.geno$Sig,]$FlyBaseID,], 
                                  permed_All.geno_Sig) +
  geom_signif(y_position = c(-1.18, -1.1, -1.25), xmin = c(2, 3, 4), 
              xmax = c(2.3, 3.3, 4.3),
              annotation = c("*", "*", "*"), tip_length = -0.01, 
              textsize = 10, size = 0.75, vjust = 1.8, color = "darkblue")

A.f.sig.exp_geno <- binPlot_RedNR(All.geno[All.geno$FlyBaseID %in% A.f.geno[A.f.geno$Sig,]$FlyBaseID,],
                                  permed_All.geno_AfSig) +
  geom_signif(y_position = c(-1.1, -1.1, -1.25), xmin = c(2, 3, 4), 
              xmax = c(2.3, 3.3, 4.3),
              annotation = c("*", "*", "*"), tip_length = -0.01, 
              textsize = 10, size = 0.75, vjust = 1.8, color = "darkblue")

A.m.sig.exp_geno <- binPlot_RedNR(All.geno[All.geno$FlyBaseID %in% A.m.geno[A.m.geno$Sig,]$FlyBaseID,],
                                  permed_All.geno_AmSig) +
  geom_signif(y_position = c(-1.15, -1, -1.2), xmin = c(2, 3, 4), 
              xmax = c(2.3, 3.3, 4.3),
              annotation = c("*", "*", "*"), tip_length = -0.01, 
              textsize = 10, size = 0.75, vjust = 1.8, color = "darkblue")


pdf(file = "~/Desktop/UofT/SSAV_RNA/Plots/finals/Fig3_main.pdf",  # The directory you want to save the file in
    width = 15, # 15 The width of the plot in inches
    height = 8 ) # 8 20 The height of the plot in inches
# ggarrange(All.sig.exp_geno + theme(axis.title.x = element_blank(), axis.text.x = element_blank()),
#           NA,
#           A.m.sig.exp_geno + theme(axis.title.x = element_blank(), axis.text.x = element_blank()),
#           NA, A.f.sig.exp_geno, NA,
#           labels = c("A)", NA, "B)", NA, "C)", NA),
#           heights = c(1, 0.05, 1, 0.05, 1, 0.01), ncol =1, nrow = 6,
#           font.label = list(size = 25),
#           common.legend = TRUE, legend = "bottom")
All.exp_geno
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



# Hist X binned dotplot
########

# for number of genes per category
# dim(A.f.geno_ASE[!is.na(A.f.geno_ASE$exp_SBGE_ase) & A.f.geno_ASE$SBGE_simp == "c.mbg",])
# A.f.geno = "Female-biased (3,493)", "Unbiased (3,087)", "Male-biased (2,501)"
# A.m.geno = "Female-biased (3,448)", "Unbiased (3,226)", "Male-biased (5,499)"
# C.m.geno = "Female-biased (3,451)", "Unbiased (3,224)", "Male-biased (5,485)"

# histogram
hist_A.f.geno
hist_A.m.geno
hist_C.m.geno <- ggplot() + 
  geom_blank(data=C.m.geno_ASE[!is.na(C.m.geno_ASE$exp_SBGE_ase),],
                 aes(x=exp_geno, fill = SBGE_simp, color = SBGE_simp)) +
  geom_histogram(data=C.m.geno_ASE[C.m.geno_ASE$SBGE_simp != "b.ubg",], 
                 aes(x=exp_geno, y=..count../sum(..count..),
                     fill=factor(SBGE_simp, levels = c("c.mbg", "b.ubg", "a.fbg")), 
                     color = factor(SBGE_simp, levels = c("c.mbg", "b.ubg", "a.fbg"))),
                 alpha = 0.7, position = "identity", binwidth = 0.05) +
  labs(y="proportion of genes") + 
  scale_colour_manual(values = c("red3", "#888888", "steelblue3"), # "purple3", "chartreuse3", "orange2"
                      labels = c("Female-biased (3,451)", "Unbiased (3,224)", "Male-biased (5,485)")) +  # "Chr-2", "Chr-3", "X-Chr"
  scale_fill_manual(values = c("red3", "#888888", "steelblue3"), # "purple3", "chartreuse3", "orange2"
                    labels = c("Female-biased (3,451)", "Unbiased (3,224)", "Male-biased (5,485)")) +  # "Chr-2", "Chr-3", "X-Chr"
  geom_vline(xintercept = 0, size = 0.5, linetype= "dashed", color = "black", alpha = 0.7) +
  coord_cartesian(xlim = c(-2,2)) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = c(0.2, 0.85),
        legend.text = element_text(size = 20, color = "black"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=20, margin = margin(0,5,0,0), color = "black"), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=40, margin = margin(0,10,0,0), color = "black"),
        plot.margin = margin(6,6,6,6)
  )

binned_A.m.geno
binned_A.f.geno
binned_C.m.geno <- ggplot(C.m.geno_ASE[!is.na(C.m.geno_ASE$SBGE_simp),], aes(SBGE_simp, exp_geno)) +
  geom_point(aes(SBGE_simp, exp_geno, color = SBGE_simp),size = 1, shape = 16, alpha = 0.3, show.legend = FALSE) +
  geom_boxplot(outlier.shape = NA, fill= c("red3", "#888888", "steelblue3"), alpha = 0.7) + 
  labs( # "omegaA_MK" = expression(italic("\u03c9A")[MK]); "alpha_MK" = expression(italic("\u03b1")[MK])
    y = "Difference in expression (log2FC Red/NonRed)") + 
  scale_color_manual(values=c("red3", "#888888", "steelblue3")) +
  guides(color = guide_legend(override.aes = list(shape = c(16, 16, 16),
                                                  size = c(4, 4, 4),
                                                  alpha = 1))) +
  # geom_text(data = permed_t_tests_RvNR %>% mutate(sig1 = if_else(significance == 1, "*", "")),
  #           aes(x =SBGE_comp, y = 0.95, label = sig1), size = 9, position = position_dodge(width = 0.8), vjust = -4) +
  geom_abline(intercept = 0, slope = 0,  size = 0.5, linetype="dashed", color = "black") +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, by = 0.5)) +
  coord_cartesian(xlim = c(-2,2)) +
  coord_flip() + 
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size=20, margin = margin(5,0,0,0), color = "black"),
        axis.title.x = element_text(size=40, margin = margin(10,0,0,0), color = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin = margin(6,6,6,6)
  ) 

hist_bin_A.m <- plot_grid(hist_A.m.geno, binned_A.m.geno,
                          ncol=1, align="v", rel_heights=c(4,1), axis = 'lrbt')

hist_bin_A.f <- plot_grid(hist_A.f.geno, binned_A.f.geno,
                          ncol=1, align="v", rel_heights=c(4,1), axis = 'lrbt')

hist_bin_C.m <- plot_grid(hist_C.m.geno, binned_C.m.geno,
                          ncol=1, align="v", rel_heights=c(4,1), axis = 'lrbt')

pdf(file = "~/Desktop/UofT/SSAV_RNA/Plots/HollisLike_graph.pdf",  # The directory you want to save the file in
    width = 28, # The width of the plot in inches
    height = 24) # The height of the plot in inches
ggarrange(hist_bin_A.f, NA,  hist_bin_A.m, 
          NA,NA,NA,
          NA, NA, hist_bin_C.m,    
          labels = c("A)", NA, "B)", NA, NA, NA, NA, NA, "C)"),
          heights = c(1, 0.05, 1), widths = c(1, 0.05, 1),
          ncol =3, nrow = 3, 
          font.label = list(size = 30)) 
dev.off()
########



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


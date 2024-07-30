###################################
#
#                             Grieshop et al. 2024
#                             Author: Michelle Liu
#             DsRed experimental evolution - transcriptomics analysis
#             Direction of Red/NR regulatory changes by SBGE catergory
#                         Figure 4, Figure S7, Table S5
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
library(cowplot)
#########

# Get Mishra et al. 2022's data 
##########

ASE <- read.delim(file="~/Desktop/UofT/SSAV_RNA/Data/DifferentialGeneExpression.whole.bodies.tsv", sep="\t", header=TRUE)
# Grab the desired variables from there
ASE <- data.frame(cbind(ASE$log2FoldChange,
                        ASE$lfcSE,
                        ASE$FlyBaseID))
colnames(ASE) <- c("exp_SBGE_ase", "se_SBGE_ase", "FlyBaseID")
# fix formatting
ASE$exp_SBGE_ase <- as.numeric(ASE$exp_SBGE_ase)
ASE$se_SBGE_ase <- as.numeric(ASE$se_SBGE_ase)
str(ASE)


# Define three levels of SBGE categorization 
x1 = 1 # first cut-off (FBG < -1, MBG > 1 , -1 < UBG < 1)
x2 = 5 # second cut-off (extreme FBG < -5, extreme MBG > 5)
y0 = 1 # tolerance of middle bins ### I DONT GET THIS
# xmid1 = (x1 + x2)/2

# Simple (3 levels)
# one level of female-biased gene expression
fbg.keep <- ASE$exp_SBGE_ase < -x1 & (ASE$exp_SBGE_ase + ASE$se_SBGE_ase) < 0
fbg <- ASE[fbg.keep,]
fbg$SBGE_simp <- rep(c("a.fbg"), dim(fbg)[1])
# one unbiased category
ubg.keep <- ASE$exp_SBGE_ase < x1 & ASE$exp_SBGE_ase > -x1 & (ASE$exp_SBGE_ase - ASE$se_SBGE_ase) > -(x1+y0) & (ASE$exp_SBGE_ase + ASE$se_SBGE_ase) < (x1+y0)
ubg <- ASE[ubg.keep,]
ubg$SBGE_simp <- rep(c("b.ubg"), dim(ubg)[1])
# two levels of male-biased gene expression
mbg.keep <- ASE$exp_SBGE_ase > x1 & (ASE$exp_SBGE_ase - ASE$se_SBGE_ase) > 0
mbg <- ASE[mbg.keep,]
mbg$SBGE_simp <- rep(c("c.mbg"), dim(mbg)[1])
# 1 gene is tossed out b/c the uncertainty in its estimate breaches a cutoff boundary 
ASE <- rbind(fbg, mbg, ubg)
str(ASE)

# Complex (5 levels)
more.fbg.keep <- ASE$exp_SBGE_ase < -x2 & (ASE$exp_SBGE_ase + ASE$se_SBGE_ase) < -(x2-y0) # extreme FBG < -5
more.fbg <- ASE[more.fbg.keep,]
more.fbg$SBGE_comp <- rep(c("a.more.fbg"), dim(more.fbg)[1])
#
fbg.keep <- ASE$exp_SBGE_ase < -x1 & ASE$exp_SBGE_ase > -x2 & (ASE$exp_SBGE_ase + ASE$se_SBGE_ase) < 0 & (ASE$exp_SBGE_ase - ASE$se_SBGE_ase) > -(x2+y0)
fbg <- ASE[fbg.keep,]
fbg$SBGE_comp <- rep(c("b.fbg"), dim(fbg)[1])
# one unbiased category
ubg.keep <- ASE$exp_SBGE_ase < x1 & ASE$exp_SBGE_ase > -x1 & (ASE$exp_SBGE_ase - ASE$se_SBGE_ase) > -(x1+y0) & (ASE$exp_SBGE_ase + ASE$se_SBGE_ase) < (x1+y0)
ubg <- ASE[ubg.keep,]
ubg$SBGE_comp <- rep(c("c.ubg"), dim(ubg)[1])
# two levels of male-biased gene expression
mbg.keep <- ASE$exp_SBGE_ase > x1 & ASE$exp_SBGE_ase < x2 & (ASE$exp_SBGE_ase - ASE$se_SBGE_ase) > 0 & (ASE$exp_SBGE_ase + ASE$se_SBGE_ase) < (x2+y0)
mbg <- ASE[mbg.keep,]
mbg$SBGE_comp <- rep(c("d.mbg"), dim(mbg)[1])
#
more.mbg.keep <- ASE$exp_SBGE_ase > x2 & (ASE$exp_SBGE_ase - ASE$se_SBGE_ase) > (x2-y0)
more.mbg <- ASE[more.mbg.keep,]
more.mbg$SBGE_comp <- rep(c("e.more.mbg"), dim(more.mbg)[1])
# 3 genes tossed out  b/c the uncertainty in their estimate breaches a cutoff boundary
ASE <- rbind(more.fbg, fbg, ubg, mbg, more.mbg)
str(ASE)

##########

# Get chromosome locations  
##########

all.genes <- read.delim(file="~/Desktop/UofT/SSAV_RNA/Data/all.genes.tsv", sep="\t", header=FALSE)
colnames(all.genes) = c("FlyBaseID")

Xchr <- read.delim(file="~/Desktop/UofT/SSAV_RNA/Data/X.chromosome.genes.tsv", sep="\t", header=TRUE)
colnames(Xchr) = c("FlyBaseID")
Xchr$Chr <- rep("X", dim(Xchr)[1])

Ychr <- read.delim(file="~/Desktop/UofT/SSAV_RNA/Data/Y.chromosome.genes.tsv", sep="\t", header=TRUE)
colnames(Ychr) = c("FlyBaseID")
Ychr$Chr <- rep("Y", dim(Ychr)[1])


chr2L <- read.delim(file="~/Desktop/UofT/SSAV_RNA/Data/2L.chromosome.genes.tsv", sep="\t", header=FALSE)
colnames(chr2L) = c("FlyBaseID")
chr2L$Chr <- rep("2L", dim(chr2L)[1])
chr2R <- read.delim(file="~/Desktop/UofT/SSAV_RNA/Data/2R.chromosome.genes.tsv", sep="\t", header=FALSE)
colnames(chr2R) = c("FlyBaseID")
chr2R$Chr <- rep("2R", dim(chr2R)[1])
# 
chr2 <- rbind(chr2L, chr2R)
chr2$Chr <- rep("2", dim(chr2)[1])

chr3L <- read.delim(file="~/Desktop/UofT/SSAV_RNA/Data/3L.chromosome.genes.tsv", sep="\t", header=FALSE)
colnames(chr3L) = c("FlyBaseID")
chr3L$Chr <- rep("3L", dim(chr3L)[1])
chr3R <- read.delim(file="~/Desktop/UofT/SSAV_RNA/Data/3R.chromosome.genes.tsv", sep="\t", header=FALSE)
colnames(chr3R) = c("FlyBaseID")
chr3R$Chr <- rep("3R", dim(chr3R)[1])
#
chr3 <- rbind(chr3L, chr3R)
chr3$Chr <- rep("3", dim(chr3)[1])

Chrs <- rbind(Xchr, Ychr, chr2, chr3) # Not all genes; just X, Y, 2, and 3.
Chrs$Chr <- as.factor(Chrs$Chr)

Chrs_All <- rbind(Xchr, Ychr, chr2L, chr2R, chr3L, chr3R)
Chrs_All$Chr <- as.factor(Chrs_All$Chr)
##########


# Prepare plotting dataset containing all samples
#########
# load datasets
A.f.geno <- read.delim("Results/A.f.geno_candidates.tsv")
A.m.geno <- read.delim("Results/A.m.geno_candidates.tsv")
C.m.geno <- read.delim("Results/C.m.geno_candidates.tsv")

# merge dataset with out population log2FC M/F estimates from Mishra et al. 2022
A.m.geno <- merge(A.m.geno, ASE, by = "FlyBaseID", all = T)
A.f.geno <- merge(A.f.geno, ASE, by = "FlyBaseID", all = T)
C.m.geno <- merge(C.m.geno, ASE, by = "FlyBaseID", all = T)
A.m.geno$trt2 = "Am"
A.f.geno$trt2 = "Af"
C.m.geno$trt2 = "Cm"

# only keep genes with data available in both datasets
A.m.geno <- A.m.geno[!is.na(A.m.geno$exp_geno) &
                       !is.na(A.m.geno$exp_SBGE_ase),]
A.f.geno <- A.f.geno[!is.na(A.f.geno$exp_geno) &
                       !is.na(A.f.geno$exp_SBGE_ase),]
C.m.geno <- C.m.geno[!is.na(C.m.geno$exp_geno) &
                       !is.na(C.m.geno$exp_SBGE_ase),]
dim(C.m.geno) # check how many genes are cut off due to not having SBGE info in Mishra et al. 2022

# merge all to one dataframe object
All.geno <- rbind(A.m.geno[-5], A.f.geno, C.m.geno)

# set as factors
All.geno$SBGE_simp <- as.factor(All.geno$SBGE_simp)
All.geno$SBGE_comp <- as.factor(All.geno$SBGE_comp)
All.geno$trt2 <- as.factor(All.geno$trt2)

# combine info for chromosome locations
# All.geno <- merge(All.geno, Chrs, by = "FlyBaseID")
str(All.geno)
#########


# Table S5. 
# Using all genes, permute log2FC Red/NR against 0 for each SBGE category
########
## One sample permutation test for differences n.e to 0 (see boot_permute.R for function)
## generates Table S5
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

# save table (Table S5)
# write.table(permed_All.geno, file = "~/Desktop/UofT/SSAV_RNA/Results/permed_All.geno.tsv", sep = "\t", # Fix file name accordingly
#             row.names = FALSE, col.names = TRUE)


# Additionally, permute the difference between log2FC Red/NR in Experimental males and in Control males
perm_A.m_C.m <- TwoPerm_SBGE(perm_dat = All.geno[All.geno$trt2 != "Af",] %>%
                               mutate(Am = ifelse(trt2 == "Am", TRUE, FALSE)), 
                             x_col = "exp_geno", groupBy = "Am", SBGE_cat = "SBGE_comp")
########


# plotting function for Fig.4 & S7
binPlot_RedNR <- function(dat, perm_dat){
  
  # perm_dat contains the p-values from the one group permutation test (against 0)
  perm_dat$pval <- as.numeric(perm_dat$pval)
  
  ggplot(dat, aes(SBGE_comp, exp_geno, color = trt2)) +
    geom_point(aes(color = trt2), size = 1, shape = 16, 
               alpha = 0.3, position = position_jitterdodge(jitter.width = 0.35)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8) + 
    labs(x = "Sex-biased Gene Expression", # "omegaA_MK" = expression(italic("\u03c9A")[MK]); "alpha_MK" = expression(italic("\u03b1")[MK])
         y = expression(Log["2"]*"FC (Red/NonRed)")) +
    # title = "C-males") +
    scale_colour_manual(values = c("#D55E00", "#0072B2", "#888888"), # "red3", "steelblue3", "#888888" # "purple3", "chartreuse3", "orange2"
                        labels = c("Exp. Females", "Exp. Males", "Ctrl. Males")) + # "Chr-2", "Chr-3", "X-Chr"
    guides(color = guide_legend(override.aes = list(shape = c(18, 18, 18),
                                                    size = c(5, 5, 5),
                                                    alpha = 1))) +
    
    # add star do signify significant difference from 0
    geom_text(data = perm_dat %>% mutate(sig1 = if_else(n < 30, "",
                                                        ifelse(pval < 0.001, "***", 
                                                               ifelse(pval < 0.01, "**" ,
                                                                      ifelse(pval < 0.05, "*", "ns"))))),
              aes(x =SBGE_comp, y = 1.6, label = sig1), size = 7.5,
              position = position_dodge(width = 0.8), show.legend = FALSE) +
    
    # add number of genes per category
    geom_text(data = perm_dat, aes(label = n, y = Inf, group = trt2), 
              position = position_dodge(width = 0.8), vjust = 5, size = 6, show.legend = FALSE) +
    
    # line at y = 0
    geom_abline(intercept = 0, slope = 0,  size = 0.5, linetype= "solid", color = "black") +
    
    # segment lines between SBGE categories
    geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5), color = "grey") +
    scale_x_discrete(labels = c("Highly FB", "Female-Biased", "Unbiased", "Male-Biased", "Highly MB")) + # "Highly ant.", "Antagonistic", "Uncorrelated", "Concordant", "Highly con." .... "Strong pur.", "Purifying sel.", "Neutral", "Positive sel.", "Strong pos."
    # scale_y_continuous(limits = c(-1.6, 1.6), breaks = c(-1.5, -1.0, 0, 1.0, 1.5)) +
    
    # default theme settings:
    theme_classic() +
    theme(plot.title.position = c("panel"),
          legend.title = element_blank(),
          legend.position = c("bottom"),
          legend.text = element_text(size = 25, color = "black"),
          axis.text.x = element_text(size=20, margin = margin(5,0,0,0), color = "black"),
          axis.text.y = element_text(size=20, margin = margin(0,5,0,0), color = "black"),
          axis.title.x = element_text(size=30, margin = margin(10,0,0,0), color = "black"),
          axis.title.y = element_text(size=30, margin = margin(0,20,0,0), color = "black"),
          plot.title = element_text(size=40, margin = margin(0,0,0,0), color = "black"),
          plot.margin = margin(6,6,6,6),
          panel.border = element_rect(colour = "black", fill=NA, size=1)
    ) 
}


Fig_4 <- binPlot_RedNR(All.geno, permed_All.geno) + 
  # add permutation result from comparison b/t Exp. and Ctrl. males
  geom_signif(y_position = c(-1.2, -1.2, -1.3, -1.4, -1.35), xmin = c(1, 2, 3, 4, 5), 
              xmax = c(1.3, 2.3, 3.3, 4.3, 5.3),
              annotation = c("***", "***", "***", "***", "***"), tip_length = -0.01, 
              textsize = 10, size = 0.75, vjust = 1.85, color = "darkblue")

# Comment in to save plot
# pdf(file = "~/Desktop/UofT/SSAV_RNA/Plots/final_2/Fig4_main.pdf",  # The directory you want to save the file in
#     width = 15, # 15 The width of the plot in inches
#     height = 8 ) # 8 20 The height of the plot in inches
# Fig_4
# dev.off()



# Repeat the above analysis but only for DE genes
### first, permute against 0
######
#------- all DE genes, regardless of the sex it is DE in
all.DE.genes <- unique(c(A.m.geno$FlyBaseID[A.m.geno$Sig], A.f.geno$FlyBaseID[A.f.geno$Sig]))
permed_A.f.geno_Sig <- OnePerm_SBGE(perm_dat = All.geno[All.geno$trt2 == "Af" &
                                                          All.geno$FlyBaseID %in% all.DE.genes,],
                                    x_col = "exp_geno", SBGE_cat = "SBGE_comp")
permed_A.m.geno_Sig <- OnePerm_SBGE(perm_dat = All.geno[All.geno$trt2 == "Am" &
                                                          All.geno$FlyBaseID %in% all.DE.genes,],
                                    x_col = "exp_geno", SBGE_cat = "SBGE_comp")
permed_C.m.geno_Sig <- OnePerm_SBGE(perm_dat = All.geno[All.geno$trt2 == "Cm" &
                                                          All.geno$FlyBaseID %in% all.DE.genes,],
                                    x_col = "exp_geno", SBGE_cat = "SBGE_comp")
permed_A.f.geno_Sig$trt2 <- "Af"
permed_A.m.geno_Sig$trt2 <- "Am"
permed_C.m.geno_Sig$trt2 <- "Cm"
permed_All.geno_Sig <- rbind(permed_A.f.geno_Sig, permed_A.m.geno_Sig, permed_C.m.geno_Sig)
permed_All.geno_Sig <- permed_All.geno_Sig %>%
  mutate(holm_padj = p.adjust(permed_All.geno_Sig$pval, method = "holm")) %>% 
  mutate(holm_Sig = ifelse(holm_padj < 0.005 & n > 30, TRUE, ifelse(n < 30, NA, FALSE)))
rm(permed_A.f.geno_Sig, permed_A.m.geno_Sig, permed_C.m.geno_Sig) # remove clutter



#------- only genes DE in females
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
  mutate(holm_Sig = ifelse((n > 30 & holm_padj < 0.05), TRUE, ifelse(n < 30, NA, FALSE)))
rm(permed_A.f.geno_AfSig, permed_A.m.geno_AfSig, permed_C.m.geno_AfSig) # remove clutter



#------- only candidate genes DE in males
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
  mutate(holm_Sig = ifelse((n > 30 & holm_padj < 0.05), TRUE, ifelse(n < 30, NA, FALSE)))
rm(permed_A.f.geno_AmSig, permed_A.m.geno_AmSig, permed_C.m.geno_AmSig) # remove clutter

######

### for each analysis using the subset of genes, permute difference between Exp. and Ctrl. males
######
#------- all DE genes, regardless of the sex it is DE in
perm_A.m_C.m_Sig <- TwoPerm_SBGE(perm_dat = All.geno[All.geno$trt2 != "Af" &
                                                       All.geno$FlyBaseID %in% all.DE.genes,] %>%
                                   mutate(Am = ifelse(trt2 == "Am", TRUE, FALSE)), 
                                 x_col = "exp_geno", groupBy = "Am", SBGE_cat = "SBGE_comp")

#------- only candidate genes DE in females
perm_A.m.C.m_FemSig <- TwoPerm_SBGE(perm_dat = All.geno[All.geno$trt2 != "Af" &
                                                          All.geno$FlyBaseID %in% A.f.geno[A.f.geno$Sig,]$FlyBaseID,] %>%
                                      mutate(Am = ifelse(trt2 == "Am", TRUE, FALSE)), 
                                    x_col = "exp_geno", groupBy = "Am", SBGE_cat = "SBGE_comp")

#------- only candidate genes DE in males
perm_A.m.C.m_MaleSig <- TwoPerm_SBGE(perm_dat = All.geno[All.geno$trt2 != "Af" &
                                                           All.geno$FlyBaseID %in% A.m.geno[A.m.geno$Sig,]$FlyBaseID,] %>%
                                       mutate(Am = ifelse(trt2 == "Am", TRUE, FALSE)), 
                                     x_col = "exp_geno", groupBy = "Am", SBGE_cat = "SBGE_comp")

######



### generate each plot
#------- all DE genes, regardless of the sex it is DE in
FigS7_A <- binPlot_RedNR(All.geno[All.geno$FlyBaseID %in% SSAV.geno[SSAV.geno$Sig,]$FlyBaseID,], 
                                  permed_All.geno_Sig) +
  # add permutation result from comparison b/t Exp. and Ctrl. males
  geom_signif(y_position = c(-1.18, -1.4, -1.3, -1.3), xmin = c(2, 3, 4, 5), 
              xmax = c(2.3, 3.3, 4.3, 5.3),
              annotation = c("***", "***", "***", "***"), tip_length = -0.01, 
              textsize = 10, size = 0.75, vjust = 1.8, color = "darkblue")

#------- only candidate genes DE in males
FigS7_B <- binPlot_RedNR(All.geno[All.geno$FlyBaseID %in% A.m.geno[A.m.geno$Sig,]$FlyBaseID,],
                                  permed_All.geno_AmSig) +
  # add permutation result from comparison b/t Exp. and Ctrl. males
  geom_signif(y_position = c(-1.27, -1.255, -1.22), xmin = c(3, 4, 5), 
              xmax = c(3.3, 4.3, 5.3),
              annotation = c("***", "***", "***"), tip_length = -0.01, 
              textsize = 10, size = 0.75, vjust = 1.8, color = "darkblue")

#------- only candidate genes DE in females
FigS7_C <- binPlot_RedNR(All.geno[All.geno$FlyBaseID %in% A.f.geno[A.f.geno$Sig,]$FlyBaseID,],
                                  permed_All.geno_AfSig) +
  # add permutation result from comparison b/t Exp. and Ctrl. males
  geom_signif(y_position = c(-1.1, -1.1, -1.25), xmin = c(2, 3, 4), 
              xmax = c(2.3, 3.3, 4.3),
              annotation = c("***", "***", "***"), tip_length = -0.01, 
              textsize = 10, size = 0.75, vjust = 1.8, color = "darkblue")

# comment in to save plot
# pdf(file = "~/Desktop/UofT/SSAV_RNA/Plots/final_2/Fig_S7.pdf",  # The directory you want to save the file in
#     width = 15, height = 20 )
# ggarrange(FigS7_A + theme(axis.title.x = element_blank(), axis.text.x = element_blank()),
#           NA,
#           FigS7_B + theme(axis.title.x = element_blank(), axis.text.x = element_blank()),
#           NA, 
#           FigS7_C, 
#           NA,
#           labels = c("A)", NA, "B)", NA, "C)", NA),
#           heights = c(1, 0.05, 1, 0.05, 1, 0.01), ncol =1, nrow = 6,
#           font.label = list(size = 25),
#           common.legend = TRUE, legend = "bottom")
# dev.off()

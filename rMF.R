###################################
#
#                             Grieshop et al. 2023
#                             Author: Michelle Liu
#             DsRed experimental evolution - transcriptomics analysis
#                Comparing rMF between candidates vs non-candidates
# 
# 
###################################

# load results if not loaded in env.
A.f.geno <- read.delim("Results/A.f.geno_candidates.tsv")
A.m.geno <- read.delim("Results/A.m.geno_candidates.tsv")

# load rMF calculations by Aneil (also in Singh & Agrawal 2023)
Aneil <- read.delim(file="~/Desktop/UofT/SSAV_RNA/Data/rmf.For.Karl.To.Use.v2.txt", header=TRUE)
colnames(Aneil) <- c("FlyBaseID", "rMF", "Var.F", "pVal.F", "Var.M", "pVal.M")
Aneil <- na.omit(Aneil) # There's an N/A for rMF

# Combine results into one data frame
##########
# Genes present in both SSAV males and SSAV females data
SSAV.geno <- merge(A.m.geno, A.f.geno, by = "FlyBaseID", all = TRUE)
colnames(SSAV.geno) <- c("FlyBaseID", "A.m.exp_geno", "A.m.se_geno", "A.m.padj", "A.m.TopSig", "A.m.Sig",
                                 "A.f.exp_geno", "A.f.se_geno", "A.f.padj", "A.f.Sig")
SSAV.geno <- SSAV.geno %>% mutate(Sig = A.m.Sig | A.f.Sig) # column denotes genes that are candidates in males or females
SSAV.geno <- merge(SSAV.geno, Aneil, by = "FlyBaseID", all = TRUE)
SSAV.geno <- SSAV.geno[!is.na(SSAV.geno$Sig) & !is.na(SSAV.geno$rMF),]
##########

## plotting function
pointSEplot <- function(boot_dat, perm_dat, x_col, SBGE_cat = NA){
  # set y-axis value above each error bar
  if(!is.na(SBGE_cat)){
    y_count_FALSE <- boot_dat[!boot_dat$Sig,] %>% 
      group_by(SBGE_comp) %>%
      summarise(max = q95 + 0.05)
    y_count_TRUE <- boot_dat[boot_dat$Sig,] %>%
      group_by(SBGE_comp) %>%
      summarise(max = q95 + 0.05)
  }
  else{
    y_count_FALSE <- boot_dat[!boot_dat$Sig,] %>% 
      summarise(max = q95 + 0.05)
    y_count_TRUE <- boot_dat[boot_dat$Sig,] %>% 
      summarise(max = q95 + 0.05)
    perm_dat <- as.data.frame(t(perm_dat))
    colnames(perm_dat) <- c("obs_diff", "n_TRUE", "n_FALSE", "pval", "Sig")
    str(perm_dat)
  }
  # join y-axis coordinate with text dataframes
  perm_dat$y_count_FALSE <- y_count_FALSE$max
  perm_dat$y_count_TRUE <- y_count_TRUE$max
  perm_dat$y_star <- perm_dat$y_count_FALSE + 0.1
  
  if(!is.na(SBGE_cat)){
    pointSE_plot <- ggplot(boot_dat) + 
      geom_point(aes_string(SBGE_cat, "q50", color = "Sig"), 
                 position = position_jitterdodge(jitter.width = 0.001), size = 6) +
      geom_errorbar(aes_string(SBGE_cat, "q50", color = "Sig"),
                    ymin = boot_dat$q05, ymax = boot_dat$q95, 
                    position = position_jitterdodge(jitter.width = 0.01), width = 0.5) +
      labs(x = "SBGE (ASE)", # "omegaA_MK" = expression(italic("\u03c9A")[MK]); "alpha_MK" = expression(italic("\u03b1")[MK])
           y = print(x_col)) +
      geom_text(data = perm_dat %>% mutate(sig1 = if_else(Sig==1, "*", "")),
                aes(x =SBGE_comp, y = y_star, label = sig1), color = "black", size = 10) +
      geom_text(data = perm_dat,aes(x =SBGE_comp, y = y_count_FALSE, label = n_FALSE), color = "black", size = 6,
                hjust = 1) +
      geom_text(data = perm_dat,aes(x =SBGE_comp, y = y_count_TRUE, label = n_TRUE), color = "black", 
                size = 5.5, hjust = -1) +
      scale_x_discrete(labels = c("Highly FB","Female-Biased", "Unbiased", "Male-Biased", "Highly MB")) 
  }
  else{
    pointSE_plot <- ggplot(boot_dat) +
      geom_point(aes(x = Sig, y = q50, color = Sig),
                 position = position_jitterdodge(jitter.width = 0.001), size = 6) +
      geom_errorbar(aes(x = Sig, y = q50, color = Sig),
                    ymin = boot_dat$q05, ymax = boot_dat$q95, 
                    position = position_jitterdodge(jitter.width = 0.01), width = 0.5) +
      labs(y = print(x_col)) +
      geom_text(data = perm_dat %>% mutate(sig1 = if_else(Sig==1, "*", "")),
                aes(x = 1.5, y = y_star, label = sig1), color = "black", size = 10) +
      geom_text(data = perm_dat,aes(x = 1.04, y = y_count_FALSE, label = n_FALSE), 
                color = "black", size = 6, hjust = 1) +
      geom_text(data = perm_dat,aes(x = 1.9, y = y_count_TRUE, label = n_TRUE), color = "black", 
                size = 5.5, hjust = -1) 
  }
  
  pointSE_plot <- pointSE_plot +
    scale_colour_manual(values = c("orange3", "forestgreen"), # "red3", "steelblue3", "#888888" # "purple3", "chartreuse3", "orange2"
                        labels = c("background", "candidates")) + # "Chr-2", "Chr-3", "X-Chr"
    guides(color = guide_legend(override.aes = list(shape = c(16, 16),
                                                    size = c(4, 4),
                                                    alpha = 1))) +
    scale_y_continuous(limits = c(0,1)) +
    theme_classic() +
    theme(plot.title.position = c("panel"),
          legend.title = element_blank(),
          legend.position = c("bottom"),
          legend.box.background = element_rect(),
          legend.text = element_text(size = 20, color = "black"),
          axis.text.x = element_text(size=20, margin = margin(5,0,0,0), color = "black"),
          axis.text.y = element_text(size=20, margin = margin(0,5,0,0), color = "black"),
          axis.title.x = element_text(size=40, margin = margin(10,0,0,0), color = "black"),
          axis.title.y = element_text(size=40, margin = margin(0,10,0,0), color = "black"),
          plot.title = element_text(size=40, margin = margin(0,0,0,0), color = "black"),
          plot.margin = margin(6,6,6,6)
    )
  
  return(pointSE_plot)
}


# mean rMF candidates vs non-candidates
########
# use the boot and permute functions in boot_permute.R
perm_rMF <- TwoPerm(SSAV.geno, x_col = "rMF", groupBy = "Sig")
boot_rMF <- TwoBoot(SSAV.geno, x_col = "rMF", groupBy = "Sig")
rMF_all <- pointSEplot(boot_dat = boot_rMF, perm_dat = perm_rMF, x_col = "rMF")
########



# rMF by SBGE categories
#########
# for all candidates (combined male and female candidates)
SSAV.geno_ASE <- merge(SSAV.geno, ASE, by = "FlyBaseID", all = TRUE)
SSAV.geno_ASE <- SSAV.geno_ASE[!is.na(SSAV.geno_ASE$Sig) &
                                 !is.na(SSAV.geno_ASE$rMF) &
                                 !is.na(SSAV.geno_ASE$exp_SBGE_ase),]
perm_rMF_SBGE <- TwoPerm_SBGE(SSAV.geno_ASE, x_col = "rMF", groupBy = "Sig", SBGE_cat = "SBGE_comp")
boot_rMF_SBGE <- TwoBoot_SBGE(SSAV.geno_ASE, x_col = "rMF", groupBy = "Sig", SBGE_cat = "SBGE_comp")
rMF_SBGE <- pointSEplot(boot_dat = boot_rMF_SBGE, perm_dat = perm_rMF_SBGE, 
                        x_col = "rMF", SBGE_cat = "SBGE_comp")


# only female candidates
perm_rMF_SBGE_A.f <- TwoPerm_SBGE(SSAV.geno_ASE[!is.na(SSAV.geno_ASE$A.f.Sig),], 
                                  x_col = "rMF", groupBy = "A.f.Sig", SBGE_cat = "SBGE_comp")
boot_rMF_SBGE_A.f <- TwoBoot_SBGE(SSAV.geno_ASE[!is.na(SSAV.geno_ASE$A.f.Sig),],
                                  x_col = "rMF", groupBy = "A.f.Sig", SBGE_cat = "SBGE_comp")
rMF_SBGE_A.f <- pointSEplot(boot_dat = boot_rMF_SBGE_A.f, perm_dat = perm_rMF_SBGE_A.f, 
                        x_col = "rMF", SBGE_cat = "SBGE_comp") # might just cut off the Highly SB genes


# only male candidates
perm_rMF_SBGE_A.m <- TwoPerm_SBGE(SSAV.geno_ASE[!is.na(SSAV.geno_ASE$A.m.Sig),], 
                                  x_col = "rMF", groupBy = "A.m.Sig", SBGE_cat = "SBGE_comp")
boot_rMF_SBGE_A.m <- TwoBoot_SBGE(SSAV.geno_ASE[!is.na(SSAV.geno_ASE$A.m.Sig),],
                                  x_col = "rMF", groupBy = "A.m.Sig", SBGE_cat = "SBGE_comp")
rMF_SBGE_A.m <- pointSEplot(boot_dat = boot_rMF_SBGE_A.m, perm_dat = perm_rMF_SBGE_A.m, 
                            x_col = "rMF", SBGE_cat = "SBGE_comp") # might just cut off the Highly SB genes

#########


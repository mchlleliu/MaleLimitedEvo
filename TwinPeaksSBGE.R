###################################
#
#                             Grieshop et al. 2023
#                             Author: Michelle Liu
#             DsRed experimental evolution - transcriptomics analysis
#        Distribution of SA candidates by sex-bias - Cheng & Kirkpatrick model
# 
# 
###################################

# This uses Mishra et al.'s data, 
# Or can be used with any data which have log2FC value of SBGE
# Bins for SBGE categories should be defined first. 
#     -> See External_data.R under "Get Mishra et al.'s data" for how to prepare the external dataset.

rm(list=ls()) # Clears the environment
setwd("~/Desktop/UofT/SSAV_RNA/")

# packages required
#########
library(tidyr)
library(plyr)
library(dplyr)
library(broom)
library(ggplot2)
library(lmtest)
library(splines)
library(mgcv)
########

# Get Mishra et al.'s data 
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
rm(more.fbg, fbg, ubg, mbg, more.mbg)
##########


# Prepare plotting dataset 
# (Run this first before anything else!!, make sure you have the correct files/path to them)
########
# load results if not loaded in env.
A.f.geno <- read.delim("Results/A.f.geno_candidates.tsv")
A.m.geno <- read.delim("Results/A.m.geno_candidates.tsv")
SSAV.geno <- read.delim("Results/All.geno_candidates.tsv")

# include SBGE categories (using Mishra et al. dataset. Look at External_data.R)
A.m.geno_ASE <- merge(A.m.geno, ASE, by = "FlyBaseID", all = TRUE)
A.m.geno_ASE <- A.m.geno_ASE[!is.na(A.m.geno_ASE$Sig) & !is.na(A.m.geno_ASE$exp_SBGE_ase),]
A.m.geno_ASE$SBGE_comp <- as.factor(A.m.geno_ASE$SBGE_comp)
A.m.geno_ASE$SBGE_simp <- as.factor(A.m.geno_ASE$SBGE_simp)
str(A.m.geno_ASE)

A.f.geno_ASE <- merge(A.f.geno, ASE, by = "FlyBaseID", all = TRUE)
A.f.geno_ASE <- A.f.geno_ASE[!is.na(A.f.geno_ASE$Sig) & !is.na(A.f.geno_ASE$exp_SBGE_ase),]
A.f.geno_ASE$SBGE_comp <- as.factor(A.f.geno_ASE$SBGE_comp)
A.f.geno_ASE$SBGE_simp <- as.factor(A.f.geno_ASE$SBGE_simp)
str(A.f.geno_ASE)

# Genes present in both SSAV males and SSAV females data
SSAV.geno_ASE <- merge(SSAV.geno, ASE, by = "FlyBaseID", all = TRUE)
SSAV.geno_ASE <- SSAV.geno_ASE[!is.na(SSAV.geno_ASE$Sig) & !is.na(SSAV.geno_ASE$exp_SBGE_ase),]
SSAV.geno_ASE$SBGE_comp <- as.factor(SSAV.geno_ASE$SBGE_comp)
SSAV.geno_ASE$SBGE_simp <- as.factor(SSAV.geno_ASE$SBGE_simp)
str(SSAV.geno_ASE)

########



# Cheng & Kirkpatrick 2016 logistic regression & Twin Peaks model
# convert Log2FC(M/F) to delta
convert.Log2Sex.To.CKdelta<-function(x) (1-2^(-x))/(1+2^(-x))


hist(SSAV.geno_ASE$exp_SBGE_ase, breaks = 100)
hist(A.f.geno_ASE$exp_SBGE_ase, breaks = 100)
hist(A.m.geno_ASE$exp_SBGE_ase, breaks = 100)

# splines with GAM
#######
plotSplinesDelta <- function(df){
  gam_mod <- gam(I(Sig) ~ s(sapply(exp_SBGE_ase, FUN = convert.Log2Sex.To.CKdelta), k = 5, bs = "cs"),
               data = df, method = "REML", family = "binomial")
  print(summary(gam_mod))
  print(AIC(gam_mod))
  print(gam.check(gam_mod))

  gam.plot <- plot(gam_mod, residuals = TRUE)
  gam.plot <- gam.plot[[1]] # just one smooth so select the first component
  sm_df <- as.data.frame(gam.plot[c("x", "se", "fit")])

  ## plot
  return.plot <- ggplot(sm_df, aes(x = x, y = fit)) +
    geom_ribbon(aes(ymin = fit - se, ymax = fit + se, y = NULL),
              alpha = 0.3) +
    geom_line()  +
    labs(y = "Log Likelihood of Antagonism", x = "") +
    theme_classic() +
    theme(plot.title.position = c("panel"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size=18, margin = margin(0,5,0,0), color = "black"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=24, margin = margin(0,10,0,0), color = "black")) 

return(return.plot)
}
#######


# bootstrap regression (using CK delta)
boot_reg_delta <- function(dat, boot_N, degree){

  # convert to Cheng and Kirpatrick delta SBGE estimate
  dat <- dat %>% 
    dplyr::mutate(delta = convert.Log2Sex.To.CKdelta(exp_SBGE_ase))
  
  boot_pred <- NULL
  
  
  # bootstrapping
  for (i in 1:boot_N) {
    # Resample data 
    sub_sample = dat %>% 
      dplyr::sample_frac(size = 1, replace = T)
    #Running the regression on these data
    mod_boot <- glm(Sig ~ poly(delta, degree = degree, raw = TRUE), 
                    data = sub_sample, family = "binomial")
    
    # predict using current model
    boot_pred <- rbind(boot_pred, 
                      cbind(i, dat$delta, predict(mod_boot, newdata = dat, type = "response")))
  }
  
  colnames(boot_pred) <- c("bs", "delta", "Prob")
  
  out_boot <- as_tibble(boot_pred) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(delta) %>%
    dplyr::summarise(q05 = quantile(Prob, 0.025, na.rm = T),
              q50 = quantile(Prob, 0.5, na.rm = T),
              q95 = quantile(Prob, 0.975, na.rm =T))
  
  return(out_boot)
}

# bootstrap regression (using CK delta)
boot_reg_Log2SB <- function(dat, boot_N, degree){
  
  boot_pred <- NULL
  
  # bootstrapping
  for (i in 1:boot_N) {
    # Resample data 
    sub_sample = dat %>% 
      dplyr::sample_frac(size = 1, replace = T)
    #Running the regression on these data
    mod_boot <- glm(Sig ~ poly(exp_SBGE_ase, degree = degree, raw = TRUE), 
                    data = sub_sample, family = "binomial")
    
    # predict using current model
    boot_pred <- rbind(boot_pred, 
                       cbind(i, dat$exp_SBGE_ase, predict(mod_boot, newdata = dat, type = "response")))
  }
  
  colnames(boot_pred) <- c("bs", "exp_SBGE_ase", "Prob")
  
  out_boot <- as_tibble(boot_pred) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(exp_SBGE_ase) %>%
    dplyr::summarise(q05 = quantile(Prob, 0.025),
              q50 = quantile(Prob, 0.5),
              q95 = quantile(Prob, 0.975))
  
  return(out_boot)
}



A.f.geno.reg_Log2SB.noSL <- boot_reg_Log2SB(A.f.geno_ASE.noSL, 1000, 4)
A.f.geno.reg_delta.noSL <- boot_reg_delta(A.f.geno_ASE.noSL, 1000, 4)

A.f.geno.reg_Log2SB <- boot_reg_Log2SB(A.f.geno_ASE, 1000, 4)
A.f.geno.reg_delta <- boot_reg_delta(A.f.geno_ASE, 1000, 4)

# plotting female CK plot with delta
######
# SA_dist <- ggplot(A.f.reg_delta) + 
#   geom_ribbon(aes(x= delta, ymin=q05, ymax=q95),
#               size = 0.25, alpha = 0.7, fill = "grey")  +
#   geom_line(aes(x=delta, y=q50), size = 1) +
#   labs(y="P(Sexually Antagonistic)") +
#   theme_classic() +
#   theme(plot.title.position = c("panel"),
#         axis.text.x = element_blank(),
#         axis.text.y = element_text(size=18, margin = margin(0,5,0,0), color = "black"),
#         axis.title.x = element_blank(),
#         axis.title.y = element_text(size=24, margin = margin(0,10,0,0), color = "black")) 

SA_dist <- plotSplinesDelta(A.f.geno_ASE) 

SBGE_dist <- ggplot(A.f.geno_ASE) + 
  geom_density(aes(sapply(exp_SBGE_ase, convert.Log2Sex.To.CKdelta)), fill = "grey", color = "grey") +
  labs(x=expression(Delta),y="Density") +
  theme_classic() +
  theme(axis.text.x = element_text(size=18, margin = margin(5,0,5,0), color = "black"),
        axis.text.y = element_text(size=12, margin = margin(0,5,0,0), color = "black"),
        axis.title.x = element_text(size=30, margin = margin(10,0,10,0), color = "black"),
        axis.title.y = element_text(size=15, margin = margin(0,25,0,10), color = "black")) 

datSBGE <- ggplot_build(SBGE_dist)
x1 <- min(which(datSBGE$data[[1]]$x >= -0.75))
x2 <- max(which(datSBGE$data[[1]]$x <= -0.25))
x3 <- min(which(datSBGE$data[[1]]$x >= .25))
x4 <- max(which(datSBGE$data[[1]]$x <= .75))

SBGE_dist <- SBGE_dist +  
  geom_area(data=data.frame(x=datSBGE$data[[1]]$x[x1:x2],
                                                    y=datSBGE$data[[1]]$y[x1:x2]),
                                    aes(x=x, y=y), fill="grey36") +
  geom_area(data=data.frame(x=datSBGE$data[[1]]$x[x3:x4],
                            y=datSBGE$data[[1]]$y[x3:x4]),
            aes(x=x, y=y), fill="grey36")


Fig2_CKsuppl_fem <- ggarrange(SA_dist + scale_x_continuous(expand = c(0.0005, 0.0005)), NA,
                              SBGE_dist + scale_x_continuous(expand = c(0.0005, 0.0005)), NA,
                              nrow = 2, heights = c(1, 0.30),
                              ncol = 2, widths = c(1, 0.01)) 
######

# plot female data with log2FC sex-bias
######
SA_dist <- ggplot(A.f.reg_Log2SB) + 
  geom_ribbon(aes(x= exp_SBGE_ase, ymin=q05, ymax=q95),
              size = 0.25, alpha = 0.7, fill = "grey")  +
  geom_line(aes(x=exp_SBGE_ase, y=q50), size = 1) +
  labs(y="P(Sexually Antagonistic)") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=18, margin = margin(0,5,0,0), color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=24, margin = margin(0,10,0,0), color = "black")) 


SBGE_dist <- ggplot(check_test) + 
  geom_density(aes(exp_SBGE_ase), fill = "grey", color = "grey") +
  labs(x=expression("Log2FC sex-bias"),y="Number of loci") +
  theme_classic() +
  theme(axis.text.x = element_text(size=18, margin = margin(5,0,5,0), color = "black"),
        axis.text.y = element_text(size=12, margin = margin(0,5,0,0), color = "black"),
        axis.title.x = element_text(size=24, margin = margin(10,0,10,0), color = "black"),
        axis.title.y = element_text(size=15, margin = margin(0,10,0,0), color = "black"),
        plot.margin = margin(6,6,6,6)) 

datSBGE <- ggplot_build(SBGE_dist)
x1 <- min(which(datSBGE$data[[1]]$x >= -5))
x2 <- max(which(datSBGE$data[[1]]$x <= -1))
x3 <- min(which(datSBGE$data[[1]]$x >= 1))
x4 <- max(which(datSBGE$data[[1]]$x <= 5))

SBGE_dist <- SBGE_dist +  
  geom_area(data=data.frame(x=datSBGE$data[[1]]$x[x1:x2],
                            y=datSBGE$data[[1]]$y[x1:x2]),
            aes(x=x, y=y), fill="grey36") +
  geom_area(data=data.frame(x=datSBGE$data[[1]]$x[x3:x4],
                            y=datSBGE$data[[1]]$y[x3:x4]),
            aes(x=x, y=y), fill="grey36")


Fig2_CKsuppl_fem_log2FC <- ggarrange(NA, SA_dist + scale_x_continuous(expand = c(0, 0)) + 
                                       coord_cartesian(xlim = c(-10,10)), NA,
                              NA, SBGE_dist + scale_x_continuous(expand = c(0, 0)) + 
                                coord_cartesian(xlim = c(-10,10)), NA,
                              nrow = 2, heights = c(1, 0.30),
                              ncol = 3, widths = c(0.005, 1, 0.005)) 
######

A.m.geno.reg_Log2SB.noSL <- boot_reg_Log2SB(A.m.geno_ASE.noSL, 1000, 4)
A.m.geno.reg_delta.noSL <- boot_reg_delta(A.m.geno_ASE.noSL, 1000, 4)

A.m.geno.reg_Log2SB <- boot_reg_Log2SB(A.m.geno_ASE, 1000, 4)
A.m.geno.reg_delta <- boot_reg_delta(A.m.geno_ASE, 1000, 4)

# plotting male CK plot delta
######
# SA_dist <- ggplot(A.m.geno.reg_delta) + 
#   geom_ribbon(aes(x= delta, ymin=q05, ymax=q95),
#               size = 0.25, alpha = 0.7, fill = "grey")  +
#   geom_line(aes(x=delta, y=q50), size = 1) +
#   labs(y="P(Sexually Antagonistic)") +
#   theme_classic() +
#   theme(plot.title.position = c("panel"),
#         axis.text.x = element_blank(),
#         axis.text.y = element_text(size=18, margin = margin(0,5,0,0), color = "black"),
#         axis.title.x = element_blank(),
#         axis.title.y = element_text(size=24, margin = margin(0,10,0,0), color = "black")) 

SA_dist <- plotSplinesDelta(A.m.geno_ASE)

SBGE_dist <- ggplot(A.m.geno_ASE) + 
  geom_density(aes(sapply(exp_SBGE_ase, convert.Log2Sex.To.CKdelta)), fill = "grey", color = "grey") +
  labs(x=expression(Delta),y="Density") +
  theme_classic() +
  theme(axis.text.x = element_text(size=18, margin = margin(5,0,5,0), color = "black"),
        axis.text.y = element_text(size=12, margin = margin(0,5,0,0), color = "black"),
        axis.title.x = element_text(size=30, margin = margin(10,0,10,0), color = "black"),
        axis.title.y = element_text(size=15, margin = margin(0,20,0,9), color = "black")) 

datSBGE <- ggplot_build(SBGE_dist)
x1 <- min(which(datSBGE$data[[1]]$x >= -0.75))
x2 <- max(which(datSBGE$data[[1]]$x <= -0.25))
x3 <- min(which(datSBGE$data[[1]]$x >= .25))
x4 <- max(which(datSBGE$data[[1]]$x <= .75))

SBGE_dist <- SBGE_dist +  
  geom_area(data=data.frame(x=datSBGE$data[[1]]$x[x1:x2],
                            y=datSBGE$data[[1]]$y[x1:x2]),
            aes(x=x, y=y), fill="grey36") +
  geom_area(data=data.frame(x=datSBGE$data[[1]]$x[x3:x4],
                            y=datSBGE$data[[1]]$y[x3:x4]),
            aes(x=x, y=y), fill="grey36")


Fig2_CKsuppl_male <- ggarrange(SA_dist + scale_x_continuous(expand = c(0.0005, 0.0005)), NA,
                              SBGE_dist + scale_x_continuous(expand = c(0.0005, 0.0005)), NA,
                              nrow = 2, heights = c(1, 0.30),
                              ncol = 2, widths = c(1, 0.01)) 
######

# plot male data with log2FC sex-bias
######
SA_dist <- ggplot(A.m.geno.reg_Log2SB) + 
  geom_ribbon(aes(x= exp_SBGE_ase, ymin=q05, ymax=q95),
              size = 0.25, alpha = 0.7, fill = "grey")  +
  geom_line(aes(x=exp_SBGE_ase, y=q50), size = 1) +
  labs(y="P(Sexually Antagonistic)") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=18, margin = margin(0,5,0,0), color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=24, margin = margin(0,10,0,0), color = "black")) 


SBGE_dist <- ggplot(A.m.geno.reg_Log2SB) + 
  geom_density(aes(exp_SBGE_ase), fill = "grey", color = "grey") +
  labs(x=expression("Log2FC sex-bias"),y="Number of loci") +
  theme_classic() +
  theme(axis.text.x = element_text(size=18, margin = margin(5,0,5,0), color = "black"),
        axis.text.y = element_text(size=12, margin = margin(0,5,0,0), color = "black"),
        axis.title.x = element_text(size=24, margin = margin(10,0,10,0), color = "black"),
        axis.title.y = element_text(size=15, margin = margin(0,30,0,10), color = "black"),
        plot.margin = margin(6,6,6,6)) 

datSBGE <- ggplot_build(SBGE_dist)
x1 <- min(which(datSBGE$data[[1]]$x >= -5))
x2 <- max(which(datSBGE$data[[1]]$x <= -1))
x3 <- min(which(datSBGE$data[[1]]$x >= 1))
x4 <- max(which(datSBGE$data[[1]]$x <= 5))

SBGE_dist <- SBGE_dist +  
  geom_area(data=data.frame(x=datSBGE$data[[1]]$x[x1:x2],
                            y=datSBGE$data[[1]]$y[x1:x2]),
            aes(x=x, y=y), fill="grey36") +
  geom_area(data=data.frame(x=datSBGE$data[[1]]$x[x3:x4],
                            y=datSBGE$data[[1]]$y[x3:x4]),
            aes(x=x, y=y), fill="grey36") 


Fig2_CKsuppl_male_log2FC <- ggarrange(NA, SA_dist + scale_x_continuous(expand = c(0, 0)) + 
                                       coord_cartesian(xlim = c(-10,15), ylim = c(0, 0.1)), NA,
                                     NA, SBGE_dist + scale_x_continuous(expand = c(0, 0)) + 
                                       coord_cartesian(xlim = c(-10,15)), NA,
                                     nrow = 2, heights = c(1, 0.30),
                                     ncol = 3, widths = c(0.005, 1, 0.005)) 
######



SSAV.geno.reg_Log2SB.noSL <- boot_reg_Log2SB(SSAV.geno_ASE.noSL, 1000, 4)
SSAV.geno.reg_delta.noSL <- boot_reg_delta(SSAV.geno_ASE.noSL, 1000, 4)

SSAV.geno.reg_Log2SB <- boot_reg_Log2SB(SSAV.geno_ASE, 10, 4)
SSAV.geno.reg_delta <- boot_reg_delta(SSAV.geno_ASE, 1000, 4)


# plotting all CK plot delta
######
# SA_dist <- ggplot(SSAV.geno.reg_delta) + 
#   geom_ribbon(aes(x= delta, ymin=q05, ymax=q95),
#               size = 0.25, alpha = 0.7, fill = "grey")  +
#   geom_line(aes(x=delta, y=q50), size = 1) +
#   labs(y="P(Sexually Antagonistic)") +
#   theme_classic() +
#   theme(plot.title.position = c("panel"),
#         axis.text.x = element_blank(),
#         axis.text.y = element_text(size=20, margin = margin(0,5,0,0), color = "black"),
#         axis.title.x = element_blank(),
#         axis.title.y = element_text(size=30, margin = margin(0,10,0,0), color = "black")) 

SA_dist <- plotSplinesDelta(SSAV.geno_ASE)

SBGE_dist <- ggplot(SSAV.geno_ASE) + 
  geom_density(aes(sapply(exp_SBGE_ase, convert.Log2Sex.To.CKdelta)), fill = "grey", color = "grey") +
  labs(x=expression(Delta),y="Density") +
  theme_classic() +
  theme(axis.text.x = element_text(size=20, margin = margin(5,0,5,0), color = "black"),
        axis.text.y = element_text(size=15, margin = margin(0,5,0,0), color = "black"),
        axis.title.x = element_text(size=40, margin = margin(10,0,10,0), color = "black"),
        axis.title.y = element_text(size=18, margin = margin(0,24,0,10), color = "black")) 

datSBGE <- ggplot_build(SBGE_dist)
x1 <- min(which(datSBGE$data[[1]]$x >= -0.75))
x2 <- max(which(datSBGE$data[[1]]$x <= -0.25))
x3 <- min(which(datSBGE$data[[1]]$x >= .25))
x4 <- max(which(datSBGE$data[[1]]$x <= .75))

SBGE_dist <- SBGE_dist +  
  geom_area(data=data.frame(x=datSBGE$data[[1]]$x[x1:x2],
                            y=datSBGE$data[[1]]$y[x1:x2]),
            aes(x=x, y=y), fill="grey36") +
  geom_area(data=data.frame(x=datSBGE$data[[1]]$x[x3:x4],
                            y=datSBGE$data[[1]]$y[x3:x4]),
            aes(x=x, y=y), fill="grey36") 


Fig2_CKsuppl_all <-  ggarrange(SA_dist + scale_x_continuous(expand = c(0.0005, 0.0005)), NA,
                              NA, NA,
                               SBGE_dist + scale_x_continuous(expand = c(0.0005, 0.0005)) +
                                coord_cartesian(xlim=c(-1,1)),
                               nrow = 3, heights = c(1, 0.005, 0.30), ncol = 2, widths = c(1, 0.05)) 

######

  
# plot all with log2FC sex-bias
######
SA_dist <- ggplot(SSAV.geno.reg_Log2SB) + 
  geom_ribbon(aes(x= exp_SBGE_ase, ymin=q05, ymax=q95),
              size = 0.25, alpha = 0.7, fill = "grey")  +
  geom_line(aes(x=exp_SBGE_ase, y=q50), size = 1) +
  labs(y="P(Sexually Antagonistic)") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=18, margin = margin(0,5,0,0), color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=24, margin = margin(0,10,0,0), color = "black")) 


SBGE_dist <- ggplot(SSAV.geno.reg_Log2SB) + 
  geom_density(aes(sapply(exp_SBGE_ase, convert.Log2Sex.To.CKdelta)), fill = "grey", color = "grey") +
  labs(x=expression("Log2FC sex-bias"),y="Number of loci") +
  geom_vline(xintercept = 0) +
  theme_classic() +
  theme(axis.text.x = element_text(size=18, margin = margin(5,0,5,0), color = "black"),
        axis.text.y = element_text(size=12, margin = margin(0,5,0,0), color = "black"),
        axis.title.x = element_text(size=24, margin = margin(10,0,10,0), color = "black"),
        axis.title.y = element_text(size=15, margin = margin(0,30,0,10), color = "black"),
        plot.margin = margin(6,6,6,6)) 

datSBGE <- ggplot_build(SBGE_dist)
x1 <- min(which(datSBGE$data[[1]]$x >= -5))
x2 <- max(which(datSBGE$data[[1]]$x <= -1))
x3 <- min(which(datSBGE$data[[1]]$x >= 1))
x4 <- max(which(datSBGE$data[[1]]$x <= 5))

SBGE_dist <- SBGE_dist +  
  geom_area(data=data.frame(x=datSBGE$data[[1]]$x[x1:x2],
                            y=datSBGE$data[[1]]$y[x1:x2]),
            aes(x=x, y=y), fill="grey36") +
  geom_area(data=data.frame(x=datSBGE$data[[1]]$x[x3:x4],
                            y=datSBGE$data[[1]]$y[x3:x4]),
            aes(x=x, y=y), fill="grey36")


Fig2_CKsuppl_all_log2FC <- ggarrange(NA, SA_dist + scale_x_continuous(expand = c(0, 0)) + 
                                        coord_cartesian(xlim = c(-15,10)), NA,
                                      NA, SBGE_dist + scale_x_continuous(expand = c(0, 0)) + 
                                        coord_cartesian(xlim = c(-15,10)), NA,
                                      nrow = 2, heights = c(1, 0.30),
                                      ncol = 3, widths = c(0.005, 1, 0.005)) 
######


pdf(file = "~/Desktop/UofT/SSAV_RNA/Plots/finals/Fig2_supplCK.pdf",   # The directory you want to save the file in
    width = 24, # 14 18 24 The width of the plot in inches
    height = 20) # 10 20 20 The height of the plot in inches

ggarrange(NA, Fig2_CKsuppl_all, NA,
          NA, NA, NA,
          NA, ggarrange(Fig2_CKsuppl_male, NA, Fig2_CKsuppl_fem,
                        labels = c("B)", NA, "C)"),
                        ncol = 3,
                        widths = c(1, 0.05, 1), font.label = list(size = 30), hjust = -0.01), NA,
          labels = c(NA, "A)", NA, NA, NA, NA, NA),
          nrow = 3, ncol = 3,
          heights = c(1.5, 0.05, 1), widths = c(0.05, 1, 0.05),
          font.label = list(size = 30), hjust = -0.01)

# ggarrange(Fig2_CKsuppl_all, Fig2_CKsuppl_all.noSL,
#           labels = c("A)", "B)"),
#           nrow =2, heights = c(1, 1),
#           font.label = list(size = 30), hjust = -0.01)


dev.off()


# testing models for males
#####
# Fit fourth-order polynomial logistic regression
SSAV.geno_ASE <- SSAV.geno_ASE[!is.na(SSAV.geno_ASE$A.m.exp_geno) &
                                 !is.na(SSAV.geno_ASE$A.f.exp_geno),]
model6 <- glm(Sig ~ poly(exp_SBGE_ase, degree = 6, raw = TRUE), data = SSAV.geno_ASE, family = "binomial")
model5 <- glm(Sig ~ poly(exp_SBGE_ase, degree = 5, raw = TRUE), data = SSAV.geno_ASE, family = "binomial")
model4 <- glm(Sig ~ poly(exp_SBGE_ase, degree = 4, raw = TRUE), data = SSAV.geno_ASE, family = "binomial")
model3 <- glm(Sig ~ poly(exp_SBGE_ase, degree = 3, raw = TRUE), data = SSAV.geno_ASE, family = "binomial")
model2 <- glm(Sig ~ poly(exp_SBGE_ase, degree = 2, raw = TRUE), data = SSAV.geno_ASE, family = "binomial")
model1 <- glm(Sig ~ poly(exp_SBGE_ase, degree = 1, raw = TRUE), data = SSAV.geno_ASE, family = "binomial")
model0 <- glm(Sig ~ 1, data = SSAV.geno_ASE, family = "binomial") # null model

lrtest(model2, model4)
AIC.vec = AIC(model0, model1, model2, model3, model4, model5, model6)
plot(x=seq(0,6,by=1), y=AIC.vec$AIC)

# requirement "0" for Twin Peaks
#Does model 4 have the best AIC?
TP.requirement0 = AIC.vec[5,2] == min(AIC.vec$AIC)

# test quadratic vs quartic models
lrtest(model2, model4)



ggplot(SSAV.geno_ASE %>% mutate(prob=ifelse(Sig, 1, 0))) + 
  stat_smooth(aes(exp_SBGE_ase, prob), formula = y ~ poly(x, degree = 4, raw = TRUE), 
              method = "glm", method.args = list(family = "binomial")) 

ggplot(SSAV.geno_ASE %>% mutate(prob=ifelse(Sig, 1, 0))) + 
  stat_smooth(aes(sapply(exp_SBGE_ase, convert.Log2Sex.To.CKdelta), prob), 
              formula = y ~ poly(x, degree = 4, raw = TRUE), 
              method = "glm", method.args = list(family = "binomial"), se = T )


######


# permute SBGE to test for twinPeaks pattern
# function to test for twinPeaks
twinPeaks <- function(dat){
  dat <- dat %>%
    dplyr::mutate(delta = convert.Log2Sex.To.CKdelta(exp_SBGE_ase))
  
  # Fit fourth-order polynomial logistic regression
  model4 <- glm(Sig ~ poly(delta, degree = 4, raw = TRUE), data = dat, family = "binomial")
  model3 <- glm(Sig ~ poly(delta, degree = 3, raw = TRUE), data = dat, family = "binomial")
  model2 <- glm(Sig ~ poly(delta, degree = 2, raw = TRUE), data = dat, family = "binomial")
  model1 <- glm(Sig ~ poly(delta, degree = 1, raw = TRUE), data = dat, family = "binomial")
  model0 <- glm(Sig ~ 1, data = dat, family = "binomial") # null model
  
  AIC.vec = AIC(model0, model1, model2, model3, model4)
  # plot(x=seq(0,5,by=1), y=AIC.vec)
  
  # requirement "0" for Twin Peaks
  #Does model 4 have the best AIC?
  TP.requirement0 = AIC.vec[5,2] == min(AIC.vec$AIC)
  
  # requirement 1 for Twin Peaks
  # is quartic model significant over null model
  lrtest(model0, model4) # likelihood ratio test
  TP.requirement1 = lrtest(model0, model4)$`Pr(>Chisq)`[2] < 0.05
  
  
  # requirement 2 for Twin Peaks
  # test limits of graph are towards -Inf
  # Is the coefficient negative for the 4th order negative?
  TP.requirement2 = model4$coefficients[5] < 0
  # print(as.numeric(model4$coefficients))
  
  
  ## use this function to check if values are all within the range of interest
  ## the default range is -10, 15
  withinBounds<-function(x, lower = min(dat$delta), 
                         upper = max(dat$delta)) (Re(x )> lower) & ( Re(x) < upper)
  
  # requirement 3 for Twin Peaks
  # 3a) first check if roots are real
  # Find roots of derivative of polynomial to locate maxima and minima
  root.locations = polyroot( (1:4)*as.numeric(model4$coefficients)[2:5] )
  
  # Are the roots real? check that imaginary part is very small
  TP.requirement3a = all(  abs( Im(root.locations) ) < 10^-13  )
  
  # 3b) real roots/maxima & minima fall within bounds
  withinBounds((1:4)*as.numeric(model4$coefficients)[2:5])
  
  TP.requirement3b = all( withinBounds(Re(root.locations)) )
  
  AllRequirementsMet = all (c(TP.requirement0, TP.requirement1, TP.requirement2, TP.requirement3a, TP.requirement3b))
  return(AllRequirementsMet)
}


PermuteTP <- function(dat, perm_N){
  
  print(twinPeaks(dat))
  
  fracTP <- NULL # initialize data frame to store results
    
  n1 <- dim(dat[dat$Sig,])[1] # number of significant/candidate genes
  n2 <- dim(dat[!dat$Sig,])[1] # number of non-significant/non-candidate genes
  
  # permute N times
  for (i in 1:perm_N) {
    
    # Randomly re-assign significance
    perm_sig <- sample(dat$FlyBaseID, n1, replace = F)
    perm_dt <- dat %>% mutate(Sig = ifelse(FlyBaseID %in% perm_sig, TRUE, FALSE))
    
    # test model based on current permutation
    fracTP <- c(fracTP, twinPeaks(perm_dt))
  }
  
  return((sum(fracTP))/(perm_N))
}

perm_A.f <- PermuteTP(A.f.geno_ASE, 10)
perm_A.m <- PermuteTP(A.m.geno_ASE, 10)
perm_All <- PermuteTP(SSAV.geno_ASE, 10)


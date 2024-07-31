###################################
#
#                             Grieshop et al. 2024
#                             Author: Michelle Liu
#             DsRed experimental evolution - transcriptomics analysis
#                  Distribution of SA candidates by sex-bias
#                         Figure 2, Figure S3, Figure S4
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
########

# Get Mishra et al.'s data 
source("Mishra_et.al_SBGE.R")

# Prepare plotting dataset 
# (Run this first before anything else!!, make sure you have the correct files/path to them)
source("DE_Plotting_data.R")


# Figure 2 plotting
# Dot & whisker plot for proportion of candidates vs non-candidates per SBGE bin
########
# function to count the number of DE genes each SBGE category
propSBGE <- function(dat, SBGE_cat){
  total_All <- dim(dat)[1]
  total_Sig <- dim(dat[dat$Sig,])[1]
  total_NS <- total_All - total_Sig
  
  total_SBGE <- dat %>% 
    dplyr::group_by(.[[SBGE_cat]]) %>%
    dplyr::summarise(SBGE_count = n()) %>%
    dplyr::rename({{SBGE_cat}} := 1)
  
  fract <- dat %>% 
    dplyr::group_by(.[[SBGE_cat]], Sig) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::mutate(frac_Sig = ifelse(Sig, count/total_Sig, count/total_NS)) %>%
    dplyr::rename({{SBGE_cat}} := 1) 
  
  fract <- merge(fract, total_SBGE, by = SBGE_cat)
  prop_SIG <- sum(fract[fract$Sig,]$count)/sum(fract[fract$Sig,]$SBGE_count)
  prop_NS <- sum(fract[!fract$Sig,]$count)/sum(fract[!fract$Sig,]$SBGE_count)
  
  fract <- fract %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(frac_by_SBGE_bin = count/SBGE_count,
                  lower = list(binom.test(x=count, n=SBGE_count, p=ifelse(Sig, prop_SIG, prop_NS))), 
                  upper = lower$conf.int[2], 
                  lower = lower$conf.int[1])
  
  return(fract)
}

# plotting function
plotSBGEprop <- function(dat, SBGE_cat, xlab){
  
  # calculate the proportion of genes for each SBGE category
  plot_dat <- propSBGE(dat, SBGE_cat)
  
  # initialize new data frame to store confidence interval
  prop_all <- data.frame(dim(dat[dat$Sig,])[1]/dim(dat)[1])
  # calculate 95% confidence interval
  prop_all$lower <- binom.test(dim(dat[dat$Sig,])[1], dim(dat)[1])$conf.int[1]
  prop_all$upper <- binom.test(dim(dat[dat$Sig,])[1], dim(dat)[1])$conf.int[2]
  
  
  plot_vec <- ggplot(plot_dat[plot_dat$Sig,], aes_string(x = SBGE_cat, y = "frac_by_SBGE_bin")) + 
    geom_hline(yintercept = prop_all[,1],
               linetype = "dashed", show.legend = TRUE, size = 0.75) +
    annotate('ribbon', x = c(-Inf, Inf), ymin = prop_all$lower, ymax = prop_all$upper, 
             alpha = 0.20, fill = 'grey30') +
    geom_errorbar(ymin = plot_dat[plot_dat$Sig,]$lower, ymax = plot_dat[plot_dat$Sig,]$upper,
                  width = 0.5, size = 1, color = "#009E73") +
    geom_point(fill = "#009E73", color = "#009E73", size = 7) + 
    labs(x = "Sex-biased Gene Expression", 
         y = "Proportion of DE Genes") +
    scale_x_discrete(labels = c("H.FB", "FB", "UB", "MB", "H.MB")) +
    scale_y_continuous(limits = c(0, 0.20)) +
    geom_text(data = plot_dat[plot_dat$Sig,], aes_string(x = SBGE_cat, y = "upper" , label = "count" ), 
              size = 7.5, vjust = -0.7) +
    theme_classic() +
    theme(plot.title.position = c("panel"),
          legend.title = element_blank(),
          legend.position = c("bottom"),
          legend.box.background = element_rect(),
          legend.text = element_text(size = 20, color = "black"),
          axis.text.x = element_text(size=20, margin = margin(5,0,5,0), color = "black"),
          axis.text.y = element_text(size=20, margin = margin(0,5,0,0), color = "black"),
          axis.title.x = element_text(size=30, margin = margin(10,0,10,0), color = "black"),
          axis.title.y = element_text(size=30, margin = margin(0,10,0,0), color = "black"),
          plot.title = element_text(size=40, margin = margin(0,0,0,0), color = "black"),
          plot.margin = margin(6,6,6,6))
  return(plot_vec)
}


# Print tables for the proportion DE genes in each SBGE category
propSBGE(SSAV.geno, "SBGE_comp")
propSBGE(A.m.geno, "SBGE_comp")
propSBGE(A.f.geno, "SBGE_comp")

# Fig 2
bin_All <- plotSBGEprop(SSAV.geno, "SBGE_comp", "SBGE (ASE)")

# Fig S3.A
bin_A.m <- plotSBGEprop(A.m.geno, "SBGE_comp", "SBGE (ASE)") + coord_cartesian(ylim = c(0, 0.1))

# Fig S3.B
bin_A.f <- plotSBGEprop(A.f.geno, "SBGE_comp", "SBGE (ASE)")


########


# comment in to save plots
pdf(file = "~/Desktop/UofT/SSAV_RNA/Plots/final_2/Fig2_main_filtered.pdf",   # The directory you want to save the file in
    width = 14, # 14 24 The width of the plot in inches
    height = 10) # 10 20 The height of the plot in inches

# main Figure 2
bin_All + coord_cartesian(ylim = c(0,0.15))


# Figure S3
# ggarrange(bin_A.m, NA, bin_A.f,
#           labels = c("A)", NA, "B)"),
#           widths = c(1, 0.05, 1),
#           ncol = 3,
#           font.label = list(size = 30), hjust = -0.01)

dev.off()




# Repeat analysis, but using Cheng & Kirkpatrick 2016 logistic regression & Twin Peaks model
# convert Log2FC(M/F) to delta
convert.Log2Sex.To.CKdelta <- function(x) (1-2^(-x))/(1+2^(-x))

# fit cubic splines with GAM
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
    labs(y = "Likelihood of Antagonism", x = "") +
    theme_classic() +
    theme(plot.title.position = c("panel"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size=18, margin = margin(0,5,0,0), color = "black"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=24, margin = margin(0,10,0,0), color = "black")) 
  
  return(return.plot)
}

# bootstrap fourth degree polynomial regression
boot_CK_polynomial <- function(dat, boot_N, degree){
  
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


# only with data from Exp. female comparison
######
# bootstrapped 4th degree polynomial version
# A.f.geno.reg_delta <- boot_CK_polynomial(A.f.geno_ASE, 1000, 4)
# SA_dist <- ggplot(A.f.geno.reg_delta) +
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

# cubic spline version (final version displayed in manuscript)
SA_dist <- plotSplinesDelta(A.f.geno) 

SBGE_dist <- ggplot(A.f.geno) + 
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


# only with data from Exp. male comparison
######
# bootstrapped 4th degree polynomial version
# A.m.geno.reg_delta <- boot_reg_delta(A.m.geno_ASE, 1000, 4)
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

# cubic spline version
SA_dist <- plotSplinesDelta(A.m.geno)

SBGE_dist <- ggplot(A.m.geno) + 
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


# combine genes from both assays
######
# bootstrapped 4th degree polynomal version
# SSAV.geno.reg_delta <- boot_reg_delta(SSAV.geno_ASE, 1000, 4)
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

# cubic splines version
SA_dist <- plotSplinesDelta(SSAV.geno)

SBGE_dist <- ggplot(SSAV.geno) + 
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


# Comment in to save plot (Fig. S4)
# pdf(file = "~/Desktop/UofT/SSAV_RNA/Plots/final_2/Fig_S4_CK.pdf",   # The directory you want to save the file in
#     width = 24, # 14 18 24 The width of the plot in inches
#     height = 20) # 10 20 20 The height of the plot in inches
# 
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
# 
# dev.off()



# Function to test for twinPeaks requirements as stated in Cheng and Kirkpatrick (2016)
# (0. Fit with 4th degree polynomial should be better than with a quadratic model)
# 1. The quartic model fits the data better than the null model
# 2. Is the coefficient negative for the 4th order negative?
# 3. Check if maxima and minima are real and within bounds
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

# permute the data and test for twinPeaks pattern n number of times
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


# perm_A.f <- PermuteTP(A.f.geno, 1000)
# perm_A.m <- PermuteTP(A.m.geno, 1000)
# perm_All <- PermuteTP(SSAV.geno, 1000)


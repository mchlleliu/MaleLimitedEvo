###################################
#
#                             Grieshop et al. 2023
#             DsRed experimental evolution - transcriptomics analysis
#             Data were prepared using file DsRed_transcriptomics.txt
# 
# 
###################################

# rm(list=ls()) # Clears the environment 

# set path to data files
Data_path <- "~/Desktop/UofT/PRJ1/ReadCount_KG_renamed/"

# Packages
##########

# if needed
# install.packages("BiocManager")
# BiocManager::install("DESeq2")
# BiocManager::install("apeglm")
# BiocManager::install("ashr")
# BiocManager::install("mltools")

library(DESeq2)
library(apeglm)
library(ashr)
library(vsn)
library(hexbin)
library(ggplot2)
library(ggblend)
library(ggpubr)
library(gridExtra)
library(plyr)
library(dplyr)
library(RColorBrewer)
library(pheatmap)
library(car)
library(lmtest)
library(lmPerm)
library(Rmisc)
library(Hmisc)
library(mltools)
library(boot)
library(clipr) 
library(broom)
library(tidyr)
library(extrafont)
library(scales)
library(readxl)
library(cowplot)
library(ggstatsplot)

##########


# make dataframe for DESeq2
##########

# list the files by category
# the ordering is based on treatment(A/C)_rep(1:6)_sex(F/M)
files <- list.files(Data_path)

# chop off the ".tsv" from those file names
samplename <- gsub('.{4}$', '', files)

# set up factor for condition "sex" (must match order that files are listed; see "samplename" and "files")
sex <- factor(c(rep(c("Female", "Male"), each = 2, times = 6), # Af, Am
                rep("Male", times = 12)))  # Cm

# set up factor for each genotype (Red or NR) "geno" (must match order that files are listed; see "samplename" and "files")
geno <- factor(c(rep(c("NR", "Red"), each = 1, times = 18)))

# set up factor for each treatment (A or C) "trt" (must match order that files are listed; see "samplename" and "files")
trt <- factor(c(rep("A", 24), 
                rep("C", 12)))

# set up factor for each treatment (Af, Am, or Cm) "trt2" 
trt2 <- factor(c(rep(c("Af", "Am"), each = 2, times = 6), 
                 rep("Cm", 12)))

# set up factor for batch "rep" (must match order that files are listed; see "samplename" and "files")
rep <- factor(c(as.character(c(rep(1:6, each = 4, times = 1), # 6 replicates for each of the 4 sex-geno in A
                               rep(1:6, each = 2, times = 1)))))  # 6 replicates for each of the male-geno in C


# put those into a dataframe 
sampleTable <- data.frame(sampleName = samplename, 
                          fileName = files, 
                          sex = sex,
                          geno = geno,
                          trt = trt,
                          trt2 = trt2,
                          rep = rep)

str(sampleTable) # Check that factors are factors, else use e.g.: sampleTable$rep <- factor(sampleTable$rep)

##########
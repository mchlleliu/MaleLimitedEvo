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

# Split data frame up into different subsets 
##########

# trt A data (list of sample names, and files)
A <- sampleTable[(sampleTable$trt == "A"),] 
A$trt <- droplevels(A$trt)

# A.f data
A.f <- A[(A$sex == "Female"),]
A.f$sex <- droplevels(A.f$sex)

# A.m data
A.m <- A[(A$sex == "Male"),]
A.m$sex <- droplevels(A.m$sex)

# C.m
C.m <- sampleTable[(sampleTable$trt == "C"),] # There are only Male samples for C
C.m$trt <- droplevels(C.m$trt)

# A.f.C.m (contrast A females and C males)
A.f.C.m <- rbind(A.f, C.m) # For contrasting A.f to C.m
A.f.C.m_Red <- A.f.C.m[(A.f.C.m$geno == "Red"),] # within Red
A.f.C.m_Red$geno <- droplevels(A.f.C.m_Red$geno)
A.f.C.m_NR <- A.f.C.m[(A.f.C.m$geno == "NR"),] # within NR
A.f.C.m_NR$geno <- droplevels(A.f.C.m_NR$geno)

# A.f.nr_A.m.r (contrast A nonRed females and A red males)
A.f.nr <- A.f[(A.f$geno == "NR"),]
A.m.r <- A.m[(A.m$geno == "Red"),]
A.f.nr_A.m.r <- rbind(A.f.nr, A.m.r)

# A.Red
A.Red <- A[(A$geno == "Red"),]
A.Red$geno <- droplevels(A.Red$geno)

# A.NR
A.NR <- A[(A$geno == "NR"),]
A.NR$geno <- droplevels(A.NR$geno)

# Just Males (A and C)
Males <- sampleTable[(sampleTable$sex == "Male"),] 
Males$sex <- droplevels(Males$sex)

# Just Females (A)
Females <- sampleTable[(sampleTable$sex == "Female"),] 
Females$sex <- droplevels(Females$sex)

# A.Red.m
A.Red.m <- Males[(Males$geno == "Red") & Males$trt == "A",]
A.Red.m$geno <- droplevels(A.Red.m$geno)

# A.NR.m
A.NR.m <- Males[(Males$geno == "NR") & Males$trt == "A",]
A.NR.m$geno <- droplevels(A.NR.m$geno)

# A.Red.f
A.Red.f <- Females[(Females$geno == "Red"),]
A.Red.f$geno <- droplevels(A.Red.f$geno)

# A.NR.f
A.NR.f <- Females[(Females$geno == "NR"),]
A.NR.f$geno <- droplevels(A.NR.f$geno)

# A.Red.m.NR.f
A.Red.m.NR.f <- rbind(A.Red.m, A.NR.f)
A.Red.m.NR.f$geno <- droplevels(A.Red.m.NR.f$geno)

# A.NR.m.Red.f
A.NR.m.Red.f <- rbind(A.NR.m, A.Red.f)
A.NR.m.Red.f$geno <- droplevels(A.NR.m.Red.f$geno)

##########



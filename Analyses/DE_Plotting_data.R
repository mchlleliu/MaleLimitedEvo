###################################
#
#                             Grieshop et al. 2024
#                             Author: Michelle Liu
#             DsRed experimental evolution - transcriptomics analysis
#                         Set up plotting data set
# 
# 
###################################

# packages required
#########
library(tidyr)
library(plyr)
library(dplyr)
library(broom)
########

# load SBGE info from Mishra et
source("Mishra_et.al_SBGE.R")

# Prepare plotting dataset 
# (make sure you have the correct files/path to them)
#########
# load results if not loaded in env.
A.f.geno <- read.delim("Results/A.f.geno_DE.candidates.tsv")
A.m.geno <- read.delim("Results/A.m.geno_DE.candidates.tsv")
C.m.geno <- read.delim("Results/C.m.geno_DE.candidates.tsv")
SSAV.geno <- read.delim("Results/All.geno_DE.candidates.tsv")

# include SBGE categories (using Mishra et al. 2022 dataset)
A.m.geno <- merge(A.m.geno, Mishra, by = "FlyBaseID", all = TRUE)
A.m.geno <- A.m.geno[!is.na(A.m.geno$Sig) & !is.na(A.m.geno$exp_SBGE_ase),]
A.m.geno$SBGE_comp <- as.factor(A.m.geno$SBGE_comp)
A.m.geno$SBGE_simp <- as.factor(A.m.geno$SBGE_simp)
str(A.m.geno)

A.f.geno <- merge(A.f.geno, Mishra, by = "FlyBaseID", all = TRUE)
A.f.geno <- A.f.geno[!is.na(A.f.geno$Sig) & !is.na(A.f.geno$exp_SBGE_ase),]
A.f.geno$SBGE_comp <- as.factor(A.f.geno$SBGE_comp)
A.f.geno$SBGE_simp <- as.factor(A.f.geno$SBGE_simp)
str(A.f.geno)


C.m.geno <- merge(C.m.geno, Mishra, by = "FlyBaseID", all = TRUE)
C.m.geno <- C.m.geno[!is.na(C.m.geno$Sig) & !is.na(C.m.geno$exp_SBGE_ase),]
C.m.geno$SBGE_comp <- as.factor(C.m.geno$SBGE_comp)
C.m.geno$SBGE_simp <- as.factor(C.m.geno$SBGE_simp)
str(C.m.geno)


# Genes present in both SSAV males and SSAV females data
SSAV.geno <- merge(SSAV.geno, Mishra, by = "FlyBaseID", all = TRUE)
SSAV.geno <- SSAV.geno[!is.na(SSAV.geno$Sig) & !is.na(SSAV.geno$exp_SBGE_ase),]
SSAV.geno$SBGE_comp <- as.factor(SSAV.geno$SBGE_comp)
SSAV.geno$SBGE_simp <- as.factor(SSAV.geno$SBGE_simp)
str(SSAV.geno)

#########
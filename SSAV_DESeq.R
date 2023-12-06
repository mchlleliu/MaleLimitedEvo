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
files <- list.files(Data_path)

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


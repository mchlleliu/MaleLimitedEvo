###################################
#
#                             Grieshop et al. 2023
#                             Author: Karl Grieshop
#             DsRed experimental evolution - transcriptomics analysis
#                      Load Mishra et al. 2022 data to assign SBGE
# 
# 
###################################

# packages
##########
library(readxl)
library(plyr)
library(dplyr)
##########

# Get Mishra et al.'s data 
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
# # xmid1 = (x1 + x2)/2

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


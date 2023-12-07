###################################
#
#                             Grieshop et al. 2023
#             DsRed experimental evolution - Analysis into dimorphism
#                     Michelle's code - dds.A.Red.m.NR.f.sex
#                                     - dds.A.NR.m.Red.f.sex    
# 
# 
###################################

# \\||// #
# Separate SBGE categories
##########

SBGE_log2FC <- read.delim(file="~/Desktop/UofT/SSAV_RNA/Results/Results.A.Red.m.NR.fem.sex.txt", sep="\t", header=TRUE)

# Define three levels of SBGE categorization 
x1 = 1 # first cut-off (FBG < -1, MBG > 1 , -1 < UBG < 1)
x2 = 5 # second cut-off (extreme FBG < -5, extreme MBG > 5)
y0 = 1 # tolerance of middle bins ### I DONT GET THIS 
# xmid1 = (x1 + x2)/2

# Simple (3 levels)
# one level of female-biased gene expression
fbg.keep <- SBGE_log2FC$log2FC < -x1 & (SBGE_log2FC$log2FC + SBGE_log2FC$lfcSE) < 0
fbg <- SBGE_log2FC[fbg.keep,]
fbg$SBGE_simp <- rep(c("a.fbg"), dim(fbg)[1])
# one unbiased category
ubg.keep <- SBGE_log2FC$log2FC < x1 & SBGE_log2FC$log2FC > -x1 & 
  (SBGE_log2FC$log2FC - SBGE_log2FC$lfcSE) > -(x1+y0) & (SBGE_log2FC$log2FC + SBGE_log2FC$lfcSE) < (x1+y0)
ubg <- SBGE_log2FC[ubg.keep,]
ubg$SBGE_simp <- rep(c("b.ubg"), dim(ubg)[1])
# two levels of male-biased gene expression
mbg.keep <- SBGE_log2FC$log2FC > x1 & (SBGE_log2FC$log2FC - SBGE_log2FC$lfcSE) > 0
mbg <- SBGE_log2FC[mbg.keep,]
mbg$SBGE_simp <- rep(c("c.mbg"), dim(mbg)[1])
# 1 gene is tossed out b/c the uncertainty in its estimate breaches a cutoff boundary 
SBGE_log2FC <- rbind(fbg, mbg, ubg)
str(SBGE_log2FC)

# Complex (5 levels)
# two levels of female-biased gene expression 
more.fbg.keep <- SBGE_log2FC$log2FC < -x2 & (SBGE_log2FC$log2FC + SBGE_log2FC$lfcSE) < -(x2-y0) # extreme FBG < -5
more.fbg <- SBGE_log2FC[more.fbg.keep,]
more.fbg$SBGE_comp <- rep(c("a.more.fbg"), dim(more.fbg)[1])
#
fbg.keep <- SBGE_log2FC$log2FC < -x1 & SBGE_log2FC$log2FC > -x2 & 
  (SBGE_log2FC$log2FC + SBGE_log2FC$lfcSE) < 0 & (SBGE_log2FC$log2FC - SBGE_log2FC$lfcSE) > -(x2+y0)
fbg <- SBGE_log2FC[fbg.keep,]
fbg$SBGE_comp <- rep(c("b.fbg"), dim(fbg)[1])
# one unbiased category
ubg.keep <- SBGE_log2FC$log2FC < x1 & SBGE_log2FC$log2FC > -x1 & 
  (SBGE_log2FC$log2FC - SBGE_log2FC$lfcSE) > -(x1+y0) & (SBGE_log2FC$log2FC + SBGE_log2FC$lfcSE) < (x1+y0)
ubg <- SBGE_log2FC[ubg.keep,]
ubg$SBGE_comp <- rep(c("c.ubg"), dim(ubg)[1])
# two levels of male-biased gene expression
mbg.keep <- SBGE_log2FC$log2FC > x1 & SBGE_log2FC$log2FC < x2 & 
  (SBGE_log2FC$log2FC - SBGE_log2FC$lfcSE) > 0 & (SBGE_log2FC$log2FC + SBGE_log2FC$lfcSE) < (x2+y0)
mbg <- SBGE_log2FC[mbg.keep,]
mbg$SBGE_comp <- rep(c("d.mbg"), dim(mbg)[1])
#
more.mbg.keep <- SBGE_log2FC$log2FC > x2 & (SBGE_log2FC$log2FC - SBGE_log2FC$lfcSE) > (x2-y0)
more.mbg <- SBGE_log2FC[more.mbg.keep,]
more.mbg$SBGE_comp <- rep(c("e.more.mbg"), dim(more.mbg)[1])
# 3 genes tossed out  b/c the uncertainty in their estimate breaches a cutoff boundary
SBGE_log2FC <- rbind(more.fbg, fbg, ubg, mbg, more.mbg)
str(SBGE_log2FC)

## Save to object name that makes sense.
# A.NR.m.Red.f.SBGE <- SBGE_log2FC
A.Red.m.NR.f.SBGE <- SBGE_log2FC

# write.table(SBGE_log2FC, file = "~/Desktop/UofT/SSAV_RNA/Results/SBGE.A.Red.m.NR.f.tsv", sep = "\t", # Fix file name accordingly
#             row.names = FALSE, col.names = TRUE)
##########


# Get Mishra et al.'s data 
##########

ASE <- read.delim(file="~/Desktop/UofT/PRJ1/data/DifferentialGeneExpression.whole.bodies.tsv", sep="\t", header=TRUE)
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
# y0 = 1 # tolerance of middle bins ### I DONT GET THIS 
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
# two levels of female-biased gene expression (NOTE: I don't get this end term: -(x2-y0))
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


# Define three levels of SBGE categorization based on Singh and Agrawal 2023
x1 = 0.5 # first cut-off 
x2 = 2 # second cut-off 
x3 = 5 # third cut-off
y0 = 1 # tolerance of middle bins ### I DONT GET THIS 
# xmid1 = (x1 + x2)/2

# More Complex (7 levels)
# three levels of female-biased gene expression
ext.fbg.keep <- ASE$exp_SBGE_ase < -x3 & 
  (ASE$exp_SBGE_ase - ASE$se_SBGE_ase) < -(x3-y0) # extreme FBG < -5
ext.fbg <- ASE[ext.fbg.keep,]
ext.fbg$SBGE_singh <- rep(c("a.ext.fbg"), dim(ext.fbg)[1])
#
more.fbg.keep <- ASE$exp_SBGE_ase < -x2 & ASE$exp_SBGE_ase > -x3 & 
  (ASE$exp_SBGE_ase + ASE$se_SBGE_ase) < -(x2+y0) &
  (ASE$exp_SBGE_ase - ASE$se_SBGE_ase) > -(x3+y0) 
more.fbg <- ASE[more.fbg.keep,]
more.fbg$SBGE_singh <- rep(c("b.more.fbg"), dim(more.fbg)[1])
#
fbg.keep <- ASE$exp_SBGE_ase < -x1 & ASE$exp_SBGE_ase > -x2 & 
  (ASE$exp_SBGE_ase + ASE$se_SBGE_ase) < 0 &
  (ASE$exp_SBGE_ase - ASE$se_SBGE_ase) > -(x2+y0)
fbg <- ASE[fbg.keep,]
fbg$SBGE_singh <- rep(c("c.fbg"), dim(fbg)[1])
# one unbiased category
ubg.keep <- ASE$exp_SBGE_ase < x1 & ASE$exp_SBGE_ase > -x1 & 
  (ASE$exp_SBGE_ase - ASE$se_SBGE_ase) > -(x1+y0) & (ASE$exp_SBGE_ase + ASE$se_SBGE_ase) < (x1+y0)
ubg <- ASE[ubg.keep,]
ubg$SBGE_singh <- rep(c("d.ubg"), dim(ubg)[1])
# three levels of male-biased gene expression
mbg.keep <- ASE$exp_SBGE_ase > x1 & ASE$exp_SBGE_ase < x2 & 
  (ASE$exp_SBGE_ase - ASE$se_SBGE_ase) > 0 & 
  (ASE$exp_SBGE_ase + ASE$se_SBGE_ase) < (x2+y0)
mbg <- ASE[mbg.keep,]
mbg$SBGE_singh <- rep(c("e.mbg"), dim(mbg)[1])
#
more.mbg.keep <- ASE$exp_SBGE_ase > x2 & ASE$exp_SBGE_ase < x3 &
  (ASE$exp_SBGE_ase - ASE$se_SBGE_ase) > (x2-y0) & (ASE$exp_SBGE_ase - ASE$se_SBGE_ase) < (x3-y0)
more.mbg <- ASE[more.mbg.keep,]
more.mbg$SBGE_singh <- rep(c("f.more.mbg"), dim(more.mbg)[1])
#
ext.mbg.keep <- ASE$exp_SBGE_ase > x3 & (ASE$exp_SBGE_ase - ASE$se_SBGE_ase) > (x3-y0)
ext.mbg <- ASE[ext.mbg.keep,]
ext.mbg$SBGE_singh <- rep(c("g.ext.mbg"), dim(ext.mbg)[1])
# 3 genes tossed out  b/c the uncertainty in their estimate breaches a cutoff boundary
ASE <- rbind(ext.fbg, more.fbg, fbg, ubg, mbg, more.mbg, ext.mbg)
str(ASE)

##########


# Compare with ASE population
##########
# Load Prashastha's ASE data
ASE <- 

# Merge Natural and Unnatural SBGE estimates.
Sex_Bias_Nat_UnNat <- merge(A.Red.m.NR.f.SBGE, A.NR.m.Red.f.SBGE, by = "FlyBaseID", all = T)
colnames(Sex_Bias_Nat_UnNat) <- c("FlyBaseID", "log2FC_Nat", "lfcSE_Nat", "p_Nat", "SBGE_simp.Nat", "SBGE_comp.Nat",
                                  "log2FC_UnNat", "lfcSE_UnNat", "p_UnNat", "SBGE_simp.UnNat", "SBGE_comp.UnNat")
Sex_Bias_Nat_UnNat$diff <- Sex_Bias_Nat_UnNat$log2FC_Nat - Sex_Bias_Nat_UnNat$log2FC_UnNat
Sex_Bias_Nat_UnNat <- merge(Sex_Bias_Nat_UnNat, ASE, by = "FlyBaseID", all = T)


# compare with ASE population
Sex_Bias_Nat_UnNat$diff_Nat.ASE <- abs(Sex_Bias_Nat_UnNat$log2FC_Nat) - abs(Sex_Bias_Nat_UnNat$exp_SBGE_ase)
Sex_Bias_Nat_UnNat$diff_UnNat.ASE <- abs(Sex_Bias_Nat_UnNat$log2FC_UnNat) - abs(Sex_Bias_Nat_UnNat$exp_SBGE_ase)

# compare with OSADA population
Sex_Bias_Nat_UnNat <- merge(Sex_Bias_Nat_UnNat, SDIU, by = "FlyBaseID", all = T)
Sex_Bias_Nat_UnNat <- subset(Sex_Bias_Nat_UnNat, select=c(1:18,23:25))
Sex_Bias_Nat_UnNat$diff_Nat.OSADA <- abs(Sex_Bias_Nat_UnNat$log2FC_Nat) - abs(Sex_Bias_Nat_UnNat$Whole.SBGE.Osada)
Sex_Bias_Nat_UnNat$diff_UnNat.OSADA <- abs(Sex_Bias_Nat_UnNat$log2FC_UnNat) - abs(Sex_Bias_Nat_UnNat$Whole.SBGE.Osada)


# UnNat_ASE, Nat_ASE, Nat_UnNat, Nat_OSADA, UnNat_OSADA
Nat_UnNat <- ggplot(Sex_Bias_Nat_UnNat[!is.na(Sex_Bias_Nat_UnNat$log2FC_Nat) & !is.na(Sex_Bias_Nat_UnNat$log2FC_UnNat) &
                                         !is.na(Sex_Bias_Nat_UnNat$exp_SBGE_ase),], # & 
                    #Sex_Bias_Nat_UnNat$FlyBaseID %in% All.geno[!is.na(All.geno$Sig) & All.geno$Sig,]$FlyBaseID,],
                    aes(x = exp_SBGE_ase, y = diff)) + 
  geom_point(size = 2, shape = 16, alpha = 0.7, colour = "grey") +
  geom_smooth() +
  scale_x_continuous(breaks = seq(-15, 15, 2.5)) +
  labs(y = "Red.m/NR.f - NR.m/Red.f", x="SBGE(ASE)") +
  geom_hline(yintercept = 0, size = 0.5, linetype= "dashed", color = "black", alpha = 0.7) +
  geom_vline(xintercept = 0, size = 0.5, linetype= "dashed", color = "black", alpha = 0.7)

##########

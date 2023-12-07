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

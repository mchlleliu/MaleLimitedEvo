###################################
#
#                             Grieshop et al. 2023
#                             Altered from: JunctionSeq.Script.R by 
#                       Amardeep Singh -- amardeep.singh[at]utoronto.ca
#                  DsRed experimental evolution - transcriptomics analysis
# 
# 
###################################

require(VennDiagram)
require(grDevices)
require(dplyr)
require(tidyr)
require(ggstatsplot)

# Set global variables
FDRThreshold = 0.1
mappedReadsThreshold = 50


jseq.A.f.geno = read.table("JunctionSeq/A.f.geno/A.f.geno.splicingallGenes.results.txt", 
                           sep = "\t", header = TRUE)
jseq.A.m.geno = read.table("JunctionSeq/A.m.geno/A.m.geno.splicingallGenes.results.txt", 
                           sep = "\t", header = TRUE)
jseq.C.m.geno = read.table("JunctionSeq/C.m.geno/C.m.geno.splicingallGenes.results.txt", 
                           sep = "\t", header = TRUE)
jseq.ASE = read.table("JunctionSeq/SDIU_ase/JSresults/SDIU_ASEallGenes.results.txt",
                      sep = "\t", header = TRUE)

### prep data
########
# Assign SDIU and non-SDIU genes
jseq.C.m.geno <- jseq.C.m.geno %>% 
  dplyr::mutate(sig.hit = ifelse(geneWisePadj <= FDRThreshold, TRUE, FALSE))

jseq.A.f.geno <- jseq.A.f.geno %>% 
  dplyr::mutate(sig.hit = ifelse(geneWisePadj <= FDRThreshold & 
                          !(FlyBaseID %in% jseq.C.m.geno[jseq.C.m.geno$sig.hit,]$FlyBaseID),
                          TRUE, FALSE))
jseq.A.m.geno <- jseq.A.m.geno %>% mutate(sig.hit = ifelse(geneWisePadj <= FDRThreshold &
                                                             !(FlyBaseID %in% jseq.C.m.geno[jseq.C.m.geno$sig.hit,]$FlyBaseID),
                                                           TRUE, FALSE))
jseq.ASE <- jseq.ASE %>% mutate(SSS = ifelse(geneWisePadj <= 0.01, TRUE, FALSE))
# fisher's test for genes identified as SSS in Osada et al. and ASE data
colnames(jseq.ASE)[2]="FlyBaseID"
test<-merge(SDIU, jseq.ASE, by = "FlyBaseID")
fisher.test(test$SDIU.body.sig, test$SSS)


#######


# compare with SDIU data
# cut down the columns not needed
jseq.A.f.geno = jseq.A.f.geno[,c(2,14,25,26)]
jseq.A.m.geno = jseq.A.m.geno[,c(2,14,25,26)]
jseq.C.m.geno = jseq.C.m.geno[,c(2,14,25,26)]
jseq.ASE = jseq.ASE[,c(2,14,25,26)]

jseq.A.f.geno = jseq.A.f.geno[!duplicated(jseq.A.f.geno[1:4]),]
jseq.A.m.geno = jseq.A.m.geno[!duplicated(jseq.A.m.geno[1:4]),]
jseq.C.m.geno = jseq.C.m.geno[!duplicated(jseq.C.m.geno[1:4]),]
jseq.ASE = jseq.ASE[!duplicated(jseq.ASE[1:4]),]

colnames(jseq.A.m.geno)[1] <- "FlyBaseID"
colnames(jseq.A.f.geno)[1] <- "FlyBaseID"
colnames(jseq.C.m.geno)[1] <- "FlyBaseID"
colnames(jseq.ASE)[1]="FlyBaseID"


# SSS in body?
# Osada data
jseq.A.f.geno = jseq.A.f.geno %>% mutate(SSS = ifelse(FlyBaseID %in% SDIU[SDIU$SDIU.body.sig,]$FlyBaseID, TRUE, FALSE))
jseq.A.m.geno = jseq.A.m.geno %>% mutate(SSS = ifelse(FlyBaseID %in% SDIU[SDIU$SDIU.body.sig,]$FlyBaseID, TRUE, FALSE))
jseq.C.m.geno = jseq.C.m.geno %>% mutate(SSS = ifelse(FlyBaseID %in% SDIU[SDIU$SDIU.body.sig,]$FlyBaseID, TRUE, FALSE))

# ASE data
jseq.A.f.geno = jseq.A.f.geno %>% mutate(SSS.ase = ifelse(FlyBaseID %in% jseq.ASE[jseq.ASE$SSS,]$FlyBaseID, TRUE, FALSE))
jseq.A.m.geno = jseq.A.m.geno %>% mutate(SSS.ase = ifelse(FlyBaseID %in% jseq.ASE[jseq.ASE$SSS,]$FlyBaseID, TRUE, FALSE))
jseq.C.m.geno = jseq.C.m.geno %>% mutate(SSS.ase = ifelse(FlyBaseID %in% jseq.ASE[jseq.ASE$SSS,]$FlyBaseID, TRUE, FALSE))

# Significant DE?
jseq.A.f.geno = jseq.A.f.geno %>% mutate(DE.geno = ifelse(FlyBaseID %in% A.f.geno[A.f.geno$Sig,]$FlyBaseID, TRUE, FALSE))
jseq.A.m.geno = jseq.A.m.geno %>% mutate(DE.geno = ifelse(FlyBaseID %in% A.f.geno[A.f.geno$Sig,]$FlyBaseID, TRUE, FALSE))
jseq.C.m.geno = jseq.C.m.geno %>% mutate(DE.geno = ifelse(FlyBaseID %in% A.f.geno[A.f.geno$Sig,]$FlyBaseID, TRUE, FALSE))



dim(na.omit(jseq.A.m.geno[jseq.A.m.geno$sig.hit,]))
dim(na.omit(jseq.A.m.geno[jseq.A.m.geno$sig.hit & jseq.A.m.geno$SSS,]))
dim(na.omit(jseq.A.f.geno[jseq.A.f.geno$sig.hit & jseq.A.f.geno$DE.geno,]))



test <- merge(jseq.All.geno, jseq.ASE, by = "FlyBaseID", all )
test.fisher.results <- fisher.test(test$sig.spliced, test$SSS)

SSS.ase <- ggbarstats(
  test[test$FlyBaseID %in% filter.low.exp.genes.list,], SSS, sig.spliced,
  results.subtitle = FALSE,
  subtitle = paste0(
    "Fisher's exact test", ", p-value = ",
    ifelse(test.fisher.results$p.value < 0.001, "< 0.001", round(test.fisher.results$p.value, 3))
  ), 
  xlab = NULL
) +
  scale_fill_manual(labels = c("SSS", "not SSS"),
                    values = c("darkorchid4", "darkgrey")) + # "Chr-2", "Chr-3", "X-Chr"
  scale_x_discrete(labels = c("ns spliced", "sig spliced")) +
  theme(plot.title.position = c("panel"),
        legend.title = element_blank(),
        legend.position = c("bottom"),
        legend.box.background = element_rect(),
        legend.text = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size=20, margin = margin(5,0,0,0), color = "black"),
        axis.text.y = element_text(size=10, margin = margin(0,5,0,0), color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=40, margin = margin(0,10,0,0), color = "black"),
        plot.title = element_text(size=40, margin = margin(0,0,0,0), color = "black"),
        plot.margin = margin(6,6,6,6)
  )



pdf(file = "~/Desktop/UofT/SSAV_RNA/Plots/splicing/ASE.SSS.v.RedNR.spliced.pdf",  # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 7) # The height of the plot in inches
SSS.ase # mos_plot_InnoMorr_Chr2, mos_plot_Ruz, mos_plot_Ruz_Chr2, mos_plot_Chloe_Chr2, mos_plot_Wong, mos_plot_Wong_Chr2
dev.off()


# Venn Diagrams 
#######
# separate significant and non-significant genes
sig.genes.A.f.geno = unique(jseq.A.f.geno$FlyBaseID[jseq.A.f.geno$geneWisePadj < FDRThreshold &
                                                      !is.na(jseq.A.f.geno$geneWisePadj)])
nonsig.genes.A.f.geno = unique(jseq.A.f.geno$FlyBaseID[jseq.A.f.geno$geneWisePadj > FDRThreshold &
                                                         !is.na(jseq.A.f.geno$geneWisePadj)])

sig.genes.A.m.geno = unique(jseq.A.m.geno$FlyBaseID[jseq.A.m.geno$geneWisePadj < FDRThreshold &
                                                      !is.na(jseq.A.m.geno$geneWisePadj)])
nonsig.genes.A.m.geno = unique(jseq.A.m.geno$FlyBaseID[jseq.A.m.geno$geneWisePadj > FDRThreshold &
                                                         !is.na(jseq.A.m.geno$geneWisePadj)])

sig.genes.C.m.geno = unique(jseq.C.m.geno$FlyBaseID[jseq.C.m.geno$geneWisePadj < FDRThreshold &
                                                      !is.na(jseq.C.m.geno$geneWisePadj)])
nonsig.genes.C.m.geno = unique(jseq.C.m.geno$FlyBaseID[jseq.C.m.geno$geneWisePadj > FDRThreshold &
                                                         !is.na(jseq.C.m.geno$geneWisePadj)])

venndiagram = venn.diagram(x = list(sig.genes.A.f.geno, sig.genes.A.m.geno, nonsig.genes.A.f.geno, nonsig.genes.A.m.geno),
                           category.names = c("A.f Sig", "A.m Sig", "A.f NonSig", "A.m NonSig"),
                           filename = NULL,
                           # Circles
                           lwd = 2, lty = 'blank', fill = c("#7294d4", "#72D481", "#D4B272", "#D472C5"),
                           # Numbers
                           cex = 1.5, fontface = "bold", fontfamily = "Helvetica",
)
pdf(file="Plots/vennJS_FDR01.pdf", height = 10, width = 10)
grid.draw(venndiagram)
dev.off()
#######


# SSAV FPKM check

dds.A.f.geno
load("JunctionSeq/SDIU_ase/Drosophila_melanogaster.BDGP6.28.102.exonLengths.RData", verbose = T)

exonic <- GRangesList(exonic)
rowRanges(dds.A.f.geno) <- exonic # add gene length data to DESeq2 object
rowRanges(dds.A.m.geno) <- exonic 
rowRanges(dds.C.m.geno) <- exonic

FPKM.A.f <- fpkm(dds.A.f.geno) # get FPKM measures
FPKM.A.f <- data.frame(FPKM.A.f) %>%
  dplyr::mutate(totalCounts = rowSums(select_if(., is.numeric)))
test <- as.matrix(FPKM.A.f[,1:12])
testVarAf <- matrixStats::rowVars(test)
testVarAf <- data.frame(testVarAf)
testVarAf$trt2 <- "Af"
colnames(testVarAf)[1] <- "Var"
FPKM.A.f$trt2 <- "Af"


FPKM.A.m <- fpkm(dds.A.m.geno) # get FPKM measures
FPKM.A.m <- data.frame(FPKM.A.m) %>%
  dplyr::mutate(totalCounts = rowSums(select_if(., is.numeric)))
test <- as.matrix(FPKM.A.m[,1:12])
testVarAm <- matrixStats::rowVars(test)
testVarAm <- data.frame(testVarAm)
testVarAm$trt2 <- "Am"
colnames(testVarAm)[1] <- "Var"
FPKM.A.m$trt2 <- "Am"


FPKM.C.m <- fpkm(dds.C.m.geno) # get FPKM measures
FPKM.C.m <- data.frame(FPKM.C.m) %>%
  dplyr::mutate(totalCounts = rowSums(select_if(., is.numeric)))
test <- as.matrix(FPKM.C.m[,1:12])
testVar <- matrixStats::rowVars(test)
testVar <- data.frame(testVar)
testVar$trt2 <- "Cm" 
colnames(testVar)[1] <- "Var"
FPKM.C.m$trt2 <- "Cm"

quantile(FPKM.C.m$totalCounts)
quantile(FPKM.A.m$totalCounts)
quantile(FPKM.A.f$totalCounts)

# test <- rbind(FPKM.A.f[,c("totalCounts", "trt2")], 
#               FPKM.A.m[,c("totalCounts", "trt2")],
#               FPKM.C.m[,c("totalCounts", "trt2")])
test <- rbind(testVarAf, testVarAm, testVar)
ggplot(test,  aes(trt2, Var)) +  geom_boxplot() + 
  coord_cartesian(ylim = c(0, 10000))



library(GEOquery)
library(affy)
library(limma)
library(tidyverse)


# load series and platform data from GEO

gset <- getGEO("GSE17013", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GSE17013", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# filePaths = getGEOSuppFiles("GSE17013")
# varLabels(gset)
# untar("GSE17013/GSE17013_RAW.tar", exdir = "GSE17013/")
cels = list.files("GSE17013/", pattern = "CEL")
# sapply(paste("GSE17013/", cels, sep = "/"), gunzip)

cels = paste0("GSE17013/",cels)

raw.data = ReadAffy(verbose = FALSE, filenames = cels)
data.rma.norm = rma(raw.data)


# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
sex <- rep(c(rep("female", 4), rep("male", 4)), 15)
line <- c("H07", "H09", "H10", "H12", "H13", "H14", 
          "P07", "P13", "P15", "P18", "P22", "P26", "P48", "P50", "P56")
rep <- rep(c(1:4), 30)
for(i in line) lines = c(lines, rep(i, 8))
lines = lines[-1]
design_exp <- data.frame(sex = sex, line = unlist(lines), rep = rep)


design <- model.matrix(~0+sex+line, design_exp)

dupcor <- duplicateCorrelation(gset, design, block=rep)
fit <- lmFit(gset, design, block = rep, correlation = dupcor$consensus)  # fit linear model


# set up contrasts of interest and recalculate model coefficients
cts <- paste("sexmale", "sexfemale", sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=18952)
summary(decideTests(fit2,lfc=1))


InnoMorrow_SBGE <- subset(tT, select=c("Platform_CLONEID", "Chromosome.annotation", "logFC", "AveExpr", "adj.P.Val"))
str(InnoMorrow_SBGE)
colnames(InnoMorrow_SBGE)[1] = "FlyBaseID"

InnoMorrow_SBGE <- InnoMorrow_SBGE[!(InnoMorrow_SBGE$FlyBaseID %in% Chrs$FlyBaseID[Chrs$Chr == "X"]) &
                                     !(InnoMorrow_SBGE$FlyBaseID %in% Chrs$FlyBaseID[Chrs$Chr == "Y"]),] %>%
  dplyr::mutate(Sig = ifelse(FlyBaseID %in% Inno_Morrow$FlyBaseID, TRUE, FALSE)) 

write.table(InnoMorrow_SBGE, "Data/InnoMorrow.SBGE.tsv", sep = "\t", quote = F)


### compare with SSAV DE analysis dataset
SSAV.geno <- read.delim("Results/All.geno_candidates.tsv")
A.f.geno <- read.delim("Results/A.f.geno_candidates.tsv")
A.m.geno <- read.delim("Results/A.m.geno_candidates.tsv")

SSAV.geno_IM <- merge(Results.df, InnoMorrow_SBGE, by = "FlyBaseID", all = TRUE)
A.f.geno_IM <- merge(A.f.geno, InnoMorrow_SBGE, by = "FlyBaseID", all = TRUE)
A.m.geno_IM <- merge(A.m.geno, InnoMorrow_SBGE, by = "FlyBaseID", all = TRUE)


SSAV.geno_IM <- InnoMorrow_SBGE %>%
  mutate(Sig = ifelse(FlyBaseID %in% SSAV.geno$FlyBaseID[SSAV.geno$Sig], TRUE, FALSE))

SSAV.geno_IM <- SSAV.geno_IM[!is.na(SSAV.geno_IM$Sig) & !is.na(SSAV.geno_IM$logFC),]
str(SSAV.geno_IM)
A.f.geno_IM <- A.f.geno_IM[!is.na(A.f.geno_IM$Sig) & !is.na(A.f.geno_IM$logFC),]
A.m.geno_IM <- A.m.geno_IM[!is.na(A.m.geno_IM$Sig) & !is.na(A.m.geno_IM$logFC),]

# Cheng & Kirkpatrick 2016 logistic regression & Twin Peaks model
# convert Log2FC(M/F) to delta
convert.Log2Sex.To.CKdelta<-function(x) (1-2^(-x))/(1+2^(-x))


SSAV.geno.IM_delta <- test_bootReg(SSAV.geno_IM, 50, 4)


A.f.geno.IM_delta <- test_bootReg(A.f.geno_IM, 1000, 4)
A.m.geno.IM_delta <- test_bootReg(A.m.geno_IM, 1000, 4)
test_delta <- test_bootReg(test_IM, 10, 4)
dim(test_IM)
# plotting all CK plot delta
######
SA_dist <- ggplot(test_delta) + 
  geom_ribbon(aes(x= delta, ymin=q05, ymax=q95),
              size = 0.25, alpha = 0.7, fill = "grey")  +
  geom_line(aes(x=delta, y=q50), size = 1) +
  labs(y="P(Sexually Antagonistic)") +
  theme_classic() +
  theme(plot.title.position = c("panel"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=18, margin = margin(0,5,0,0), color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=24, margin = margin(0,10,0,0), color = "black")) 


SBGE_dist <- ggplot(InnoMorrow_SBGE[InnoMorrow_SBGE$AveExpr > 3,]) + 
  geom_density(aes(x=sapply(logFC, convert.Log2Sex.To.CKdelta), y= ..count..), fill = "grey", color = "grey") +
  labs(x=expression(Delta),y="Number of loci") +
  theme_classic() +
  theme(axis.text.x = element_text(size=18, margin = margin(5,0,5,0), color = "black"),
        axis.text.y = element_text(size=12, margin = margin(0,5,0,0), color = "black"),
        axis.title.x = element_text(size=30, margin = margin(10,0,10,0), color = "black"),
        axis.title.y = element_text(size=15, margin = margin(0,20,0,12.5), color = "black")) 

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

Fig2_CKsuppl_IM
Af_CKsuppl_IM
Am_CKsuppl_IM <- ggarrange(SA_dist + scale_x_continuous(expand = c(0.0005, 0.0005)), NA,
                              NA, NA,
                              SBGE_dist + scale_x_continuous(expand = c(0.0005, 0.0005), 
                                                             labels = c("-1.0", "-0.5", "0", "0.5", "1.0")), NA,
                              nrow = 3, heights = c(1, 0.005, 0.30), ncol = 2, widths = c(1, 0.05)) 


compareCK <- ggarrange(Fig2_CKsuppl_all + ggtitle("ASE SBGE") + theme(plot.title = element_text(size = 18, margin=margin(5,5,5,5))), 
          Fig2_CKsuppl_IM + ggtitle("IM SBGE") + theme(plot.title = element_text(size = 18, margin=margin(5,5,5,5))))
######


pdf(file = "~/Desktop/UofT/SSAV_RNA/Plots/AfAm_CK.pdf",   # The directory you want to save the file in
    width = 17, # 14 24 The width of the plot in inches
    height = 10) # 10 20 The height of the plot in inches
ggarrange(Am_CKsuppl_IM, NA, Af_CKsuppl_IM, 
            labels = c("B)", NA, "C)"),
           ncol = 3,
           widths = c(1, 0.05, 1), font.label = list(size = 30), hjust = -0.01)
dev.off()

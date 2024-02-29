###################################
#
#                                   Grieshop et al. 2023
#                  DsRed experimental evolution - transcriptomics analysis
#                 Differential splicing analysis for Red vs Non-Red samples
# 
# 
###################################


# required packages
#######
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(DESeq2)
library(GenomicFeatures)
require(VennDiagram)
require(grDevices)
require(ggstatsplot)
library(broom)
library(ggblend)
#######

# Set global variables
FDRThreshold = 0.01
mappedReadsThreshold = 50


# load data from JunctionSeq analysis
# to run JunctionSeq, see: Amardeep Singh's JunctionSeq.Script.R,
# or Michelle's modified version for these (SSAV) populations
#######
# ASE reference
jseq.ASE = read.table("JunctionSeq/SDIU_ase/JSresults/SDIU_ASEallGenes.results.txt",
                      sep = "\t", header = TRUE)
colnames(jseq.ASE)[2]="FlyBaseID" # change column name to reflect the rest of the dataset

# SSAV samples
jseq.A.f.geno = read.table("JunctionSeq/A.f.geno/A.f.geno.splicingallGenes.results.txt", 
                           sep = "\t", header = TRUE)
jseq.A.m.geno = read.table("JunctionSeq/A.m.geno/A.m.geno.splicingallGenes.results.txt", 
                           sep = "\t", header = TRUE)
jseq.C.m.geno = read.table("JunctionSeq/C.m.geno/C.m.geno.splicingallGenes.results.txt", 
                           sep = "\t", header = TRUE)
colnames(jseq.A.f.geno)[2]="FlyBaseID" 
colnames(jseq.A.m.geno)[2]="FlyBaseID" 
colnames(jseq.C.m.geno)[2]="FlyBaseID" 
# make a list of samples. This is just so we can automate 
# it easier in subsequent steps bcs i don't want to write everything thrice lol T-T.
SSAV.sample.types <- data.frame(sampleType = c("A.m", "A.f", "C.m")) %>% 
  mutate(raw.JS.table = paste0("jseq.",sampleType,".geno"))
######


# SSAV sample types data table:
SSAV.sample.types


# prepare data: filtering un-assayed genes, assign significance
#######
# ------ First for ASE reference 
# remove genes that were not assayed in males and females
jseq.ASE = jseq.ASE[!is.na(jseq.ASE$expr_F) & !is.na(jseq.ASE$expr_M),]

# remove novel counting bins
#   for the purpose of comparing ASE with the SSAV populations, we want to make sure that the counting bins
#   used are the same.
jseq.ASE = jseq.ASE %>% filter(!str_detect(countbinID, "N"))
dim(jseq.ASE) # check how many novel splice sites were removed

# asign genes with significant sex-specific splicing
jseq.ASE <- jseq.ASE %>% mutate(SSS = ifelse(geneWisePadj < 0.01, TRUE, FALSE))

ASE.sig.SSS <- unique(jseq.ASE$FlyBaseID[jseq.ASE$geneWisePadj < 0.01])
ASE.nonsig.SSS <-  unique(jseq.ASE$FlyBaseID[jseq.ASE$geneWisePadj >= 0.01])
length(ASE.sig.SSS)
length(ASE.nonsig.SSS)

# fisher's test for genes identified as SSS in Osada et al. and P. Mishra's ASE population
SinghAgrawalSDIU <- read.csv(file="~/Desktop/UofT/SSAV_RNA/Data/SBGEandSSSdataForMBE.csv", sep=",", header=TRUE)
colnames(SinghAgrawalSDIU)[2] <- "FlyBaseID"
test<-SinghAgrawalSDIU %>% mutate(ASE.SSS = ifelse(FlyBaseID %in% ASE.sig.SSS, TRUE, FALSE))
fisher.test(test$SDIU.body.sig, test$ASE.SSS)
rm(test, SinghAgrawalSDIU)


# ------ For SSAV samples
# set FDRThreshold to a less stringent one for the SSAV data
FDRThreshold = 0.1
for(i in SSAV.sample.types$raw.JS.table){
  tmp.JS.table <- get(paste0(i))
  # remove untested genes
  tmp.JS.table = tmp.JS.table[!is.na(tmp.JS.table$expr_Red) & !is.na(tmp.JS.table$expr_NR),]
  tmp.JS.table = tmp.JS.table[!is.na(tmp.JS.table$geneWisePadj),]
  
  print(dim(tmp.JS.table)) # check number of counting bins
  # remove novel splice sites
  tmp.JS.table = tmp.JS.table %>% filter(!str_detect(countbinID, "N"))
  print(dim(tmp.JS.table)) # check how many novel splice sites were removed
  
  tmp.JS.table <- tmp.JS.table %>% mutate(sig.hit = ifelse(geneWisePadj <= FDRThreshold, TRUE, FALSE))
  
  assign(i, tmp.JS.table)
  rm(tmp.JS.table)
}

# Assign significantly spliced genes 
jseq.C.m.geno <- jseq.C.m.geno %>% mutate(sig.hit = ifelse(geneWisePadj <= FDRThreshold, TRUE, FALSE))
# do not assign as significantly spliced if the controls also differ
jseq.A.f.geno <- jseq.A.f.geno %>% mutate(sig.hit = ifelse(geneWisePadj <= FDRThreshold & 
                                                             !(FlyBaseID %in% jseq.C.m.geno[jseq.C.m.geno$sig.hit,]$FlyBaseID),
                                                           TRUE, FALSE))
jseq.A.m.geno <- jseq.A.m.geno %>% mutate(sig.hit = ifelse(geneWisePadj <= FDRThreshold &
                                                             !(FlyBaseID %in% jseq.C.m.geno[jseq.C.m.geno$sig.hit,]$FlyBaseID),
                                                           TRUE, FALSE))
jseq.All.geno <- merge(jseq.A.m.geno, jseq.A.f.geno, by = c("FlyBaseID", "countbinID"), all = T) %>%
  mutate(sig.hit = ifelse((!(is.na(sig.hit.x) & is.na(sig.hit.y)) & (sig.hit.x | sig.hit.y)), 
                          TRUE, ifelse(is.na(sig.hit.x) & is.na(sig.hit.y), NA, FALSE)))
SSAV.sample.types <- rbind(SSAV.sample.types, c("A.all", "jseq.All.geno"))

#######


# filtering based on ASE reference Male and Female lowest FPKM
######
# load metadata file
decoder.ASE <- read.table("JunctionSeq/SDIU_ase/QoRTs.decoder.file.for.JunctionSeq.txt", header=T, stringsAsFactors = F)
colnames(decoder.ASE)[1]="unique.ID"
decoder.ASE$rep = rep(seq(1:3),4)

# make list containing the paths to sample count files
countFiles.ASE <- paste0("JunctionSeq/SDIU_ase/count.files/",decoder.ASE$unique.ID,"/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")

# use this function to read and combine all sample counts and create a counts matrix
# data.files is a list containing the paths to sample count files
makeCountMatrix <- function(data.files, decoder.file){
  count.file <- NULL
  # load and combine sample counts 
  for(i in data.files){
    tmp <- read.delim(i, header=F, sep = "\t")
    tmp <- tmp %>% separate(V1, into = c("FlyBaseID", "countBin"), sep = ":") %>%
      filter(str_detect(countBin, "A000"))
    if(!is.null(count.file)){
      count.file <- merge(count.file, tmp, by = "FlyBaseID", all = T)}
    else
      count.file <- tmp
  }
  
  colnames(count.file) <- c("FlyBaseID", decoder.file$unique.ID)
  
  # remove bins that overlap multiple genes
  count.file <- count.file %>% filter(!str_detect(FlyBaseID, "\\+")) 
  
  # change geneIDs to rownames
  count.file <- count.file %>% remove_rownames %>% column_to_rownames(var="FlyBaseID")
  
  return(count.file)
}

ASE.count.matrix <- makeCountMatrix(countFiles.ASE, decoder.ASE)

# get column metadata from decoder file
ASE.colData <- decoder.ASE %>% remove_rownames %>% column_to_rownames(var="unique.ID")
# make DESeq2 object to calculate fpkm
dds.ASE <- DESeqDataSetFromMatrix(countData = ASE.count.matrix, 
                                  colData = ASE.colData, 
                                  design = ~ rep + sex)


# get gene lengths from GTF file 
# (done on server, scp to local so load this to env if you also did it that way)
# if not, just use the "exonic" GRanges object
# txdb <- makeTxDbFromGFF(gtfFile, format="gtf")
# exonic <- exonsBy(txdb, by="gene")
# save(exonic, file="/plas1/michelle.liu/Dmel_BDGP6.28/Drosophila_melanogaster.BDGP6.28.102.exonLengths.RData")
load("JunctionSeq/SDIU_ase/Drosophila_melanogaster.BDGP6.28.102.exonLengths.RData", verbose = T)

exonic <- GRangesList(exonic)
exonic <- exonic[names(exonic) %in% rownames(ASE.count.matrix)] # keep only gene lengths that are in the count matrix
rowRanges(dds.ASE) <- exonic # add gene length data to DESeq2 object
FPKM.ASE <- fpkm(dds.ASE) # get FPKM measures

# separate the counts to males and females
FPKM.ASE.fem <- data.frame(FPKM.ASE) %>% select(., contains("F_")) %>%
  mutate(totalCounts = rowSums(select_if(., is.numeric))) # get total counts for females
FPKM.ASE.fem[FPKM.ASE.fem==0] <- NA # set genes not present as NAs

FPKM.ASE.male <- data.frame(FPKM.ASE) %>% select(., contains("M_")) %>%
  mutate(totalCounts = rowSums(select_if(., is.numeric))) # get total counts for males
FPKM.ASE.male[FPKM.ASE.male==0] <- NA # set genes not present as NAs

# list of genes with FPKM > the 25% cut-off in females
filter.low.exp.genes.fem.q25 <- rownames(FPKM.ASE.fem[!is.na(FPKM.ASE.fem$totalCounts) &
                                                        FPKM.ASE.fem$totalCounts  > quantile(FPKM.ASE.fem$totalCounts, 0.25, na.rm=T),])
# list of genes with FPKM > the 25% cut-off in males
filter.low.exp.genes.male.q25 <- rownames(FPKM.ASE.male[!is.na(FPKM.ASE.male$totalCounts) &
                                                          FPKM.ASE.male$totalCounts  > quantile(FPKM.ASE.male$totalCounts, 0.25, na.rm=T),])
# combine list of genes that passed filtering
filter.low.exp.genes.q25 <- unique(c(filter.low.exp.genes.fem.q25, filter.low.exp.genes.male.q25))
# remove Y chr genes that somehow gets there(?) could be from spermatheca in females?
filter.low.exp.genes.q25 <- filter.low.exp.genes.q25[!(filter.low.exp.genes.q25 %in% Chrs[Chrs$Chr == "Y",]$FlyBaseID)]
length(filter.low.exp.genes.q25) # check how many genes are left


# list of genes with FPKM > the 10% cut-off in females
filter.low.exp.genes.fem.q10 <- rownames(FPKM.ASE.fem[!is.na(FPKM.ASE.fem$totalCounts) &
                                                        FPKM.ASE.fem$totalCounts  > quantile(FPKM.ASE.fem$totalCounts, 0.10, na.rm=T),])
# list of genes with FPKM > the 10% cut-off in males
filter.low.exp.genes.male.q10 <- rownames(FPKM.ASE.male[!is.na(FPKM.ASE.male$totalCounts) &
                                                          FPKM.ASE.male$totalCounts  > quantile(FPKM.ASE.male$totalCounts, 0.10, na.rm=T),])
# combine list of genes that passed filtering
filter.low.exp.genes.q10 <- unique(c(filter.low.exp.genes.fem.q10, filter.low.exp.genes.male.q10))
# remove Y chr genes that somehow gets there(?) could be from spermatheca in females?
filter.low.exp.genes.q10 <- filter.low.exp.genes.q10[!(filter.low.exp.genes.q10 %in% Chrs[Chrs$Chr == "Y",]$FlyBaseID)]
length(filter.low.exp.genes.q10) # check how many genes are left
######


# results:
length(filter.low.exp.genes.q10) # list of genes that passed the 10% filter
length(filter.low.exp.genes.q25) # list of genes that passed the 25% filter




# differential splicing analysis for SSAV samples
#######
# downstream checks for significantly spliced genes

# check how many genes were spliced differently between Red and NR, 
# do some Fisher's exact tests for association with sex-specific splicing in SSS data or in ASE data:
# load external Singh & Agrawal 2023 list of SSS genes (Osada, 2017 populations)
SinghAgrawalSDIU <- read.csv(file="~/Desktop/UofT/SSAV_RNA/Data/SBGEandSSSdataForMBE.csv", sep=",", header=TRUE)
colnames(SinghAgrawalSDIU)[2] <- "FlyBaseID"
# results from differential expression analysis
All.geno.tmp <- read.delim("Results/All.geno_candidates.tsv", header = TRUE, sep = "\t")

# initialize dataframe object to store results (wow that variable name do be long)
fishers.test.RedNR.splice.results <- data.frame(sampleType = SSAV.sample.types$sampleType)

for(i in 1:dim(SSAV.sample.types)[1]){
  tmp.JS.table <- get(paste0(SSAV.sample.types[i, 2]))
  tmp.JS.table <- tmp.JS.table[!is.na(tmp.JS.table$sig.hit),]
  if(i < 4) 
    tmp.JS.table <- tmp.JS.table[,c(2,26)]
  else 
    tmp.JS.table <- tmp.JS.table[,c(1, 51)]
  
  tmp.JS.table = tmp.JS.table[!duplicated(tmp.JS.table$FlyBaseID),]
  # get number of significant genes
  fishers.test.RedNR.splice.results$N.all[i] = length(tmp.JS.table$FlyBaseID)
  fishers.test.RedNR.splice.results$N.sig[i] = length(tmp.JS.table$FlyBaseID[!is.na(tmp.JS.table$sig.hit) & tmp.JS.table$sig.hit])
  
  # assign overlap with SSS
  tmp.JS.table = tmp.JS.table %>% mutate(Osada.SSS = ifelse(FlyBaseID %in% SinghAgrawalSDIU[!is.na(SinghAgrawalSDIU$SDIU.body.sig) & SinghAgrawalSDIU$SDIU.body.sig,]$FlyBaseID, 
                                                            TRUE, ifelse(is.na(SinghAgrawalSDIU$SDIU.body.sig), NA, FALSE)))
  tmp.JS.table = tmp.JS.table %>% mutate(ASE.SSS = ifelse(FlyBaseID %in% jseq.ASE[!is.na(jseq.ASE$SSS) & jseq.ASE$SSS,]$FlyBaseID, 
                                                          TRUE, ifelse(is.na(jseq.ASE$SSS), NA, FALSE)))
  
  # any overlap with differentially expressed genes?
  tmp.JS.table = tmp.JS.table %>% mutate(DE.Sig = ifelse(FlyBaseID %in% All.geno.tmp[!is.na(All.geno.tmp$Sig) & All.geno.tmp$Sig,]$FlyBaseID, 
                                                         TRUE, ifelse(is.na(All.geno.tmp$Sig), NA, FALSE)))
  
  
  fishers.test.RedNR.splice.results$N.filter10[i] = length(tmp.JS.table[tmp.JS.table$FlyBaseID %in% filter.low.exp.genes.q10,]$FlyBaseID)
  fishers.test.RedNR.splice.results$N.sig.filter10[i] = length(tmp.JS.table$FlyBaseID[tmp.JS.table$FlyBaseID %in% filter.low.exp.genes.q10 &
                                                                                        !is.na(tmp.JS.table$sig.hit) & tmp.JS.table$sig.hit])
  # test for SSS overlap, 10% FPKM filter
  fishers.test.RedNR.splice.results$Osada.SSS.filter10[i] <- fisher.test(tmp.JS.table[tmp.JS.table$FlyBaseID %in% filter.low.exp.genes.q10,]$sig.hit, 
                                                                         tmp.JS.table[tmp.JS.table$FlyBaseID %in% filter.low.exp.genes.q10,]$Osada.SSS)$p.value
  fishers.test.RedNR.splice.results$ASE.SSS.filter10[i] <- fisher.test(tmp.JS.table[tmp.JS.table$FlyBaseID %in% filter.low.exp.genes.q10,]$sig.hit, 
                                                                       tmp.JS.table[tmp.JS.table$FlyBaseID %in% filter.low.exp.genes.q10,]$ASE.SSS)$p.value
  
  # test for DE overlap, 10% FPKM filter
  fishers.test.RedNR.splice.results$SSAV.DE.Sig.filter10[i] <- fisher.test(tmp.JS.table[tmp.JS.table$FlyBaseID %in% filter.low.exp.genes.q10,]$sig.hit, 
                                                                           tmp.JS.table[tmp.JS.table$FlyBaseID %in% filter.low.exp.genes.q10,]$DE.Sig)$p.value
  
  
  fishers.test.RedNR.splice.results$N.filter25[i] = length(tmp.JS.table[tmp.JS.table$FlyBaseID %in% filter.low.exp.genes.q25,]$FlyBaseID)
  fishers.test.RedNR.splice.results$N.sig.filter25[i] = length(tmp.JS.table$FlyBaseID[tmp.JS.table$FlyBaseID %in% filter.low.exp.genes.q25 &
                                                                                        !is.na(tmp.JS.table$sig.hit) & tmp.JS.table$sig.hit])
  # test for SSS overlap, 25% FPKM filter
  fishers.test.RedNR.splice.results$Osada.SSS.filter25[i] <- fisher.test(tmp.JS.table[tmp.JS.table$FlyBaseID %in% filter.low.exp.genes.q25,]$sig.hit, 
                                                                         tmp.JS.table[tmp.JS.table$FlyBaseID %in% filter.low.exp.genes.q25,]$Osada.SSS)$p.value
  fishers.test.RedNR.splice.results$ASE.SSS.filter25[i] <- fisher.test(tmp.JS.table[tmp.JS.table$FlyBaseID %in% filter.low.exp.genes.q25,]$sig.hit, 
                                                                       tmp.JS.table[tmp.JS.table$FlyBaseID %in% filter.low.exp.genes.q25,]$ASE.SSS)$p.value
  
  # test for DE overlap, 25% FPKM filter
  fishers.test.RedNR.splice.results$SSAV.DE.Sig.filter25[i] <- fisher.test(tmp.JS.table[tmp.JS.table$FlyBaseID %in% filter.low.exp.genes.q25,]$sig.hit, 
                                                                           tmp.JS.table[tmp.JS.table$FlyBaseID %in% filter.low.exp.genes.q25,]$DE.Sig)$p.value
  
  print(length(tmp.JS.table$FlyBaseID[tmp.JS.table$FlyBaseID %in% filter.low.exp.genes.q10 &
                                        !is.na(tmp.JS.table$Osada.SSS) & tmp.JS.table$Osada.SSS]))
  print(length(tmp.JS.table$FlyBaseID[tmp.JS.table$FlyBaseID %in% filter.low.exp.genes.q10 &
                                        !is.na(tmp.JS.table$ASE.SSS) & tmp.JS.table$ASE.SSS]))
  print(length(tmp.JS.table$FlyBaseID[tmp.JS.table$FlyBaseID %in% filter.low.exp.genes.q10 &
                                        !is.na(tmp.JS.table$DE.Sig) & tmp.JS.table$DE.Sig]))
  
  print(length(tmp.JS.table$FlyBaseID[tmp.JS.table$FlyBaseID %in% filter.low.exp.genes.q25 &
                                        !is.na(tmp.JS.table$Osada.SSS) & tmp.JS.table$Osada.SSS]))
  print(length(tmp.JS.table$FlyBaseID[tmp.JS.table$FlyBaseID %in% filter.low.exp.genes.q25 &
                                        !is.na(tmp.JS.table$ASE.SSS) & tmp.JS.table$ASE.SSS]))
  print(length(tmp.JS.table$FlyBaseID[tmp.JS.table$FlyBaseID %in% filter.low.exp.genes.q25 &
                                        !is.na(tmp.JS.table$DE.Sig) & tmp.JS.table$DE.Sig]))
  
  
  assign(paste0(SSAV.sample.types[i, 1],".tmp.fisher"), tmp.JS.table)
}

rm(tmp.JS.table, SinghAgrawalSDIU, All.geno.tmp) # clean environment

# save table
write.table(fishers.test.RedNR.splice.results, file = "Results/splicing.Red.NR.fisher.tests.csv", sep = ",", quote = FALSE, row.names = F)

test.fisher.results <- fisher.test(A.m.tmp.fisher$Osada.SSS, A.m.tmp.fisher$ASE.SSS)

SSS.sdiu <- ggbarstats(
  A.m.tmp.fisher, ASE.SSS, Osada.SSS,
  results.subtitle = FALSE,
  subtitle = paste0(
    "Fisher's exact test", ", p-value = ",
    ifelse(test.fisher.results$p.value < 0.001, "< 0.001", round(test.fisher.results$p.value, 3))
  ), 
  xlab = NULL
) +
  scale_fill_manual(labels = c("SSS ASE", "not SSS ASE"),
                    values = c("darkorchid4", "darkgrey")) + # "Chr-2", "Chr-3", "X-Chr"
  scale_x_discrete(labels = c("not SSS OSADA","SSS OSADA")) +
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

pdf(file = "~/Desktop/UofT/SSAV_RNA/Plots/splicing/ASE.SSS.v.RedNR.spliced.cut.pdf",  # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 7) # The height of the plot in inches
SSS.ase # mos_plot_InnoMorr_Chr2, mos_plot_Ruz, mos_plot_Ruz_Chr2, mos_plot_Chloe_Chr2, mos_plot_Wong, mos_plot_Wong_Chr2
dev.off()

#######



# for the genes with significant alt. splicing in Red and Non-Red,
# are they more masculinized or feminized?
A.f.sig.AS <- unique(jseq.A.f.geno$FlyBaseID[!is.na(jseq.A.f.geno$sig.hit) & jseq.A.f.geno$sig.hit])
A.m.sig.AS <- unique(jseq.A.m.geno$FlyBaseID[!is.na(jseq.A.m.geno$sig.hit) & jseq.A.m.geno$sig.hit])


# --------------- analyse differences in splicing profiles
# using count.files (read counts per exon, generated by QoRTs. See JunctionSeq_Linux.sh)
# this is a decoder file with $1:unique.sample.ID, $2:cdt1, $3:cdt2, $n:cdtn (See JunctionSeq_Linux.sh on how to make this)
decoder <- read.table("JunctionSeq/QoRTs.decoder.file.for.JunctionSeq.txt", header = TRUE, stringsAsFactors = FALSE)

decoder.ASE <- read.table("JunctionSeq/SDIU_ase/QoRTs.decoder.file.for.JunctionSeq.txt", header=T, stringsAsFactors = F)
colnames(decoder.ASE)[1]="unique.ID"

decoder.SinghAgrawal <- read.table("JunctionSeq/SDIU_singh/QoRTs.decoder.SinghAgrawal.body.txt", header = FALSE, sep = "\t")
colnames(decoder.SinghAgrawal) <- c("unique.ID", "sex")

decoder.MCsim <- read.table("JunctionSeq/SDIU_MC/MCsim/decoder.MCsim.txt", header = F, sep="\t")
colnames(decoder.MCsim) <- c("unique.ID", "sex")

decoder.MCcom <- read.table("JunctionSeq/SDIU_MC/MCcom/decoder.MCcom.txt", header=F, sep="\t")
colnames(decoder.MCcom) <- c("unique.ID", "sex")

decoder.MCabs <- read.table("JunctionSeq/SDIU_MC/MCabs/decoder.MCabs.txt", header=F, sep="\t")
colnames(decoder.MCabs) <- c("unique.ID", "sex")


# load JunctionSeq sizeFactors (for per sample normalization factors)
# this has to be created through JunctionSeq. 
# You can technically calculate the sizeFactors yourself, but JunctionSeq can give you the by countingBin and by gene factors so it is easier.
# look at DESeq2's geometric normalization calculation if you want to do it yourself
#######
SSAV.factors <- read.delim("JunctionSeq/RAL.size.Factors.GEO.txt", header=T, sep="\t")
A.f.factors <- SSAV.factors  %>% filter(str_detect(sample.ID, "\\_F_"))
A.m.factors <- SSAV.factors  %>% filter(str_detect(sample.ID, "A[[:digit:]]\\_M_"))
C.m.factors <- SSAV.factors  %>% filter(str_detect(sample.ID, "C[[:digit:]]\\_M_"))

ASE.factors <- read.delim("JunctionSeq/SDIU_ase/RAL.size.Factors.GEO.txt", header = T, sep="\t")
F.ASE.factors <- ASE.factors  %>% filter(str_detect(sample.ID, "\\F_"))
M.ASE.factors <- ASE.factors  %>% filter(str_detect(sample.ID, "\\M_"))

SinghAgrawal.factors <- read.delim("JunctionSeq/SDIU_singh/RAL.size.Factors.GEO.txt", header=T, sep="\t")
F.SinghAgrawal.factors <- SinghAgrawal.factors %>% filter(str_detect(sample.ID, "\\.female"))
M.SinghAgrawal.factors <- SinghAgrawal.factors %>% filter(str_detect(sample.ID, "\\.male"))


MC.factors <- read.delim("JunctionSeq/SDIU_MC/JS.GEO.size.factors.txt", header=T, sep="\t")
F.MCabs.factors <- MC.factors %>% filter(str_detect(sample.ID, "M[[:digit:]]\\_Female_body"))
M.MCabs.factors <- MC.factors %>% filter(str_detect(sample.ID, "M[[:digit:]]\\_Male_body"))
F.MCsim.factors <- MC.factors %>% filter(str_detect(sample.ID, "P[[:digit:]]\\_Female_body"))
M.MCsim.factors <- MC.factors %>% filter(str_detect(sample.ID, "P[[:digit:]]\\_Male_body"))
F.MCcom.factors <- MC.factors %>% filter(str_detect(sample.ID, "C[[:digit:]]\\_Female_body"))
M.MCcom.factors <- MC.factors %>% filter(str_detect(sample.ID, "C[[:digit:]]\\_Male_body"))

#######


# separate sample groups & merge with sizeFactors
#######
# females
A.f.decoder <- decoder[decoder$sex == "F",]
A.f.decoder <- merge(A.f.decoder, A.f.factors, by = 1)
A.f.Red.decoder <- A.f.decoder[A.f.decoder$geno == "Red",]
A.f.NR.decoder <- A.f.decoder[A.f.decoder$geno == "NR",]

# males
A.m.decoder <- decoder[decoder$sex == "M" & decoder$pop == "A",]
A.m.decoder <- merge(A.m.decoder, A.m.factors, by = 1)
A.m.Red.decoder <- A.m.decoder[A.m.decoder$geno == "Red",]
A.m.NR.decoder <- A.m.decoder[A.m.decoder$geno == "NR",]

# control males
C.m.decoder <- decoder[decoder$sex == "M" & decoder$pop == "C",]
C.m.decoder <- merge(C.m.decoder, C.m.factors, by = 1)
C.m.Red.decoder <- C.m.decoder[C.m.decoder$geno == "Red",]
C.m.NR.decoder <- C.m.decoder[C.m.decoder$geno == "NR",]


# location of count files
countFiles.A.f <- paste0("JunctionSeq/count.files/", A.f.decoder$unique.ID, "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")
countFiles.A.f.Red <- paste0("JunctionSeq/count.files/", A.f.Red.decoder$unique.ID, "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")
countFiles.A.f.NR <- paste0("JunctionSeq/count.files/", A.f.NR.decoder$unique.ID, "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")

countFiles.A.m <- paste0("JunctionSeq/count.files/", A.m.decoder$unique.ID, "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")
countFiles.A.m.Red <- paste0("JunctionSeq/count.files/", A.m.Red.decoder$unique.ID, "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")
countFiles.A.m.NR <- paste0("JunctionSeq/count.files/", A.m.NR.decoder$unique.ID, "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")

countFiles.C.m <- paste0("JunctionSeq/count.files/", C.m.decoder$unique.ID, "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")
countFiles.C.m.Red <- paste0("JunctionSeq/count.files/", C.m.Red.decoder$unique.ID, "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")
countFiles.C.m.NR <- paste0("JunctionSeq/count.files/", C.m.NR.decoder$unique.ID, "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")


# ASE data for reference
M.decoder.ASE <- decoder.ASE[decoder.ASE$sex=="M",]
M.decoder.ASE <- merge(M.decoder.ASE, M.ASE.factors, by = 1)

F.decoder.ASE <- decoder.ASE[decoder.ASE$sex=="F",]
F.decoder.ASE <- merge(F.decoder.ASE, F.ASE.factors, by =1)

countFiles.m.ASE <- paste0("JunctionSeq/SDIU_ase/count.files/",M.decoder.ASE$unique.ID,"/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")
countFiles.f.ASE <- paste0("JunctionSeq/SDIU_ase/count.files/",F.decoder.ASE$unique.ID,"/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")




# Singh & Agrawal data for another reference
M.decoder.SinghAgrawal <- decoder.SinghAgrawal[decoder.SinghAgrawal$sex=="male",]
M.decoder.SinghAgrawal <- merge(M.decoder.SinghAgrawal, M.SinghAgrawal.factors, by = 1)

F.decoder.SinghAgrawal <- decoder.SinghAgrawal[decoder.SinghAgrawal$sex=="female",]
F.decoder.SinghAgrawal <- merge(F.decoder.SinghAgrawal, F.SinghAgrawal.factors, by = 1)

countFiles.m.SinghAgrw <- paste0("JunctionSeq/SDIU_singh/body.only.replicate.1/",M.decoder.SinghAgrawal$unique.ID,"/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")
countFiles.f.SinghAgrw <- paste0("JunctionSeq/SDIU_singh/body.only.replicate.1/",F.decoder.SinghAgrawal$unique.ID,"/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")



# MC data
# monogamy
M.decoder.MCabs <- decoder.MCabs[decoder.MCabs$sex == "Male",]
M.decoder.MCabs <- merge(M.decoder.MCabs, M.MCabs.factors, by = 1)
F.decoder.MCabs <- decoder.MCabs[decoder.MCabs$sex == "Female",]
F.decoder.MCabs <- merge(F.decoder.MCabs, F.MCabs.factors, by = 1)

countFiles.m.MCabs <- paste0("JunctionSeq/SDIU_MC/MCabs/",M.decoder.MCabs$unique.ID,"/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt")
countFiles.f.MCabs <- paste0("JunctionSeq/SDIU_MC/MCabs/",F.decoder.MCabs$unique.ID,"/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt")

# simple polygamy
M.decoder.MCsim <- decoder.MCsim[decoder.MCsim$sex == "Male",]
M.decoder.MCsim <- merge(M.decoder.MCsim, M.MCsim.factors, by = 1)
F.decoder.MCsim <- decoder.MCsim[decoder.MCsim$sex == "Female",]
F.decoder.MCsim <- merge(F.decoder.MCsim, F.MCsim.factors, by = 1)[-1,]

countFiles.m.MCsim <- paste0("JunctionSeq/SDIU_MC/MCsim/",M.decoder.MCsim$unique.ID,"/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt")
countFiles.f.MCsim <- paste0("JunctionSeq/SDIU_MC/MCsim/",F.decoder.MCsim$unique.ID,"/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt")

# complex polygamy
M.decoder.MCcom <- decoder.MCcom[decoder.MCcom$sex == "Male",]
M.decoder.MCcom <- merge(M.decoder.MCcom, M.MCcom.factors, by = 1)
F.decoder.MCcom <- decoder.MCcom[decoder.MCcom$sex == "Female",]
F.decoder.MCcom <- merge(F.decoder.MCcom, F.MCcom.factors, by = 1)

countFiles.m.MCcom <- paste0("JunctionSeq/SDIU_MC/MCcom/",M.decoder.MCcom$unique.ID,"/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt")
countFiles.f.MCcom <- paste0("JunctionSeq/SDIU_MC/MCcom/",F.decoder.MCcom$unique.ID,"/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt")


#######


# ------ analysis Functions:
# data.files is a list containing the paths to sample count files
# this function is used inside GeomNormCounts
ReadCountFiles <- function(data.files, decoder.file){
  count.file <- NULL
  for(i in data.files){
    tmp <- read.delim(i, header=F, sep = "\t")
    count.file <- cbind(count.file, tmp[,2])
  }
  
  countBins <- tmp %>% separate(V1, into = c("FlyBaseID", "countBin"), sep = ":")
  count.file <- cbind(countBins[,1:2], count.file)
  colnames(count.file) <- c("FlyBaseID", "countBin", decoder.file$unique.ID)
  print(dim(count.file)) # check number of splice sites
  
  # exclude novel splice sites
  count.file <- count.file %>% filter(!str_detect(countBin, "N") & !str_detect(countBin, "A"))
  print(dim(count.file)) # check number of remaining splice sites
  
  return(count.file)
}

# normalize expression by sample, and by gene
# change the index of the column name for the right sizeFactors!!!
GeomNormCounts <- function(data.files, decoder.file){
  
  count.file <- ReadCountFiles(data.files, decoder.file)
  
  # first normalize by sample to account for technical biases due to the sequencing
  # size factor for each sample is the median of the count/geometric mean of each bin. 
  # calculated using JunctionSeq. Look at DESeq2's "geometric" normalization method
  # change the index of the column name for the right sizeFactors!!!
  for(i in 1:dim(decoder.file)[1]){
    count.file[,i+2] <- count.file[,i+2] / decoder.file$size.factor[i] 
  }
  
  # sum the total 
  count.file <- count.file %>% mutate(total=rowSums(select_if(., is.numeric), na.rm = T)) # sum per row
  print(dim(count.file)) # check number of splice sites
  
  # calculate the global gene expression across all samples
  glob.exp <- count.file %>% group_by(FlyBaseID) %>% # sum per gene
    summarise(glob.exp = sum(total, na.rm = T))
  
  # calculate proportion of reads by gene for each counting bin
  count.file <- merge(count.file, glob.exp, by = "FlyBaseID", all = T) %>%
    mutate(frac.exp.per.gene = ifelse(!is.na(glob.exp), total/glob.exp, NA))
  print(dim(count.file)) # make sure nothing gets removed
  
  # set genes or counting bins with no expression to NA instead of 0s
  
  
  return(count.file)
}

# compute euclidean distance between expression of two sample types
# this function uses EuclideanDistance and PercentSimilarity to compute 3 metrics:
#   Euclidean distance between the profiles of the two samples being compared
#   percent similarity of the profiles
#   percent dissimilarity (reciprocal of % sim)
compareSplicingProfiles <- function(d1, d2){
  tmp <- merge(d1, d2, by = c("FlyBaseID", "countBin")) # only get counting bins present in both
  tmp <- tmp %>%
    group_by(FlyBaseID) %>%
    summarise(percent.sim = PercentSimilarity(frac.exp.per.gene.x, frac.exp.per.gene.y),
              percent.dissim = 1 - percent.sim,
              dist = EuclideanDistance(frac.exp.per.gene.x, frac.exp.per.gene.y))
  return(tmp)
}

# Euclidean distance estimate per gene
# v1 = list of fraction in expression of each counting bin within a gene for sample 1
# v2 = list of fraction in expression of each counting bin within a gene for sample 2
EuclideanDistance <- function(v1, v2){
  euc.dist = NA
  if(length(v1) == length(v2) & abs(1-sum(v1, na.rm = T)) < 0.0001 & abs(1-sum(v2, na.rm = T)) < 0.0001){
    n = length(v1)
    sumV = 0
    for(i in 1:n) {
      if(is.na(v1[i])) v1[i] = 0 # if the counting bin is NA, assign 0
      if(is.na(v2[i])) v2[i] = 0 # if the counting bin is NA, assign 0
      # if(!is.na(v1[i]) & !(is.na(v2[i]))){ 
        sumV = sumV + ((v1[i] - v2[i])^2)
    }
    euc.dist = sqrt(sumV)
  }
  return(euc.dist)
}

# Percent Similarity index per gene
# v1 = list of fraction in expression of each counting bin within a gene for sample 1
# v2 = list of fraction in expression of each counting bin within a gene for sample 2
PercentSimilarity<-function(v1, v2){
  x = NA  ## will return NA if v1 and v2 are unequal lengths or either does not sum to 1
  if(length(v1) == length(v2) & abs(1-sum(v1, na.rm = T)) < 0.0001 & abs(1-sum(v2, na.rm = T)) < 0.0001) {
    n = length(v1)
    x = 0
    for(i in 1:n) {
      # counts need to exist in both samples, otherwise do not count towards the similarity index
      if(!is.na(v1[i]) & !(is.na(v2[i]))){ 
        x = x + min(v1[i],v2[i])}
    }
  }
  return(x)
}


CompareDistances <- function(dst1, dst2){ 
  tmp <- merge(dst1, dst2, by = "FlyBaseID") %>%
    mutate(M.sub.F = (dist.y - dist.x)/(dist.x+dist.y),
           M.sub.F = ifelse(is.nan(M.sub.F), 0, M.sub.F),
           M.sim.F = (percent.sim.x - percent.sim.y)/(percent.sim.x+percent.sim.y),
           M.dis.F = (percent.dissim.y - percent.dissim.x)/(percent.dissim.y + percent.dissim.x))
  return(tmp)
}



# average females expression
########
A.f.Red.norm.exp <- GeomNormCounts(countFiles.A.f.Red, A.f.Red.decoder)
A.f.NR.norm.exp <- GeomNormCounts(countFiles.A.f.NR, A.f.NR.decoder)

# average female in SSAV pop
fem.exp <- merge(A.f.Red.norm.exp, A.f.NR.norm.exp, by = c("FlyBaseID", "countBin")) %>%
  group_by(FlyBaseID, countBin) %>%
  summarise(total = (total.x + total.y)/2,
            glob.exp = (glob.exp.x + glob.exp.y)/2,
            frac.exp.per.gene = (frac.exp.per.gene.x + frac.exp.per.gene.y)/2)

########

# average (A) male expression
########
# normalize exp. per gene for each counting bin 
# in Red males
A.m.Red.norm.exp <- GeomNormCounts(countFiles.A.m.Red, A.m.Red.decoder)
# in Non-Red males
A.m.NR.norm.exp <- GeomNormCounts(countFiles.A.m.NR, A.m.NR.decoder)

# average male in SSAV pop
male.exp <- merge(A.m.Red.norm.exp, A.m.NR.norm.exp, by = c("FlyBaseID", "countBin")) %>%
  group_by(FlyBaseID, countBin) %>%
  summarise(total = (total.x + total.y)/2,
         glob.exp = (glob.exp.x + glob.exp.y)/2,
         frac.exp.per.gene = (frac.exp.per.gene.x + frac.exp.per.gene.y)/2)

#######

# average (C) male expression
########
# expression of NR C males
C.m.Red.norm.exp <- GeomNormCounts(countFiles.C.m.Red, C.m.Red.decoder)
C.m.NR.norm.exp <- GeomNormCounts(countFiles.C.m.NR, C.m.NR.decoder)

male.ctr.exp <- merge(C.m.Red.norm.exp, C.m.NR.norm.exp, by = c("FlyBaseID", "countBin")) %>%
  group_by(FlyBaseID, countBin) %>%
  summarise(total = (total.x + total.y)/2,
            glob.exp = (glob.exp.x + glob.exp.y)/2,
            frac.exp.per.gene = (frac.exp.per.gene.x + frac.exp.per.gene.y)/2)
########


# external references
######
# using ASE data
male.exp_ASE <- GeomNormCounts(countFiles.m.ASE, M.decoder.ASE)
fem.exp_ASE <- GeomNormCounts(countFiles.f.ASE, F.decoder.ASE)

# Singh & Agrawal data
male.exp_SinghAgrw <- GeomNormCounts(countFiles.m.SinghAgrw, M.decoder.SinghAgrawal)
fem.exp_SinghAgrw <- GeomNormCounts(countFiles.f.SinghAgrw, F.decoder.SinghAgrawal)

# Monogamy data
male.exp_MCabs <- GeomNormCounts(countFiles.m.MCabs, M.decoder.MCabs)
fem.exp_MCabs <- GeomNormCounts(countFiles.f.MCabs, F.decoder.MCabs)

# Simple Polygamy data
male.exp_MCsim <- GeomNormCounts(countFiles.m.MCsim, M.decoder.MCsim)
fem.exp_MCsim <- GeomNormCounts(countFiles.f.MCsim, F.decoder.MCsim)

# Complex Polygamy data
male.exp_MCcom <- GeomNormCounts(countFiles.m.MCcom, M.decoder.MCcom)
fem.exp_MCcom <- GeomNormCounts(countFiles.f.MCcom, F.decoder.MCcom)
######



# distance between Red and Non-Red estimates
######
A.m.Red.v.NR <- compareSplicingProfiles(A.m.Red.norm.exp, A.m.NR.norm.exp)
C.m.Red.v.NR <- compareSplicingProfiles(C.m.Red.norm.exp, C.m.NR.norm.exp)
A.f.Red.v.NR <- compareSplicingProfiles(A.f.Red.norm.exp, A.f.NR.norm.exp)
mean(C.m.Red.v.NR$percent.dissim, na.rm = T)
mean(A.m.Red.v.NR$percent.dissim, na.rm = T)
mean(A.f.Red.v.NR$percent.dissim, na.rm = T)
######



# comparing dimorphism between populations
# correlation sanity checks comparing references and SSAV males (females)
#######
# test distance from males of ASE population to males from SSAV population
test <- merge(male.exp_ASE, male.exp, by = c("FlyBaseID", "countBin"), all = T)
cor.test(test$frac.exp.per.gene.x, test$frac.exp.per.gene.y)
plot(test$frac.exp.per.gene.x, test$frac.exp.per.gene.y)

# test distance from females of ASE population to males from SSAV population
test <- merge(fem.exp_ASE, fem.exp, by = c("FlyBaseID", "countBin"), all = T)
cor.test(test$frac.exp.per.gene.x, test$frac.exp.per.gene.y)
plot(test$frac.exp.per.gene.x, test$frac.exp.per.gene.y)

# the same thing above, but comparing Osada population and SSAV
test <- merge(male.exp_SinghAgrw, male.exp, by = c("FlyBaseID", "countBin"), all = T)
cor.test(test$frac.exp.per.gene.x, test$frac.exp.per.gene.y)
plot(test$frac.exp.per.gene.x, test$frac.exp.per.gene.y)

test <- merge(fem.exp_SinghAgrw, fem.exp, by = c("FlyBaseID", "countBin"), all = T)
cor.test(test$frac.exp.per.gene.x, test$frac.exp.per.gene.y)
plot(test$frac.exp.per.gene.x, test$frac.exp.per.gene.y)





# sanity test for correlation between SSAV and ASE data M to F distance
# average distance SSAV males vs females
A.MvF.dist <- compareSplicingProfiles(fem.exp, male.exp)
mean(A.MvF.dist$percent.dissim, na.rm = T)
mean(A.MvF.dist[A.MvF.dist$FlyBaseID %in% subset.sss,]$percent.dissim, na.rm = T)
# average ASE males vs females
ASE.MvF.dist <- compareSplicingProfiles(fem.exp_ASE, male.exp_ASE)
mean(ASE.MvF.dist$percent.dissim, na.rm = T)
mean(ASE.MvF.dist[ASE.MvF.dist$FlyBaseID %in% subset.sss,]$percent.dissim, na.rm = T)
# average Osada males vs females
SinghAgrawal.MvF.dist <- compareSplicingProfiles(fem.exp_SinghAgrw, male.exp_SinghAgrw)
mean(SinghAgrawal.MvF.dist$percent.dissim, na.rm = T)
mean(SinghAgrawal.MvF.dist[SinghAgrawal.MvF.dist$FlyBaseID %in% subset.sss,]$percent.dissim, na.rm = T)

MCsim.MvF.dist <- compareSplicingProfiles(fem.exp_MCsim, male.exp_MCsim)
mean(MCsim.MvF.dist$percent.dissim, na.rm = T)
mean(MCsim.MvF.dist[MCsim.MvF.dist$FlyBaseID %in% subset.sss,]$percent.dissim, na.rm = T)





# test M/F correlation between populations...
test <- merge(A.MvF.dist, ASE.MvF.dist, by = "FlyBaseID")
cor.test(test$percent.sim.x, test$percent.sim.y)
plot(test$percent.sim.x, test$percent.sim.y, xlab = "SSAV", ylab ="ASE")
plot(test[test$FlyBaseID %in% test.filter.25,]$percent.sim.x, 
     test[test$FlyBaseID %in% test.filter.25,]$percent.sim.y, xlab = "SSAV", ylab ="ASE")


test <- merge(ASE.MvF.dist, SinghAgrawal.MvF.dist, by = "FlyBaseID")
cor.test(test$percent.sim.x, test$percent.sim.y)
plot(test$percent.sim.x, test$percent.sim.y, xlab = "ASE", ylab ="Osada")
plot(test[test$FlyBaseID %in% test.filter.25,]$percent.sim.x, 
     test[test$FlyBaseID %in% test.filter.25,]$percent.sim.y, xlab = "ASE", ylab ="Osada")


test <- merge(A.MvF.dist, SinghAgrawal.MvF.dist, by = "FlyBaseID")
cor.test(test$percent.sim.x, test$percent.sim.y)
plot(test$percent.sim.x, test$percent.sim.y, xlab = "SSAV", ylab ="Osada")
plot(test[test$FlyBaseID %in% test.filter.25,]$percent.dissim.x, 
     test[test$FlyBaseID %in% test.filter.25,]$percent.dissim.y, xlab = "SSAV", ylab ="Osada")






# ASE vs Osada populations
# ASE males to Osada females
Male.ASE.SinghAgrw.M.dist.f <- compareSplicingProfiles(male.exp_ASE, fem.exp_SinghAgrw)
mean(Male.ASE.SinghAgrw.M.dist.f$dist, na.rm = T)
# ASE males to Osada males
Male.ASE.SinghAgrw.M.dist.m <- compareSplicingProfiles(male.exp_ASE, male.exp_SinghAgrw)
mean(Male.ASE.SinghAgrw.M.dist.m$dist, na.rm = T)
# ASE females to Osada females
Fem.ASE.SinghAgrw.M.dist.f <- compareSplicingProfiles(fem.exp_ASE, fem.exp_SinghAgrw)
mean(Fem.ASE.SinghAgrw.M.dist.f$dist, na.rm = T)
# ASE females to Osada males
Fem.ASE.SinghAgrw.M.dist.m <- compareSplicingProfiles(fem.exp_ASE, male.exp_SinghAgrw) 
mean(Fem.ASE.SinghAgrw.M.dist.m$dist, na.rm = T)

ASE.SDIU.compare.male <- CompareDistances(Male.ASE.SinghAgrw.M.dist.m, Male.ASE.SinghAgrw.M.dist.f )
ASE.SDIU.compare.fem <- CompareDistances(Fem.ASE.SinghAgrw.M.dist.m, Fem.ASE.SinghAgrw.M.dist.f)
splicing.MF.metric.plot(RedData = ASE.SDIU.compare.male[ASE.SDIU.compare.male$FlyBaseID %in% subset.sss,], 
                        NRData = ASE.SDIU.compare.fem[ASE.SDIU.compare.fem$FlyBaseID %in% subset.sss,], 
                        plotCol = "M.dis.F",
                        colour_red = "steelblue3", colour_NR = "magenta4")


# Osada vs ASE populations
# ASE males to Osada females
Male.SinghAgrw.ASE.M.dist.f <- compareSplicingProfiles(male.exp_SinghAgrw, fem.exp_ASE)
mean(Male.SinghAgrw.ASE.M.dist.f$percent.dissim, na.rm = T)
mean(Male.SinghAgrw.ASE.M.dist.f[Male.SinghAgrw.ASE.M.dist.f$FlyBaseID %in% subset.sss,]$percent.dissim, na.rm = T)
# ASE males to Osada males
Male.SinghAgrw.ASE.M.dist.m <- compareSplicingProfiles(male.exp_SinghAgrw, male.exp_ASE)
mean(Male.SinghAgrw.ASE.M.dist.m$percent.dissim, na.rm = T)
mean(Male.SinghAgrw.ASE.M.dist.m[Male.SinghAgrw.ASE.M.dist.m$FlyBaseID %in% subset.sss,]$percent.dissim, na.rm = T)
# ASE females to Osada females
Fem.SinghAgrw.ASE.M.dist.f <- compareSplicingProfiles(fem.exp_SinghAgrw, fem.exp_ASE)
mean(Fem.SinghAgrw.ASE.M.dist.f$percent.dissim, na.rm = T)
mean(Fem.SinghAgrw.ASE.M.dist.f[Fem.SinghAgrw.ASE.M.dist.f$FlyBaseID %in% subset.sss,]$percent.dissim, na.rm = T)
# ASE females to Osada males
Fem.SinghAgrw.ASE.M.dist.m <- compareSplicingProfiles(fem.exp_SinghAgrw, male.exp_ASE)
mean(Fem.SinghAgrw.ASE.M.dist.m$percent.dissim, na.rm = T)
mean(Fem.SinghAgrw.ASE.M.dist.m[Fem.SinghAgrw.ASE.M.dist.m$FlyBaseID %in% subset.sss,]$percent.dissim, na.rm = T)

SDIU.ASE.compare.male <- CompareDistances(Male.SinghAgrw.ASE.M.dist.m, Male.SinghAgrw.ASE.M.dist.f )
SDIU.ASE.compare.fem <- CompareDistances(Fem.SinghAgrw.ASE.M.dist.m, Fem.SinghAgrw.ASE.M.dist.f)
splicing.MF.metric.plot(RedData = SDIU.ASE.compare.male[SDIU.ASE.compare.male$FlyBaseID %in% subset.sss,], 
                        NRData = SDIU.ASE.compare.fem[SDIU.ASE.compare.fem$FlyBaseID %in% subset.sss,], 
                        plotCol = "M.dis.F",
                        colour_red = "steelblue3", colour_NR = "magenta4")


Male.MCsim.ASE.M.dist.f <- compareSplicingProfiles(male.exp_MCsim, fem.exp_ASE)
mean(Male.MCsim.ASE.M.dist.f$percent.dissim, na.rm = T)
mean(Male.MCsim.ASE.M.dist.f[Male.MCsim.ASE.M.dist.f$FlyBaseID %in% subset.sss,]$percent.dissim, na.rm = T)
# ASE males to Osada males
Male.MCsim.ASE.M.dist.m <- compareSplicingProfiles(male.exp_MCsim, male.exp_ASE)
mean(Male.MCsim.ASE.M.dist.m$percent.dissim, na.rm = T)
mean(Male.MCsim.ASE.M.dist.m[Male.MCsim.ASE.M.dist.m$FlyBaseID %in% subset.sss,]$percent.dissim, na.rm = T)
# ASE females to Osada females
Fem.MCsim.ASE.M.dist.f <- compareSplicingProfiles(fem.exp_MCsim, fem.exp_ASE)
mean(Fem.MCsim.ASE.M.dist.f$percent.dissim, na.rm = T)
mean(Fem.MCsim.ASE.M.dist.f[Fem.MCsim.ASE.M.dist.f$FlyBaseID %in% subset.sss,]$percent.dissim, na.rm = T)
# ASE females to Osada males
Fem.MCsim.ASE.M.dist.m <- compareSplicingProfiles(fem.exp_MCsim, male.exp_ASE)
mean(Fem.MCsim.ASE.M.dist.m$percent.dissim, na.rm = T)
mean(Fem.MCsim.ASE.M.dist.m[Fem.MCsim.ASE.M.dist.m$FlyBaseID %in% subset.sss,]$percent.dissim, na.rm = T)

# Osada males to MCsim females
Male.SinghAgrw.MCsim.M.dist.f <- compareSplicingProfiles(male.exp_SinghAgrw, fem.exp_MCsim)
mean(Male.SinghAgrw.MCsim.M.dist.f$dist, na.rm = T)
# Osada males to MCsim males
Male.SinghAgrw.MCsim.M.dist.m <- compareSplicingProfiles(male.exp_SinghAgrw, male.exp_MCsim)
mean(Male.SinghAgrw.MCsim.M.dist.m$dist, na.rm = T)
# OSada females to MCsim females
Fem.SinghAgrw.MCsim.M.dist.f <- compareSplicingProfiles(fem.exp_SinghAgrw, fem.exp_MCsim)
mean(Fem.SinghAgrw.MCsim.M.dist.f$dist, na.rm = T)
# ASE females to Osada males
Fem.SinghAgrw.MCsim.M.dist.m <- compareSplicingProfiles(fem.exp_SinghAgrw, male.exp_MCsim)
mean(Fem.SinghAgrw.MCsim.M.dist.m$dist, na.rm = T)

SDIU.MCsim.compare.male <- CompareDistances(Male.SinghAgrw.MCsim.M.dist.m, Male.SinghAgrw.MCsim.M.dist.f )
SDIU.MCsim.compare.fem <- CompareDistances(Fem.SinghAgrw.MCsim.M.dist.m, Fem.SinghAgrw.MCsim.M.dist.f)
splicing.MF.metric.plot(RedData = SDIU.MCsim.compare.male[SDIU.MCsim.compare.male$FlyBaseID %in% subset.sss,], 
                        NRData = SDIU.MCsim.compare.fem[SDIU.MCsim.compare.fem$FlyBaseID %in% subset.sss,], 
                        plotCol = "M.dis.F",
                        colour_red = "steelblue3", colour_NR = "magenta4")


Male.MCsim.SinghAgrw.M.dist.m <- compareSplicingProfiles(male.exp_MCsim, male.exp_SinghAgrw)
mean(Male.MCsim.SinghAgrw.M.dist.m$percent.dissim, na.rm = T)
mean(Male.MCsim.SinghAgrw.M.dist.m[Male.MCsim.SinghAgrw.M.dist.m$FlyBaseID %in% subset.sss,]$percent.dissim, na.rm = T)
Male.MCsim.SinghAgrw.M.dist.f <- compareSplicingProfiles(male.exp_MCsim, fem.exp_SinghAgrw)
mean(Male.MCsim.SinghAgrw.M.dist.f$percent.dissim, na.rm = T)
mean(Male.MCsim.SinghAgrw.M.dist.f[Male.MCsim.SinghAgrw.M.dist.f$FlyBaseID %in% subset.sss,]$percent.dissim, na.rm = T)

Fem.MCsim.SinghAgrw.M.dist.m <- compareSplicingProfiles(fem.exp_MCsim, male.exp_SinghAgrw)
mean(Fem.MCsim.SinghAgrw.M.dist.m$percent.dissim, na.rm = T)
mean(Fem.MCsim.SinghAgrw.M.dist.m[Fem.MCsim.SinghAgrw.M.dist.m$FlyBaseID %in% subset.sss,]$percent.dissim, na.rm = T)
Fem.MCsim.SinghAgrw.M.dist.f <- compareSplicingProfiles(fem.exp_MCsim, fem.exp_SinghAgrw)
mean(Fem.MCsim.SinghAgrw.M.dist.f$percent.dissim, na.rm = T)
mean(Fem.MCsim.SinghAgrw.M.dist.f[Fem.MCsim.SinghAgrw.M.dist.f$FlyBaseID %in% subset.sss,]$percent.dissim, na.rm = T)



#######



# Get a subset of SSS genes that we can look at
#######
sss.all <- ASE.sig.SSS[ASE.sig.SSS %in% SinghAgrawal.sig.SSS &
                         ASE.sig.SSS %in% filter.low.exp.genes.q25]
length(sss.all)

# dot plots to check how many points fall into the "irregular" category, but for the MC populations
# "irregular" meaning that when looking at males (females), 
# the points fall closer to the reference females (males)
test <- merge(Male.ASE.abs.M.dist.m, Male.ASE.abs.M.dist.f, by = "FlyBaseID")
plot_corr(test[test$FlyBaseID %in% sss.all, ], 
          x = "percent.dissim.x", y = "percent.dissim.y", colx = "black", coly = "black",colNonCon = "black",xlab = "to males",
          ylab = "to females", lim = 1, title = "") + coord_cartesian(xlim = c(0,1), ylim = c(0,1))

# get a list of all irregular genes from each MC population treatment
fem.sim <- test[test$FlyBaseID %in% sss.all & test$percent.dissim.x < test$percent.dissim.y, ]$FlyBaseID
male.sim <- test[test$FlyBaseID %in% sss.all & test$percent.dissim.x > test$percent.dissim.y, ]$FlyBaseID

fem.com <- test[test$FlyBaseID %in% sss.all & test$percent.dissim.x < test$percent.dissim.y, ]$FlyBaseID
male.com <- test[test$FlyBaseID %in% sss.all & test$percent.dissim.x > test$percent.dissim.y, ]$FlyBaseID

fem.abs <- test[test$FlyBaseID %in% sss.all & test$percent.dissim.x < test$percent.dissim.y, ]$FlyBaseID
male.abs <- test[test$FlyBaseID %in% sss.all & test$percent.dissim.x > test$percent.dissim.y, ]$FlyBaseID


# interesting: the female ones have more points that fall in the "irregular" category
length(fem.sim)
length(male.sim)
length(fem.abs)
length(male.abs)
length(fem.com)
length(male.com)

# list of all "irregular" genes to be excluded
all.exclude <- unique(c(fem.sim, fem.abs, fem.com, male.sim, male.abs, male.com)) 

# subset of genes that are more consistently dimorphic:
# 1) SSS in both ASE and OSADA populations
# 2) fall into the right "sex-profile" (i.e., not "irregular") when seen in the MC populations
subset.sss <- (sss.all[!sss.all %in% all.exclude])



# distance from SSAV females to SinghAgrawal females
ASE.SinghAgrw.F.dist <- compareSplicingProfiles(fem.exp_ASE, fem.exp_SinghAgrw)
mean(ASE.SinghAgrw.F.dist$percent.sim, na.rm=T) 

# ASE vs Osada populations
# ASE males to Osada females
Male.ASE.com.M.dist.f <- compareSplicingProfiles(male.exp_ASE, fem.exp_MCcom)
# ASE males to Osada males
Male.ASE.com.M.dist.m <- compareSplicingProfiles(male.exp_ASE, male.exp_MCcom)
# ASE females to Osada females
Fem.ASE.com.M.dist.f <- compareSplicingProfiles(fem.exp_ASE, fem.exp_MCcom)
# ASE females to Osada males
Fem.ASE.com.M.dist.m <- compareSplicingProfiles(fem.exp_ASE, male.exp_MCcom) 

SinghAgrw.sim.compare.male <- CompareDistances(Male.SinghAgrw.sim.M.dist.m, Male.SinghAgrw.sim.M.dist.f )
SinghAgrw.sim.compare.fem <- CompareDistances(Fem.SinghAgrw.sim.M.dist.m, Fem.SinghAgrw.sim.M.dist.f)
splicing.MF.metric.plot(RedData = SinghAgrw.sim.compare.male[SinghAgrw.sim.compare.male$FlyBaseID %in% sss.all,], 
                        NRData = SinghAgrw.sim.compare.fem[SinghAgrw.sim.compare.fem$FlyBaseID %in% sss.all,], 
                        plotCol = "M.dis.F",
                        colour_red = "steelblue3", colour_NR = "magenta4")


# distance from SSAV males to ASE males
A.ASE.M.dist <- compareSplicingProfiles(male.exp_ASE, male.exp)
mean(A.ASE.M.dist$percent.sim, na.rm = T)
# distance from SSAV females to ASE females
A.ASE.F.dist <- compareSplicingProfiles(fem.exp_ASE, fem.exp)
mean(A.ASE.F.dist$percent.sim, na.rm=T)  # females further from the ASE?


# distance from SSAV males to SinghAgrawal males
A.SinghAgrw.M.dist <- compareSplicingProfiles(male.exp_SinghAgrw, male.exp)
mean(A.SinghAgrw.M.dist$percent.sim, na.rm = T)
# distance from SSAV females to SinghAgrawal females
A.SinghAgrw.F.dist <- compareSplicingProfiles(fem.exp_SinghAgrw, fem.exp)
mean(A.SinghAgrw.F.dist$percent.sim, na.rm=T)







# Osada vs ASE dimorphism?
# more dimorphic = M/F percent similarity smaller
# Osada more dimorphic = 
test.count <- test %>% 
  dplyr::count(X.more.dimorphic = percent.sim.x < percent.sim.y,
               same = percent.sim.x == percent.sim.y) # on the top side of plot?
# 2226 have the same M/F perc similarity
# 7192 ASE is less similar between M/F (more dimorphic)
# 4415 Osada less similar between M/F (more dimorphic) so no.

# ASE more dimorphic than Osada
# SSAV more dimorphic than Osada
# SSAV less dimorphic than ASE


ggplot() +
geom_histogram(data=ASE.MvF.dist[ASE.MvF.dist$FlyBaseID %in% ASE.sig.SSS & 
                                   ASE.MvF.dist$FlyBaseID %in% filter.low.exp.genes.q25,], 
               aes(x=dist, y= ..density.. *100),  fill = "red", bins = 40, alpha = 0.5) +
  geom_histogram(data=ASE.MvF.dist[ASE.MvF.dist$FlyBaseID %in% ASE.nonsig.SSS & 
                                     ASE.MvF.dist$FlyBaseID %in% filter.low.exp.genes.q25,],
                 aes(x=dist, y= ..density.. *100), fill = "blue", alpha =0.5, bins =40)




#######



# list of genes separated by quantiles of dimorphism (euc distance)
######
ASE.sss.only <- ASE.MvF.dist[ASE.MvF.dist$FlyBaseID %in% ASE.sig.SSS,]
ASE.MvF.dist_q25 <- ASE.sss.only$FlyBaseID[!is.na(ASE.sss.only$percent.sim) & 
                                             ASE.sss.only$dist <= quantile(ASE.sss.only$dist, 0.25, na.rm=T)]
ASE.MvF.dist_q50 <- ASE.sss.only$FlyBaseID[!is.na(ASE.sss.only$percent.sim) & 
                                             ASE.sss.only$dist > quantile(ASE.sss.only$dist, 0.25, na.rm=T) &
                                             ASE.sss.only$dist <= quantile(ASE.sss.only$dist, 0.50, na.rm=T)]
ASE.MvF.dist_q75 <- ASE.sss.only$FlyBaseID[!is.na(ASE.sss.only$percent.sim) & 
                                             ASE.sss.only$dist > quantile(ASE.sss.only$dist, 0.50, na.rm=T) &
                                             ASE.sss.only$dist <= quantile(ASE.sss.only$dist, 0.75, na.rm=T)]
ASE.MvF.dist_q100 <- ASE.sss.only$FlyBaseID[!is.na(ASE.sss.only$percent.sim) & 
                                              ASE.sss.only$dist >= quantile(ASE.sss.only$dist, 0.75, na.rm=T)]

SinghAgrawal.sss.only <- SinghAgrawal.MvF.dist[SinghAgrawal.MvF.dist$FlyBaseID %in% SinghAgrawal.sig.SSS,]
SinghAgrawal.MvF.dist_q25 <- SinghAgrawal.sss.only$FlyBaseID[!is.na(SinghAgrawal.sss.only$percent.sim) & 
                                                               SinghAgrawal.sss.only$dist <= quantile(SinghAgrawal.sss.only$dist, 0.25, na.rm=T)]
SinghAgrawal.MvF.dist_q50 <- SinghAgrawal.sss.only$FlyBaseID[!is.na(SinghAgrawal.sss.only$percent.sim) & 
                                                               SinghAgrawal.sss.only$dist > quantile(SinghAgrawal.sss.only$dist, 0.25, na.rm=T) &
                                                               SinghAgrawal.sss.only$dist <= quantile(SinghAgrawal.sss.only$dist, 0.50, na.rm=T)]
SinghAgrawal.MvF.dist_q75 <- SinghAgrawal.sss.only$FlyBaseID[!is.na(SinghAgrawal.sss.only$percent.sim) & 
                                                               SinghAgrawal.sss.only$dist > quantile(SinghAgrawal.sss.only$dist, 0.50, na.rm=T) &
                                                               SinghAgrawal.sss.only$dist <= quantile(SinghAgrawal.sss.only$dist, 0.75, na.rm=T)]
SinghAgrawal.MvF.dist_q100 <- SinghAgrawal.sss.only$FlyBaseID[!is.na(SinghAgrawal.sss.only$percent.sim) & 
                                                                SinghAgrawal.sss.only$dist >= quantile(SinghAgrawal.sss.only$dist, 0.75, na.rm=T)]


######

# list of genes separated by quantiles of dimorphism (% similarity)
######
ASE.sss.only <- ASE.MvF.dist[ASE.MvF.dist$FlyBaseID %in% ASE.sig.SSS,]
ASE.MvF.perc.sim_q25 <- ASE.sss.only$FlyBaseID[!is.na(ASE.sss.only$percent.sim) &
                                                 ASE.sss.only$percent.sim < quantile(ASE.sss.only$percent.sim, 0.25, na.rm=T)]
ASE.MvF.perc.sim_q50 <- ASE.sss.only$FlyBaseID[!is.na(ASE.sss.only$percent.sim) &
                                                 ASE.sss.only$percent.sim > quantile(ASE.sss.only$percent.sim, 0.25, na.rm=T) &
                                                 ASE.sss.only$percent.sim < quantile(ASE.sss.only$percent.sim, 0.50, na.rm=T)]
ASE.MvF.perc.sim_q75 <- ASE.sss.only$FlyBaseID[!is.na(ASE.sss.only$percent.sim) &
                                                 ASE.sss.only$percent.sim > quantile(ASE.sss.only$percent.sim, 0.50, na.rm=T) &
                                               ASE.sss.only$percent.sim < quantile(ASE.sss.only$percent.sim, 0.75, na.rm=T)]
ASE.MvF.perc.sim_q100 <- ASE.sss.only$FlyBaseID[!is.na(ASE.sss.only$percent.sim) &
                                                  ASE.sss.only$percent.sim > quantile(ASE.sss.only$percent.sim, 0.75, na.rm=T)]



# ASE.MvF.perc.sim_q25 <- ASE.MvF.dist$FlyBaseID[!is.na(ASE.MvF.dist$percent.sim) &
#                                              ASE.MvF.dist$percent.sim < quantile(ASE.MvF.dist$percent.sim, 0.25, na.rm=T)]
# ASE.MvF.perc.sim_q50 <- ASE.MvF.dist$FlyBaseID[!is.na(ASE.MvF.dist$percent.sim) & 
#                                              ASE.MvF.dist$percent.sim > quantile(ASE.MvF.dist$percent.sim, 0.25, na.rm=T) &
#                                              ASE.MvF.dist$percent.sim < quantile(ASE.MvF.dist$percent.sim, 0.50, na.rm=T)]
# ASE.MvF.perc.sim_q75 <- ASE.MvF.dist$FlyBaseID[!is.na(ASE.MvF.dist$percent.sim) & 
#                                              ASE.MvF.dist$percent.sim > quantile(ASE.MvF.dist$percent.sim, 0.50, na.rm=T) &
#                                              ASE.MvF.dist$percent.sim < quantile(ASE.MvF.dist$percent.sim, 0.75, na.rm=T)]
# ASE.MvF.perc.sim_q100 <- ASE.MvF.dist$FlyBaseID[!is.na(ASE.MvF.dist$percent.sim) & 
#                                               ASE.MvF.dist$percent.sim > quantile(ASE.MvF.dist$percent.sim, 0.75, na.rm=T)]



SinghAgrawal.sss.only <- SinghAgrawal.MvF.dist[SinghAgrawal.MvF.dist$FlyBaseID %in% SinghAgrawal.sig.SSS,]
SinghAgrawal.MvF.perc.sim_q25 <- SinghAgrawal.sss.only$FlyBaseID[!is.na(SinghAgrawal.sss.only$percent.sim) &
                                                                   SinghAgrawal.sss.only$percent.sim <= quantile(SinghAgrawal.sss.only$percent.sim, 0.25, na.rm=T)]
SinghAgrawal.MvF.perc.sim_q50 <- SinghAgrawal.sss.only$FlyBaseID[!is.na(SinghAgrawal.sss.only$percent.sim) &
                                                                   SinghAgrawal.sss.only$percent.sim > quantile(SinghAgrawal.sss.only$percent.sim, 0.25, na.rm=T) &
                                                                 SinghAgrawal.sss.only$percent.sim <= quantile(SinghAgrawal.sss.only$percent.sim, 0.50, na.rm=T)]
SinghAgrawal.MvF.perc.sim_q75 <- SinghAgrawal.sss.only$FlyBaseID[!is.na(SinghAgrawal.sss.only$percent.sim) &
                                                                   SinghAgrawal.sss.only$percent.sim > quantile(SinghAgrawal.sss.only$percent.sim, 0.50, na.rm=T) &
                                                                   SinghAgrawal.sss.only$percent.sim <= quantile(SinghAgrawal.sss.only$percent.sim, 0.75, na.rm=T)]
SinghAgrawal.MvF.perc.sim_q100 <- SinghAgrawal.sss.only$FlyBaseID[!is.na(SinghAgrawal.sss.only$percent.sim) &
                                                                    SinghAgrawal.sss.only$percent.sim >= quantile(SinghAgrawal.sss.only$percent.sim, 0.75, na.rm=T)]


# SinghAgrawal.MvF.perc.sim_q25 <- SinghAgrawal.MvF.dist$FlyBaseID[!is.na(SinghAgrawal.MvF.dist$percent.sim) & 
#                                                                SinghAgrawal.MvF.dist$percent.sim <= quantile(SinghAgrawal.MvF.dist$percent.sim, 0.25, na.rm=T)]
# SinghAgrawal.MvF.perc.sim_q50 <- SinghAgrawal.MvF.dist$FlyBaseID[!is.na(SinghAgrawal.MvF.dist$percent.sim) & 
#                                                                SinghAgrawal.MvF.dist$percent.sim > quantile(SinghAgrawal.MvF.dist$percent.sim, 0.25, na.rm=T) &
#                                                                SinghAgrawal.MvF.dist$percent.sim <= quantile(SinghAgrawal.MvF.dist$percent.sim, 0.50, na.rm=T)]
# SinghAgrawal.MvF.perc.sim_q75 <- SinghAgrawal.MvF.dist$FlyBaseID[!is.na(SinghAgrawal.MvF.dist$percent.sim) & 
#                                                                SinghAgrawal.MvF.dist$percent.sim > quantile(SinghAgrawal.MvF.dist$percent.sim, 0.50, na.rm=T) &
#                                                                SinghAgrawal.MvF.dist$percent.sim <= quantile(SinghAgrawal.MvF.dist$percent.sim, 0.75, na.rm=T)]
# SinghAgrawal.MvF.perc.sim_q100 <- SinghAgrawal.MvF.dist$FlyBaseID[!is.na(SinghAgrawal.MvF.dist$percent.sim) & 
#                                                                 SinghAgrawal.MvF.dist$percent.sim >= quantile(SinghAgrawal.MvF.dist$percent.sim, 0.75, na.rm=T)]
######


# -------- male comparison -------------
########
# expression of Red A males
A.m.Red.v.Males <- compareSplicingProfiles(A.m.Red.norm.exp, male.exp)
A.m.Red.v.Males_ase <- compareSplicingProfiles(A.m.Red.norm.exp, male.exp_ASE)
A.m.Red.v.Males_SinghAgrw <- compareSplicingProfiles(A.m.Red.norm.exp, male.exp_SinghAgrw)
A.m.Red.v.Males_MCsim <- compareSplicingProfiles(A.m.Red.norm.exp, male.exp_MCsim)
# expression of NR A males
A.m.Red.v.Fem <- compareSplicingProfiles(A.m.Red.norm.exp, fem.exp)
A.m.Red.v.Fem_ase <- compareSplicingProfiles(A.m.Red.norm.exp, fem.exp_ASE)
A.m.Red.v.Fem_SinghAgrw <- compareSplicingProfiles(A.m.Red.norm.exp, fem.exp_SinghAgrw)
A.m.Red.v.Fem_MCsim <- compareSplicingProfiles(A.m.Red.norm.exp, fem.exp_MCsim)

# join both distances
A.m.Red.compare <- CompareDistances(A.m.Red.v.Males, A.m.Red.v.Fem)
A.m.Red.compare_ase <- CompareDistances(A.m.Red.v.Males_ase, A.m.Red.v.Fem_ase) 
A.m.Red.compare_SinghAgrw <- CompareDistances(A.m.Red.v.Males_SinghAgrw, A.m.Red.v.Fem_SinghAgrw)
A.m.Red.compare_MCsim <- CompareDistances(A.m.Red.v.Males_MCsim, A.m.Red.v.Fem_MCsim)

# ---
# expression of NR A males
A.m.NR.v.Males <- compareSplicingProfiles(A.m.NR.norm.exp, male.exp)
A.m.NR.v.Males_ase <- compareSplicingProfiles(A.m.NR.norm.exp, male.exp_ASE) # with ASE data
A.m.NR.v.Males_SinghAgrw <- compareSplicingProfiles(A.m.NR.norm.exp, male.exp_SinghAgrw)
A.m.NR.v.Males_MCsim <- compareSplicingProfiles(A.m.NR.norm.exp, male.exp_MCsim)

A.m.NR.v.Fem <- compareSplicingProfiles(A.m.NR.norm.exp, fem.exp)
A.m.NR.v.Fem_ase <- compareSplicingProfiles(A.m.NR.norm.exp, fem.exp_ASE)
A.m.NR.v.Fem_SinghAgrw <- compareSplicingProfiles(A.m.NR.norm.exp, fem.exp_SinghAgrw)
A.m.NR.v.Fem_MCsim <- compareSplicingProfiles(A.m.NR.norm.exp, fem.exp_MCsim)

# join both distances
A.m.NR.compare <- CompareDistances(A.m.NR.v.Males, A.m.NR.v.Fem) 
A.m.NR.compare_ase <- CompareDistances(A.m.NR.v.Males_ase, A.m.NR.v.Fem_ase)
A.m.NR.compare_SinghAgrw <- CompareDistances(A.m.NR.v.Males_SinghAgrw, A.m.NR.v.Fem_SinghAgrw)
A.m.NR.compare_MCsim <- CompareDistances(A.m.NR.v.Males_MCsim, A.m.NR.v.Fem_MCsim)
#######

# A males figures
#######
A.m.splice
A.m.splice_ase

A.m.splice_ase.SSS
A.m.splice_ase.not.SSS

A.m.splice_ase.q25
A.m.splice_ase.q50
A.m.splice_ase.q75
A.m.splice_ase.q100

A.m.splice_ase.q25.filtered
A.m.splice_ase.q50.filtered
A.m.splice_ase.q75.filtered
A.m.splice_ase.q100.filtered


#######


# -------- female comparison -------------
########
# expression of Red A females
A.f.Red.v.Males <- compareSplicingProfiles(A.f.Red.norm.exp, male.exp)
A.f.Red.v.Males_ase <- compareSplicingProfiles(A.f.Red.norm.exp, male.exp_ASE)
A.f.Red.v.Males_SinghAgrw <- compareSplicingProfiles(A.f.Red.norm.exp, male.exp_SinghAgrw)
A.f.Red.v.Males_MCsim <- compareSplicingProfiles(A.f.Red.norm.exp, male.exp_MCsim)

A.f.Red.v.Fem <- compareSplicingProfiles(A.f.Red.norm.exp, fem.exp)
A.f.Red.v.Fem_ase <- compareSplicingProfiles(A.f.Red.norm.exp, fem.exp_ASE)
A.f.Red.v.Fem_SinghAgrw <- compareSplicingProfiles(A.f.Red.norm.exp, fem.exp_SinghAgrw)
A.f.Red.v.Fem_MCsim <- compareSplicingProfiles(A.f.Red.norm.exp, fem.exp_MCsim)

# join both distances
A.f.Red.compare <- CompareDistances(A.f.Red.v.Males, A.f.Red.v.Fem)
A.f.Red.compare_ase <- CompareDistances(A.f.Red.v.Males_ase, A.f.Red.v.Fem_ase)
A.f.Red.compare_SinghAgrw <- CompareDistances(A.f.Red.v.Males_SinghAgrw, A.f.Red.v.Fem_SinghAgrw)
A.f.Red.compare_MCsim <- CompareDistances(A.f.Red.v.Males_MCsim, A.f.Red.v.Fem_MCsim)

# expression of NR A females
A.f.NR.v.Males <- compareSplicingProfiles(A.f.NR.norm.exp, male.exp)
A.f.NR.v.Males_ase <- compareSplicingProfiles(A.f.NR.norm.exp, male.exp_ASE)
A.f.NR.v.Males_SinghAgrw <- compareSplicingProfiles(A.f.NR.norm.exp, male.exp_SinghAgrw)
A.f.NR.v.Males_MCsim <- compareSplicingProfiles(A.f.NR.norm.exp, male.exp_MCsim)

A.f.NR.v.Fem <- compareSplicingProfiles(A.f.NR.norm.exp, fem.exp)
A.f.NR.v.Fem_ase <- compareSplicingProfiles(A.f.NR.norm.exp, fem.exp_ASE)
A.f.NR.v.Fem_SinghAgrw <- compareSplicingProfiles(A.f.NR.norm.exp, fem.exp_SinghAgrw)
A.f.NR.v.Fem_MCsim <- compareSplicingProfiles(A.f.NR.norm.exp, fem.exp_MCsim)

A.f.NR.compare <- CompareDistances(A.f.NR.v.Males, A.f.NR.v.Fem)
A.f.NR.compare_ase <- CompareDistances(A.f.NR.v.Males_ase, A.f.NR.v.Fem_ase)
A.f.NR.compare_SinghAgrw <- CompareDistances(A.f.NR.v.Males_SinghAgrw, A.f.NR.v.Fem_SinghAgrw)
A.f.NR.compare_MCsim <- CompareDistances(A.f.NR.v.Males_MCsim, A.f.NR.v.Fem_MCsim)
#######

# A.f figures
######
A.f.splice
A.f.splice_ase 

A.f.splice_ase.SSS
A.f.splice_ase.not.SSS

A.f.splice_ase.q25
A.f.splice_ase.q50
A.f.splice_ase.q75
A.f.splice_ase.q100


A.f.splice_ase.q25.filtered
A.f.splice_ase.q50.filtered
A.f.splice_ase.q75.filtered
A.f.splice_ase.q100.filtered

######

# -------- control male comparison -------------
#######
C.m.Red.v.Males <- compareSplicingProfiles(C.m.Red.norm.exp, male.exp)
C.m.Red.v.Males_ase <- compareSplicingProfiles(C.m.Red.norm.exp, male.exp_ASE)
C.m.Red.v.Males_SinghAgrw <- compareSplicingProfiles(C.m.Red.norm.exp, male.exp_SinghAgrw)
C.m.Red.v.Males_MCsim <- compareSplicingProfiles(C.m.Red.norm.exp, male.exp_MCsim)

C.m.Red.v.Fem <- compareSplicingProfiles(C.m.Red.norm.exp, fem.exp)
C.m.Red.v.Fem_ase <- compareSplicingProfiles(C.m.Red.norm.exp, fem.exp_ASE)
C.m.Red.v.Fem_SinghAgrw <- compareSplicingProfiles(C.m.Red.norm.exp, fem.exp_SinghAgrw)
C.m.Red.v.Fem_MCsim <- compareSplicingProfiles(C.m.Red.norm.exp, fem.exp_MCsim)

C.m.Red.compare <- CompareDistances(C.m.Red.v.Males, C.m.Red.v.Fem)
C.m.Red.compare_ase <- CompareDistances(C.m.Red.v.Males_ase, C.m.Red.v.Fem_ase)
C.m.Red.compare_SinghAgrw <- CompareDistances(C.m.Red.v.Males_SinghAgrw, C.m.Red.v.Fem_SinghAgrw)
C.m.Red.compare_MCsim <- CompareDistances(C.m.Red.v.Males_ase, C.m.Red.v.Fem_MCsim)

C.m.NR.v.Males <- compareSplicingProfiles(C.m.NR.norm.exp, male.exp)
C.m.NR.v.Males_ase <- compareSplicingProfiles(C.m.NR.norm.exp, male.exp_ASE)
C.m.NR.v.Males_SinghAgrw <- compareSplicingProfiles(C.m.NR.norm.exp, male.exp_SinghAgrw)
C.m.NR.v.Males_MCsim <- compareSplicingProfiles(C.m.NR.norm.exp, male.exp_MCsim)

C.m.NR.v.Fem <- compareSplicingProfiles(C.m.NR.norm.exp, fem.exp)
C.m.NR.v.Fem_ase <- compareSplicingProfiles(C.m.NR.norm.exp, fem.exp_ASE)
C.m.NR.v.Fem_SinghAgrw <- compareSplicingProfiles(C.m.NR.norm.exp, fem.exp_SinghAgrw)
C.m.NR.v.Fem_MCsim <- compareSplicingProfiles(C.m.NR.norm.exp, fem.exp_MCsim)

C.m.NR.compare <- CompareDistances(C.m.NR.v.Males, C.m.NR.v.Fem)
C.m.NR.compare_ase <- CompareDistances(C.m.NR.v.Males_ase, C.m.NR.v.Fem_ase)
C.m.NR.compare_SinghAgrw <- CompareDistances(C.m.NR.v.Males_SinghAgrw, C.m.NR.v.Fem_SinghAgrw)
C.m.NR.compare_MCsim<- CompareDistances(C.m.NR.v.Males_MCsim, C.m.NR.v.Fem_MCsim)
#######

# C.m figures
#######
C.m.splice
C.m.splice_ase

C.m.splice_ase.SSS
C.m.splice_ase.not.SSS

C.m.splice_ase.q25
C.m.splice_ase.q50
C.m.splice_ase.q75
C.m.splice_ase.q100


C.m.splice_ase.q25.filtered
C.m.splice_ase.q50.filtered
C.m.splice_ase.q75.filtered
C.m.splice_ase.q100.filtered



########


Red.Af.NR.Af <- compareSplicingProfiles(A.f.Red.norm.exp, A.f.NR.norm.exp)
Red.Af.NR.Am <- compareSplicingProfiles(A.f.Red.norm.exp, A.m.NR.norm.exp)
Red.Af.Red.Am <- compareSplicingProfiles(A.f.Red.norm.exp, A.m.Red.norm.exp)
Red.Am.NR.Af <- compareSplicingProfiles(A.m.Red.norm.exp, A.f.NR.norm.exp)
Red.Am.NR.Am <- compareSplicingProfiles(A.m.Red.norm.exp, A.m.NR.norm.exp)
NR.Af.NR.Am <- compareSplicingProfiles(A.f.NR.norm.exp, A.m.NR.norm.exp)
Red.Am.Red.Af <- compareSplicingProfiles(A.m.Red.norm.exp, A.f.Red.norm.exp)
NR.Am.Red.Af <- compareSplicingProfiles(A.m.NR.norm.exp, A.f.Red.norm.exp)
NR.Am.NR.Af <- compareSplicingProfiles(A.m.NR.norm.exp, A.f.NR.norm.exp)

mean(NR.Am.Red.Af$percent.dissim, na.rm = T)
mean(Red.Am.NR.Am[Red.Am.NR.Am$FlyBaseID %in% ASE.sig.SSS,]$percent.dissim, na.rm = T)
mean(Red.Af.NR.Af[Red.Af.NR.Af$FlyBaseID %in% ASE.sig.SSS,]$percent.dissim, na.rm = T)


# lists for different subsets of genes that can be used to stratify downstream analyses
#######
# genes with significant differences in splicing
# only genes that could be tested in both A males and A females
All.geno.Sig.Spliced.genes <- unique(jseq.All.geno[!is.na(jseq.All.geno$sig.hit) & jseq.All.geno$sig.hit,]$FlyBaseID)
All.geno.not.Sig.Spliced.genes <- unique(jseq.All.geno[!is.na(jseq.All.geno$sig.hit) & !jseq.All.geno$sig.hit,]$FlyBaseID)

length(All.geno.Sig.Spliced.genes)
length(All.geno.not.Sig.Spliced.genes)

length(filter.low.exp.genes.q10) # list of genes that passed the 10% FPKM filter according to ASE data
length(filter.low.exp.genes.q25) # list of genes that passed the 25% FPKM filter according to ASE data

# only consider genes with significant SSS in ASE data
length(ASE.sig.SSS)
length(ASE.nonsig.SSS)
# get genes with significant SSS in Singh and Agrawal paper
SinghAgrawalSDIU <- read.csv(file="~/Desktop/UofT/SSAV_RNA/Data/SBGEandSSSdataForMBE.csv", sep=",", header=TRUE)
colnames(SinghAgrawalSDIU)[2] <- "FlyBaseID"
SinghAgrawal.sig.SSS <- SinghAgrawalSDIU$FlyBaseID[!is.na(SinghAgrawalSDIU$SDIU.body.sig) & SinghAgrawalSDIU$SDIU.body.sig]
SinghAgrawal.nonsig.SSS <- SinghAgrawalSDIU$FlyBaseID[!is.na(SinghAgrawalSDIU$SDIU.body.sig) & !SinghAgrawalSDIU$SDIU.body.sig]
length(SinghAgrawal.sig.SSS)
length(SinghAgrawal.nonsig.SSS)
rm(SinghAgrawalSDIU)


# quantiles of dimorphism (euc distance)
ASE.MvF.dist_q25
ASE.MvF.dist_q50
ASE.MvF.dist_q75
ASE.MvF.dist_q100

SinghAgrawal.MvF.dist_q25
SinghAgrawal.MvF.dist_q50
SinghAgrawal.MvF.dist_q75
SinghAgrawal.MvF.dist_q100

# quantiles of dimorphism (% similarity)
ASE.MvF.perc.sim_q25
ASE.MvF.perc.sim_q50
ASE.MvF.perc.sim_q75
ASE.MvF.perc.sim_q100

SinghAgrawal.MvF.perc.sim_q25
SinghAgrawal.MvF.perc.sim_q50
SinghAgrawal.MvF.perc.sim_q75
SinghAgrawal.MvF.perc.sim_q100
#######


# A bunch of t-tests ...
######
# Compare Red vs NR 
# to males
sampleTypes <- c("A.m", "A.f", "C.m")
t.test.table <- data.frame(sampleType = c("A.m", "A.m", "A.f", "A.f", "C.m", "C.m"),
                           againstASE = rep(c("male", "female"), 3))
p.val.list <- NULL
avg.Red <- NULL
avg.NR <- NULL
diff <- NULL

for(i in sampleTypes) {
  Fem_RedDat <- get(paste0(i,".Red.v.Fem_ase")) # .Red.v.Fem_SinghAgrw
  Fem_NRDat <- get(paste0(i,".NR.v.Fem_ase")) # .NR.v.Fem_SinghAgrw
  
  Male_RedDat <- get(paste0(i,".Red.v.Males_ase")) # .Red.v.Males_SinghAgrw
  Male_NRDat <- get(paste0(i,".NR.v.Males_ase")) # .NR.v.Males_SinghAgrw
  
  Fem_RedDat <- Fem_RedDat[Fem_RedDat$FlyBaseID %in% A.f.sig.AS,] # SinghAgrawal.sig.SSS
  Fem_NRDat <- Fem_NRDat[Fem_NRDat$FlyBaseID %in% A.f.sig.AS,]
  

  Male_RedDat <- Male_RedDat[Male_RedDat$FlyBaseID %in% A.f.sig.AS,]
  Male_NRDat <- Male_NRDat[Male_NRDat$FlyBaseID %in% A.f.sig.AS,]
  
  print(length(!is.na(Fem_RedDat$percent.dissim)))
  print(length(!is.na(Fem_NRDat$percent.dissim)))
  print(length(!is.na(Male_RedDat$percent.dissim)))
  print(length(!is.na(Male_NRDat$percent.dissim)))
  

  p.val.list <- c(p.val.list, t.test(Male_RedDat$percent.dissim, 
                                     Male_NRDat$percent.dissim, paired = T)$p.value)
  print(t.test(Male_RedDat$percent.dissim, 
               Male_NRDat$percent.dissim, paired = T))
  p.val.list <- c(p.val.list, t.test(Fem_RedDat$percent.dissim, 
                                     Fem_NRDat$percent.dissim, paired = T)$p.value)
  print(t.test(Fem_RedDat$percent.dissim, 
               Fem_NRDat$percent.dissim, paired = T))
  
  avg.Red <- c(avg.Red, mean(Male_RedDat$percent.dissim, na.rm = T))
  avg.Red <- c(avg.Red, mean(Fem_RedDat$percent.dissim, na.rm = T))
  
  avg.NR <- c(avg.NR, mean(Male_NRDat$percent.dissim, na.rm = T))
  avg.NR <- c(avg.NR, mean(Fem_NRDat$percent.dissim, na.rm = T))
  
  diff <- c(diff, t.test(Male_RedDat$percent.dissim, 
                         Male_NRDat$percent.dissim, paired = T)$estimate)
  diff <- c(diff, t.test(Fem_RedDat$percent.dissim, 
                         Fem_NRDat$percent.dissim, paired = T)$estimate)
  
}

t.test.table$pval = p.val.list
t.test.table$avg.Red = avg.Red
t.test.table$avg.NR = avg.NR
t.test.table$diff = diff

write.table(t.test.table, file = "Results/tmp.csv", sep = ",", quote = FALSE, row.names = F)

######


# figures and t-tests comparing (M-F)/(M+F) metric
######
splicing.MF.metric.plot <- function(RedData, NRData, 
                                    colour_red = "red3", colour_NR = "grey20", 
                                    plotCol){
  splice.plot <- ggplot() + 
    geom_histogram(data= NRData, 
                   aes_string(x=plotCol), color = colour_NR, fill = colour_NR, alpha = 0.5, bins = 45) +
    geom_histogram(data= RedData, 
                   aes_string(x=plotCol), color = colour_red, fill = colour_red, alpha = 0.5, bins = 45) +
    labs(x = "female <<---------------------->> male") +
    geom_vline(xintercept = 0, color = "black", linewidth = 1) +
    # coord_cartesian(xlim=c(-0.5,0.5)) +
    theme_classic() +
    theme(plot.title.position = c("panel"),
          legend.title = element_blank(),
          legend.position = c("bottom"),
          legend.box.background = element_rect(),
          legend.text = element_text(size = 20, color = "black"),
          axis.text.x = element_text(size=20, margin = margin(5,0,0,0), color = "black"),
          axis.text.y = element_text(size=20, margin = margin(0,5,0,0), color = "black"),
          axis.title.x = element_text(size=30, margin = margin(10,0,0,0), color = "black"),
          axis.title.y = element_text(size=30, margin = margin(0,10,0,0), color = "black"),
          plot.title = element_text(size=40, margin = margin(0,0,0,0), color = "black"),
          plot.margin = margin(6,6,6,6))
  return(splice.plot)
}


splicing.MFdiff.plot <- function(RedData, NRData, plotCol){
  tmp <- data.frame(FlyBaseID = RedData$FlyBaseID,
                    diff = RedData[[plotCol]] - NRData[[plotCol]])
  
  splice.plot <- ggplot() + 
    geom_histogram(data= tmp, 
                   aes(x= diff), color = "grey20", fill = "grey20", alpha = 0.5, bins = 45) +
    labs(x = "female <<---------------------->> male") +
    geom_vline(xintercept = 0, color = "black", linewidth = 1) +
    # coord_cartesian(xlim=c(-0.5,0.5)) +
    theme_classic() +
    theme(plot.title.position = c("panel"),
          legend.title = element_blank(),
          legend.position = c("bottom"),
          legend.box.background = element_rect(),
          legend.text = element_text(size = 20, color = "black"),
          axis.text.x = element_text(size=20, margin = margin(5,0,0,0), color = "black"),
          axis.text.y = element_text(size=20, margin = margin(0,5,0,0), color = "black"),
          axis.title.x = element_text(size=30, margin = margin(10,0,0,0), color = "black"),
          axis.title.y = element_text(size=30, margin = margin(0,10,0,0), color = "black"),
          plot.title = element_text(size=40, margin = margin(0,0,0,0), color = "black"),
          plot.margin = margin(6,6,6,6))
  return(splice.plot)
}

# the two functions below are used in the corr plot function below
# quad_count: counts the number of point in each of the 8 quadrants
# variables: dat = data.frame object containing the values to be plotted
#            x = values for group 1, plotted on the x-axis
#            y = values for group 2, plotted on the y-axis
#            lim = for plotting code. specificy the x-lim and y-lim of the plot
# return: a data frame with number of points, 
#         proportion relative to total points, 
#         and plotting coordinates
quad_count <- function(dat, x, y, lim = 5){
  count <- dat %>%
    # Count how many with each combination of X and Y being positive
    dplyr::count(na.rm = T,
                 right = .[[x]] > 0, # on the right side of plot?
                 top = .[[y]] > 0, # on the top side of plot?
                 UP = !top & abs(.[[x]]) > abs(.[[y]]) | # quadrant III up-left or quadrant IV up-right
                 top & abs(.[[x]]) < abs(.[[y]])) %>%  # quadrant I up-right or quadrant II up-left
    mutate(perc = n/sum(n, na.rm = T)) %>% # calculate percentage of points relative to total number of points
    
    # this is another strange one for setting up the coordinates
    mutate(conc = right & top | (!right & !top), # quadrant I and quadrant III (the concordant changes)
           dir_UP = conc & UP, # concordant changes where x > y
           dir_DOWN = conc & !UP) %>% # concordant changes where x < y
    
    # TRUE = 1, FALSE = 0
    # specificy coordinates for texts on plot
    mutate(!!x := lim/2*(2*(right - 0.5)+(UP - 0.5)+((conc-0.001)*0.5)-(dir_UP*1.5)+(dir_DOWN*0.5)), 
           !!y := lim/2*(2*(top - 0.5)+(UP - 0.5)))
  
  return(count)
}

colour_quadrant <-  function(dat, x, y, colx, coly, colNonCon){
  col <- dat %>%
    # logical columns to define where the point is locates
    mutate(right = .[[x]] > 0, # on the right part of plot?
           top = .[[y]] > 0, # on the top part of plot?
           # for the concordant quadrants (I & III)... 
           DOWN = !top & abs(.[[x]]) > abs(.[[y]]) | # quadrant III where x > y
             top & abs(.[[x]]) > abs(.[[y]])) %>% # quadrant I where x > y
    # add the colour
    mutate(quadrant = ifelse(right & top | (!right & !top), 
                             ifelse(DOWN, colx, coly), 
                             colNonCon))
  return(col)
}


splicing.MF.diff.dot.plot <- function(RedData, NRData, plotCol){
  tmp <- as_tibble(data.frame(FlyBaseID = RedData$FlyBaseID,
                    Red = RedData[[plotCol]],
                    NR = NRData[[plotCol]]))
  lim = max(abs(min(tmp$Red, na.rm = T)), 
            abs(min(tmp$NR, na.rm = T)), 
            max(tmp$Red, na.rm = T), 
            max(tmp$NR, na.rm = T), na.rm = T)
  lim=lim+(lim/2)
    
  # count the percentages
  quad_n <- quad_count(dat=tmp, x="Red", y="NR", lim=lim)
  # manage the colour of points
  quad_col <- colour_quadrant(tmp, "Red", "NR", "red3", "grey20", "grey60")

  splice.plot <- ggplot(tmp, aes(x = Red, y = NR)) +
    geom_point(size = 2, shape = 16, alpha = 0.5, color = quad_col$quadrant) +

    # add lines to separate quadrants
    geom_abline(intercept = 0, slope = 0,  size = 0.5, linetype="solid", color = "black") +
    geom_hline(yintercept = 0,  size = 0.5, linetype="solid", color = "black") +
    geom_vline(xintercept = 0,  size = 0.5, linetype="solid", color = "black") +
    geom_abline(intercept = 0, slope = 1,  size = 0.5, linetype="dashed", color = "black") +
    geom_abline(intercept = 0, slope = -1,  size = 0.5, linetype="dashed", color = "black") +

    # add percentages
    geom_text(aes(label = paste(round(perc*100,digits=0),"%",sep="")), data = quad_n, size = 10) +
    coord_cartesian(xlim=c(-lim, lim), ylim = c(-lim,lim)) +
    labs(x = "Red", y = "Non-Red") +
    guides(color = guide_legend(override.aes = list(shape = c(NA, NA), # c(16, 16)
                                                    size = c(4, 4),
                                                    alpha = 1))) +

    # some theme settings...
    theme_classic() +
    theme(legend.title = element_blank(),
          legend.position = c("None"),
          legend.box.background = element_rect(),
          legend.text = element_text(size = 20, color = "black"),
          plot.tag = element_text(size = 20, color = "black"),
          axis.text.x = element_text(size=20, margin = margin(5,0,0,0), color = "black"),
          axis.text.y = element_text(size=20, margin = margin(0,5,0,0), color = "black"),
          axis.title.x = element_text(size=30, margin = margin(10,0,0,0), color = "black"),
          axis.title.y = element_text(size=30, margin = margin(0,10,0,0), color = "black"),
          plot.title = element_text(size=40, margin = margin(0,0,0,0), color = "black"),
          plot.margin = margin(6,6,6,6)
    )
  return(splice.plot)
}


sampleTypes <- c("A.m", "A.f", "C.m")
M.v.F.test.table <- data.frame(sampleType = sampleTypes) %>% 
  mutate(Red.MvF.compare = paste0(sampleType,".Red.compare_SinghAgrw"), # .Red.compare_SinghAgrw .Red.compare_ase
         NR.MvF.compare = paste0(sampleType,".NR.compare_ase")) # .NR.compare_SinghAgrw .NR.compare_ase

for(i in 1:dim(M.v.F.test.table)[1]){
  RedData <- get(paste0(M.v.F.test.table[i,2]))
  # subset appropriately
  RedData <- RedData[ RedData$FlyBaseID %in% subset.sss,] # filter.low.exp.genes.q25
  # subset appropriately
  NRData <- get(paste0(M.v.F.test.table[i, 3]))
  NRData <- NRData[NRData$FlyBaseID %in% subset.sss,] # filter.low.exp.genes.q25
  
  M.v.F.test.table$N <- length(!is.na(RedData$M.dis.F))
  M.v.F.test.table$avg.Red[i] <- mean(RedData$M.dis.F, na.rm = T)
  
  M.v.F.test.table$avg.NR[i] <- mean(NRData$M.dis.F, na.rm = T)
  
  test <- merge(RedData, NRData, by = "FlyBaseID")
  M.v.F.test.table$pval[i] <- t.test(test$M.dis.F.x, 
                                     test$M.dis.F.y, paired = T)$p.value
  M.v.F.test.table$diff[i] <- t.test(test$M.dis.F.x, 
                                     test$M.dis.F.y, paired = T)$estimate
    
  rm(test)
  
  tmp.plot <- splicing.MF.diff.dot.plot(RedData = RedData, NRData = NRData, plotCol = "M.dis.F")
  
  assign(paste0(M.v.F.test.table[i,1],".compare.plot"), tmp.plot)
}

A.m.compare.plot 
A.f.compare.plot 
C.m.compare.plot

M.v.F.test.table

A.m.splice_ase.not.sss.cor.plot.dist <- A.m.compare.plot #+ coord_cartesian(xlim = c(-1,1))
A.f.splice_ase.not.sss.cor.plot.dist <- A.f.compare.plot #+ coord_cartesian(xlim = c(-1,1))
C.m.splice_ase.not.sss.cor.plot.dist <- C.m.compare.plot #+ coord_cartesian(xlim = c(-1,1))


A.f.splice_ase.q25.filtered.
A.m.splice_ase.q50.filtered.perc.sim
A.f.splice_ase.q75.filtered
A.f.splice_ase.q100.filtered

write.table(M.v.F.test.table, file = "Results/tmp.csv", sep = ",", quote = FALSE, row.names = F)


A.m.sig.sss.ase
A.m.nonsig.sss.ase

#######




pdf(file = "~/Desktop/UofT/SSAV_RNA/Plots/splicing/SSAVvSDIU.percentSim.pdf",   # The directory you want to save the file in
    width = 5, # 12, 24, 20 The width of the plot in inches
    height = 5) # 10, 20, 13 The height of the plot in inches

# A.m.compare.plot

# ggarrange(A.m.splice_sdiu.sss.perc.sim, NA, C.m.splice_sdiu.sss.perc.sim, NA, A.f.splice_sdiu.sss.perc.sim,
#           widths = c(1, 0.05, 1, 0.05, 1),
#           ncol = 5, labels = c("Am)", NA, "Cm)", NA, "Af)"),
#           font.label = list(size = 30))

# ggarrange(C.m.splice_ase.sss.cor.plot.dist ,
#           NA,
#           C.m.splice_ase.not.sss.cor.plot.dist,
#           heights = c(1, 0.05, 1),
#           nrow = 3, labels = c("A)", NA, "B)"),
#           font.label = list(size = 30))

# ggarrange(C.m.splice_sdiu.q25, #+ coord_cartesian(xlim = c(-1,1)),
#           NA,
#           C.m.splice_sdiu.q50,# + coord_cartesian(xlim = c(-1,1)),
#           NA, NA, NA,
#           C.m.splice_sdiu.q75,# + coord_cartesian(xlim = c(-1,1)),
#           NA,
#           C.m.splice_sdiu.q100,# + coord_cartesian(xlim = c(-1,1)),
#           heights = c(1, 0.05, 1),
#           widths = c(1, 0.05, 1),
#           nrow = 3, ncol = 3,
#           labels = c("A)", NA, "B)",
#                       NA, NA, NA,
#                       "C)", NA, "D)"),
#           font.label = list(size = 30))

plot(test[test$FlyBaseID %in% test.filter.25,]$percent.sim.x, 
     test[test$FlyBaseID %in% test.filter.25,]$percent.sim.y, xlab = "SSAV", ylab ="Osada")


dev.off()

###################################
#
#                                   Grieshop et al. 2024
#                                   Author: Michelle Liu
#                  DsRed experimental evolution - transcriptomics analysis
#            Masculinization/feminization of splicing profiles in Red vs. NonRed
# 
# 
###################################

# run JunctionSeq to generate count data for exons, 
# see: Amardeep Singh's JunctionSeq.Script.R,
# or Michelle's modified version for these populations "JunctionSeqRun.R"

rm(list=ls())
setwd("~/Desktop/UofT/SSAV_RNA/")


# required packages
#######
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(DESeq2)
# library(GenomicFeatures)
require(VennDiagram)
require(grDevices)
require(ggstatsplot)
library(broom)
library(ggblend)
#######

# Get chromosome locations  
source("Chromosome_df.R")


# set up external data from Mishra et al. 2022
########
jseq.Mishra = read.table("JunctionSeq/SDIU_ase/JSresults/SDIU_ASEallGenes.results.txt",
                      sep = "\t", header = TRUE)
colnames(jseq.Mishra)[2]="FlyBaseID" # change column name to reflect the rest of the dataset

# ------ First for the Mishra reference 
# remove genes that were not assayed in males and females
jseq.Mishra = jseq.Mishra[!is.na(jseq.Mishra$expr_F) & !is.na(jseq.Mishra$expr_M),]
dim(jseq.Mishra)

# remove novel counting bins
#   for the purpose of comparing Mishra with the SSAV populations, we want to make sure that the counting bins
#   used are the same.
jseq.Mishra = jseq.Mishra %>% filter(!str_detect(countbinID, "N"))
dim(jseq.Mishra) # check how many novel splice sites were removed

# assign genes with significant sex-specific splicing
jseq.Mishra <- jseq.Mishra %>% mutate(SSS = ifelse(geneWisePadj < 0.01, TRUE, FALSE))
########


# set list of genes to filter out based on FPKM calculated from Mishra et al's data
######
## calculate FPKM
# load metadata file for Mishra samples
decoder.Mishra <- read.table("JunctionSeq/SDIU_ase/QoRTs.decoder.file.for.JunctionSeq.txt", header=T, stringsAsFactors = F)
colnames(decoder.Mishra)[1]="unique.ID"
decoder.Mishra$rep = rep(seq(1:3),4)

# make list containing the paths to sample count files
countFiles.Mishra <- paste0("JunctionSeq/SDIU_ase/count.files/",decoder.Mishra$unique.ID,"/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")

# use this function to read and combine all sample counts and create a counts matrix
# data.files is a list containing the paths to sample count files
makeCountMatrix <- function(data.files, decoder.file){
  count.file <- NULL # initialize empty matrix to store count files
  # load and combine sample counts 
  for(i in data.files){
    tmp <- read.delim(i, header=F, sep = "\t") # read each count file
    tmp <- tmp %>% 
      separate(V1, into = c("FlyBaseID", "countBin"), sep = ":") %>% # split first column (ID:countBin)
      filter(str_detect(countBin, "A000")) # remove the global gene count
    if(!is.null(count.file)){
      count.file <- merge(count.file, tmp[,c("FlyBaseID", "V2")], by = "FlyBaseID", all = T)} # add to count matrix
    else
      count.file <- tmp[,c("FlyBaseID", "V2")] # store the counts in the new data frame
  }
  
  colnames(count.file) <- c("FlyBaseID", decoder.file$unique.ID) # set column names
  
  # remove bins that overlap multiple genes
  count.file <- count.file %>% 
    dplyr::filter(!str_detect(FlyBaseID, "\\+")) 
  
  # change geneIDs to rownames
  count.file <- count.file %>% 
    remove_rownames %>% 
    column_to_rownames(var="FlyBaseID")
  
  return(count.file)
}

Mishra.count.matrix <- makeCountMatrix(countFiles.Mishra, decoder.Mishra)

# get column metadata from decoder file
Mishra.colData <- decoder.Mishra %>% 
  remove_rownames %>% 
  column_to_rownames(var="unique.ID")

# make DESeq2 object to calculate fpkm
dds.Mishra <- DESeqDataSetFromMatrix(countData = Mishra.count.matrix, 
                                  colData = Mishra.colData, 
                                  design = ~ rep + sex)

# get gene lengths from GTF file 
# (done on server, scp to local so load this to env if you also did it that way)
# if not, just use the "exonic" GRanges object
# txdb <- makeTxDbFromGFF(gtfFile, format="gtf")
# exonic <- exonsBy(txdb, by="gene")
# save(exonic, file="/plas1/michelle.liu/Dmel_BDGP6.28/Drosophila_melanogaster.BDGP6.28.102.exonLengths.RData")
load("JunctionSeq/SDIU_ase/Drosophila_melanogaster.BDGP6.28.102.exonLengths.RData", verbose = T)

exonic <- GRangesList(exonic)
exonic <- exonic[names(exonic) %in% rownames(Mishra.count.matrix)] # keep only gene lengths that are in the count matrix
rowRanges(dds.Mishra) <- exonic # add gene length data to DESeq2 object
FPKM.Mishra <- fpkm(dds.Mishra) # get FPKM measures

# separate the counts to males and females
FPKM.Mishra.fem <- data.frame(FPKM.Mishra) %>% 
  select(., contains("F_")) %>%
  mutate(totalCounts = rowSums(select_if(., is.numeric))) # get total counts for females
FPKM.Mishra.fem[FPKM.Mishra.fem==0] <- NA # set genes not present as NAs

FPKM.Mishra.male <- data.frame(FPKM.Mishra) %>% 
  select(., contains("M_")) %>%
  mutate(totalCounts = rowSums(select_if(., is.numeric))) # get total counts for males
FPKM.Mishra.male[FPKM.Mishra.male==0] <- NA # set genes not present as NAs


###### 25% FILTER
# list of genes with FPKM > the 25% cut-off in females
filter.low.exp.genes.fem.q25 <- rownames(FPKM.Mishra.fem[!is.na(FPKM.Mishra.fem$totalCounts) &
                                                        FPKM.Mishra.fem$totalCounts  > quantile(FPKM.Mishra.fem$totalCounts, 0.25, na.rm=T),])
# list of genes with FPKM > the 25% cut-off in males
filter.low.exp.genes.male.q25 <- rownames(FPKM.Mishra.male[!is.na(FPKM.Mishra.male$totalCounts) &
                                                          FPKM.Mishra.male$totalCounts  > quantile(FPKM.Mishra.male$totalCounts, 0.25, na.rm=T),])
# combine list of genes that passed filtering
filter.low.exp.genes.q25 <- unique(c(filter.low.exp.genes.fem.q25, filter.low.exp.genes.male.q25))
# remove Y chr genes that somehow gets there(?) could be from spermatheca in females?
filter.low.exp.genes.q25 <- filter.low.exp.genes.q25[!(filter.low.exp.genes.q25 %in% Chrs[Chrs$Chr == "Y",]$FlyBaseID)]
length(filter.low.exp.genes.q25) # check how many genes are left



###### 10% FILTER
# list of genes with FPKM > the 10% cut-off in females
filter.low.exp.genes.fem.q10 <- rownames(FPKM.Mishra.fem[!is.na(FPKM.Mishra.fem$totalCounts) &
                                                        FPKM.Mishra.fem$totalCounts  > quantile(FPKM.Mishra.fem$totalCounts, 0.10, na.rm=T),])
# list of genes with FPKM > the 10% cut-off in males
filter.low.exp.genes.male.q10 <- rownames(FPKM.Mishra.male[!is.na(FPKM.Mishra.male$totalCounts) &
                                                          FPKM.Mishra.male$totalCounts  > quantile(FPKM.Mishra.male$totalCounts, 0.10, na.rm=T),])

# combine list of genes that passed filtering
filter.low.exp.genes.q10 <- unique(c(filter.low.exp.genes.fem.q10, filter.low.exp.genes.male.q10))
# remove Y chr genes that somehow gets there(?) could be from spermatheca in females?
filter.low.exp.genes.q10 <- filter.low.exp.genes.q10[!(filter.low.exp.genes.q10 %in% Chrs[Chrs$Chr == "Y",]$FlyBaseID)]
length(filter.low.exp.genes.q10) # check how many genes are left
######
# results:
length(filter.low.exp.genes.q10) # list of genes that passed the 10% filter
length(filter.low.exp.genes.q25) # list of genes that passed the 25% filter



####### for the genes with significant alt. splicing in Red and Non-Red, ##########
################# are they more masculinized or feminized? #######################

# --------------- analyse differences in splicing profiles
# using count.files (read counts per exon, generated by QoRTs. See JunctionSeq_Linux.sh)

####### set up decoder files for each sample group ########
##### load decoder files

# this is a decoder file with:
# $1:unique.sample.ID, $2:cdt1, $3:cdt2, $n:cdtn 
# (See JunctionSeq_Linux.sh on how to make this)
decoder <- read.table("JunctionSeq/QoRTs.decoder.file.for.JunctionSeq.txt", header = TRUE, stringsAsFactors = FALSE)


# decoder files for external datasets
# Mishra et al. 2022
decoder.Mishra <- read.table("JunctionSeq/SDIU_ase/QoRTs.decoder.file.for.JunctionSeq.txt", header=T, stringsAsFactors = F)
colnames(decoder.Mishra)[1]="unique.ID"

# Singh & Agrawal 2023
decoder.Osada <- read.table("JunctionSeq/SDIU_singh/QoRTs.decoder.SinghAgrawal.body.txt", header = FALSE, sep = "\t")
colnames(decoder.Osada) <- c("unique.ID", "sex")

# Mishra et al. 2023
decoder.MCsim <- read.table("JunctionSeq/SDIU_MC/MCsim/decoder.MCsim.txt", header = F, sep="\t")
colnames(decoder.MCsim) <- c("unique.ID", "sex")
decoder.MCcom <- read.table("JunctionSeq/SDIU_MC/MCcom/decoder.MCcom.txt", header=F, sep="\t")
colnames(decoder.MCcom) <- c("unique.ID", "sex")
decoder.MCabs <- read.table("JunctionSeq/SDIU_MC/MCabs/decoder.MCabs.txt", header=F, sep="\t")
colnames(decoder.MCabs) <- c("unique.ID", "sex")


######## load JunctionSeq sizeFactors (for per sample normalization factors) 
# this has to be created through JunctionSeq. 
# You can technically calculate the sizeFactors yourself, but JunctionSeq can give you the by countingBin and by gene factors so it is easier.
# look at DESeq2's geometric normalization calculation if you want to do it yourself
SSAV.factors <- read.delim("JunctionSeq/RAL.size.Factors.GEO.txt", header=T, sep="\t")
# separate by sample type
A.f.factors <- SSAV.factors  %>% filter(str_detect(sample.ID, "\\_F_"))
A.m.factors <- SSAV.factors  %>% filter(str_detect(sample.ID, "A[[:digit:]]\\_M_"))
C.m.factors <- SSAV.factors  %>% filter(str_detect(sample.ID, "C[[:digit:]]\\_M_"))

Mishra.factors <- read.delim("JunctionSeq/SDIU_ase/RAL.size.Factors.GEO.txt", header = T, sep="\t")
# separate by sex
F.Mishra.factors <- Mishra.factors  %>% filter(str_detect(sample.ID, "\\F_"))
M.Mishra.factors <- Mishra.factors  %>% filter(str_detect(sample.ID, "\\M_"))

Osada.factors <- read.delim("JunctionSeq/SDIU_singh/RAL.size.Factors.GEO.txt", header=T, sep="\t")
# separate by sex
F.Osada.factors <- Osada.factors %>% filter(str_detect(sample.ID, "\\.female"))
M.Osada.factors <- Osada.factors %>% filter(str_detect(sample.ID, "\\.male"))


MC.factors <- read.delim("JunctionSeq/SDIU_MC/JS.GEO.size.factors.txt", header=T, sep="\t")
# separate by sex & mating system treatment
F.MCabs.factors <- MC.factors %>% filter(str_detect(sample.ID, "M[[:digit:]]\\_Female_body"))
M.MCabs.factors <- MC.factors %>% filter(str_detect(sample.ID, "M[[:digit:]]\\_Male_body"))
F.MCsim.factors <- MC.factors %>% filter(str_detect(sample.ID, "P[[:digit:]]\\_Female_body"))
M.MCsim.factors <- MC.factors %>% filter(str_detect(sample.ID, "P[[:digit:]]\\_Male_body"))
F.MCcom.factors <- MC.factors %>% filter(str_detect(sample.ID, "C[[:digit:]]\\_Female_body"))
M.MCcom.factors <- MC.factors %>% filter(str_detect(sample.ID, "C[[:digit:]]\\_Male_body"))



######## separate decoder files by sample groups & merge with sizeFactors
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


# Mishra et al. 2022 Mishra data for reference
M.decoder.Mishra <- decoder.Mishra[decoder.Mishra$sex=="M",]
M.decoder.Mishra <- merge(M.decoder.Mishra, M.Mishra.factors, by = 1)

F.decoder.Mishra <- decoder.Mishra[decoder.Mishra$sex=="F",]
F.decoder.Mishra <- merge(F.decoder.Mishra, F.Mishra.factors, by =1)

countFiles.m.Mishra <- paste0("JunctionSeq/SDIU_ase/count.files/",M.decoder.Mishra$unique.ID,"/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")
countFiles.f.Mishra <- paste0("JunctionSeq/SDIU_ase/count.files/",F.decoder.Mishra$unique.ID,"/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")




# Singh & Agrawal 2023 (Osada) data for another reference
M.decoder.Osada <- decoder.Osada[decoder.Osada$sex=="male",]
M.decoder.Osada <- merge(M.decoder.Osada, M.Osada.factors, by = 1)

F.decoder.Osada <- decoder.Osada[decoder.Osada$sex=="female",]
F.decoder.Osada <- merge(F.decoder.Osada, F.Osada.factors, by = 1)

countFiles.m.Osada <- paste0("JunctionSeq/SDIU_singh/body.only.replicate.1/",M.decoder.Osada$unique.ID,"/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")
countFiles.f.Osada <- paste0("JunctionSeq/SDIU_singh/body.only.replicate.1/",F.decoder.Osada$unique.ID,"/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")



# Mishra et al. 2023 MC data
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



# ------ analysis functions:
# data.files is a list containing the paths to sample count files
# decoder.file has the metadata associated with each count file
# this function is used inside GeomNormCounts
ReadCountFiles <- function(data.files, decoder.file){
  count.file <- NULL
  
  for(i in data.files){
    tmp <- read.delim(i, header=F, sep = "\t")
    count.file <- cbind(count.file, tmp[,2]) #add new column for the other sample
  }
  
  countBins <- tmp %>% separate(V1, into = c("FlyBaseID", "countBin"), sep = ":")
  count.file <- cbind(countBins[,1:2], count.file)
  colnames(count.file) <- c("FlyBaseID", "countBin", decoder.file$unique.ID)
  print(dim(count.file)) # check number of splice sites
  
  # exclude novel splice sites, exclude global gene sites, and exclude reads that overlap two countbins/gene
  count.file <- count.file %>% 
    filter(!str_detect(countBin, "\\+") & !str_detect(countBin, "A") & !str_detect(countBin, "N"))
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
  count.file <- count.file %>% 
    dplyr::mutate(total=rowSums(select_if(., is.numeric), na.rm = T)) # sum per row
  print(dim(count.file)) # check number of splice sites
  
  # calculate the global gene expression across all samples
  glob.exp <- count.file %>% group_by(FlyBaseID) %>% # sum per gene
    dplyr::summarise(glob.exp = sum(total, na.rm = T))
  
  # calculate proportion of reads by gene for each counting bin
  count.file <- merge(count.file, glob.exp, by = "FlyBaseID", all = T) %>%
    dplyr::mutate(frac.exp.per.gene = ifelse(!is.na(glob.exp), total/glob.exp, NA))
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
    dplyr::group_by(FlyBaseID) %>%
    dplyr::summarise(percent.sim = PercentSimilarity(frac.exp.per.gene.x, frac.exp.per.gene.y),
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

# this function calculates phi,
# which essentially compares the relative "distance" between point A and point B in 
# a multidimensional gene space. (i.e., is the sample closer to point A or to point B ?)
CalculatePhi <- function(dst1, dst2){ 
  tmp <- merge(dst1, dst2, by = "FlyBaseID") %>%
    dplyr::mutate(
      phi = (percent.dissim.y - percent.dissim.x)/(percent.dissim.y + percent.dissim.x))
  return(tmp)
}



# Load count files:
# -------------------------
# average females expression in Experimental SSAV populations (combining Red and NR samples)
########
A.f.Red.norm.exp <- GeomNormCounts(countFiles.A.f.Red, A.f.Red.decoder)
A.f.NR.norm.exp <- GeomNormCounts(countFiles.A.f.NR, A.f.NR.decoder)

# average female in SSAV pop
fem.exp <- merge(A.f.Red.norm.exp, A.f.NR.norm.exp, by = c("FlyBaseID", "countBin")) %>%
  dplyr::group_by(FlyBaseID, countBin) %>%
  dplyr::summarise(total = (total.x + total.y)/2,
                   glob.exp = (glob.exp.x + glob.exp.y)/2,
                   frac.exp.per.gene = (frac.exp.per.gene.x + frac.exp.per.gene.y)/2)

########

# average male expression in Experimental SSAV populations (combining Red and NR samples)
########
# normalize exp. per gene for each counting bin 
# in Red males
A.m.Red.norm.exp <- GeomNormCounts(countFiles.A.m.Red, A.m.Red.decoder)
# in Non-Red males
A.m.NR.norm.exp <- GeomNormCounts(countFiles.A.m.NR, A.m.NR.decoder)

# average male in SSAV pop
male.exp <- merge(A.m.Red.norm.exp, A.m.NR.norm.exp, by = c("FlyBaseID", "countBin")) %>%
  dplyr::group_by(FlyBaseID, countBin) %>%
  dplyr::summarise(total = (total.x + total.y)/2,
                   glob.exp = (glob.exp.x + glob.exp.y)/2,
                   frac.exp.per.gene = (frac.exp.per.gene.x + frac.exp.per.gene.y)/2)

#######

# average male expression in Control populations (combining Red and NR samples)
########
# expression of NR C males
C.m.Red.norm.exp <- GeomNormCounts(countFiles.C.m.Red, C.m.Red.decoder)
C.m.NR.norm.exp <- GeomNormCounts(countFiles.C.m.NR, C.m.NR.decoder)

male.ctr.exp <- merge(C.m.Red.norm.exp, C.m.NR.norm.exp, by = c("FlyBaseID", "countBin")) %>%
  dplyr::group_by(FlyBaseID, countBin) %>%
  dplyr::summarise(total = (total.x + total.y)/2,
                   glob.exp = (glob.exp.x + glob.exp.y)/2,
                   frac.exp.per.gene = (frac.exp.per.gene.x + frac.exp.per.gene.y)/2)
########


# male and female profiles in external populations
######
# using Mishra data
male.exp_ASE <- GeomNormCounts(countFiles.m.Mishra, M.decoder.Mishra)
fem.exp_ASE <- GeomNormCounts(countFiles.f.Mishra, F.decoder.Mishra)

# Singh & Agrawal (Osada) data
male.exp_Osada <- GeomNormCounts(countFiles.m.Osada, M.decoder.Osada)
fem.exp_Osada <- GeomNormCounts(countFiles.f.Osada, F.decoder.Osada)

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


# Get a subset of SSS genes that are more consistently dimorphic
#######
Mishra.sig.SSS <- read.delim("JunctionSeq/SDIU_ase/JSresults/ASE.sig.SSS_genes.txt", header = F, skip = 1)
Mishra.sig.SSS <- Mishra.sig.SSS$V1

Mishra.sig.SSS_filt25 <- Mishra.sig.SSS[
  Mishra.sig.SSS %in% filter.low.exp.genes.q25]
length(Mishra.sig.SSS_filt25)


# checking for irregularites in MC population profiles
MC_data <- data.frame(MC = c("abs", "sim", "com"),
                      count_F = c("fem.exp_MCabs", "fem.exp_MCsim", "fem.exp_MCcom"),
                      count_M = c("male.exp_MCabs", "male.exp_MCsim", "male.exp_MCcom"))

for(i in 1:dim(MC_data)[1]){
  count_F <- get(paste0(MC_data[i,2]))
  count_M <- get(paste0(MC_data[i, 3]))
  
  Male.Mishra.dist.MCf <- compareSplicingProfiles(male.exp_ASE, count_F)
  Male.Mishra.dist.MCm <- compareSplicingProfiles(male.exp_ASE, count_M)
  Fem.Mishra.dist.MCf <- compareSplicingProfiles(fem.exp_ASE, count_F)
  Fem.Mishra.dist.MCm <- compareSplicingProfiles(fem.exp_ASE, count_M)
  
  # dot plots to check how many points fall into the "irregular" category, but for the MC populations
  # "irregular" meaning that when looking at males (females), 
  # the points fall closer to the reference females (males)
  test.M <- merge(Male.Mishra.dist.MCm, Male.Mishra.dist.MCf, by = "FlyBaseID")
  test.F <- merge(Fem.Mishra.dist.MCm, Fem.Mishra.dist.MCf, by = "FlyBaseID")
  # plot(test.F[test.F$FlyBaseID %in% Mishra.sig.SSS_filt25 & 
  #                  !is.na(test.F$percent.dissim.x) & !is.na(test.F$percent.dissim.y), ], 
  #         x = "percent.dissim.x", y = "percent.dissim.y", 
  #         colx = "black", coly = "black",colNonCon = "black",xlab = "to males",
  #         ylab = "to females", lim = 1, title = "") + coord_cartesian(xlim = c(0,1), ylim = c(0,1))
  
  # get a list of all irregular genes from each MC population treatment
  # list of genes in females where expression profile is more dissimilar to female ref than male ref
  fem <- test.F[test.F$FlyBaseID %in% Mishra.sig.SSS_filt25 & 
                  !is.na(test.F$percent.dissim.x) & !is.na(test.F$percent.dissim.y) &
                  test.F$percent.dissim.x <= test.F$percent.dissim.y, ]$FlyBaseID
  male <- test.M[test.M$FlyBaseID %in% Mishra.sig.SSS_filt25 & 
                   !is.na(test.M$percent.dissim.x) & !is.na(test.M$percent.dissim.y) &
                   test.M$percent.dissim.x >= test.M$percent.dissim.y, ]$FlyBaseID
  
  assign(paste0("fem.",MC_data[i, 1]), fem)
  assign(paste0("male.",MC_data[i, 1]), male)
}

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
# 1) SSS in the Mishra population
# 2) fall into the right "sex-profile" (i.e., not "irregular") when seen in the MC populations
# Extra filtering to removing genes near the DsRed marker
DsRed_genes <- read.delim(file="~/Desktop/UofT/SSAV_RNA/Data/dmel_2R_DsRed_ids.tsv", header=FALSE)
subset.sss <- (Mishra.sig.SSS_filt25[!Mishra.sig.SSS_filt25 %in% all.exclude &
                                    !Mishra.sig.SSS_filt25 %in% DsRed_genes$V1])
length(all.exclude)
length(subset.sss)
write_delim(data.frame(subset.sss), file = "JunctionSeq/dimorphic.subset.list.txt", delim = ",", col_names = F)

subset.sss <- read.delim("JunctionSeq/dimorphic.subset.list.txt", header = F)
subset.sss <- subset.sss$V1

#######



# -------- male comparison -------------
########
# expression of Red Experimental males 
# to reference Mishra males
A.m.Red.v.Males_ase <- compareSplicingProfiles(A.m.Red.norm.exp, male.exp_ASE)
# to reference Mishra females
A.m.Red.v.Fem_ase <- compareSplicingProfiles(A.m.Red.norm.exp, fem.exp_ASE)
# calculate Phi for Red males
A.m.Red.Phi <- CalculatePhi(A.m.Red.v.Males_ase, A.m.Red.v.Fem_ase) 

# ---
# expression of NR Experimental males
A.m.NR.v.Males_ase <- compareSplicingProfiles(A.m.NR.norm.exp, male.exp_ASE) # with Mishra data
A.m.NR.v.Fem_ase <- compareSplicingProfiles(A.m.NR.norm.exp, fem.exp_ASE)
A.m.NR.Phi <- CalculatePhi(A.m.NR.v.Males_ase, A.m.NR.v.Fem_ase)
#######


# -------- female comparison -------------
########
# expression of Red Experimental females
A.f.Red.v.Males_ase <- compareSplicingProfiles(A.f.Red.norm.exp, male.exp_ASE)
A.f.Red.v.Fem_ase <- compareSplicingProfiles(A.f.Red.norm.exp, fem.exp_ASE)
A.f.Red.Phi <- CalculatePhi(A.f.Red.v.Males_ase, A.f.Red.v.Fem_ase)

# expression of NR Experimental females
A.f.NR.v.Males_ase <- compareSplicingProfiles(A.f.NR.norm.exp, male.exp_ASE)
A.f.NR.v.Fem_ase <- compareSplicingProfiles(A.f.NR.norm.exp, fem.exp_ASE)
A.f.NR.Phi <- CalculatePhi(A.f.NR.v.Males_ase, A.f.NR.v.Fem_ase)
#######

# -------- control male comparison -------------
#######
# Control Red males
C.m.Red.v.Males_ase <- compareSplicingProfiles(C.m.Red.norm.exp, male.exp_ASE)
C.m.Red.v.Fem_ase <- compareSplicingProfiles(C.m.Red.norm.exp, fem.exp_ASE)
C.m.Red.Phi <- CalculatePhi(C.m.Red.v.Males_ase, C.m.Red.v.Fem_ase)

# Control NR males
C.m.NR.v.Males_ase <- compareSplicingProfiles(C.m.NR.norm.exp, male.exp_ASE)
C.m.NR.v.Fem_ase <- compareSplicingProfiles(C.m.NR.norm.exp, fem.exp_ASE)
C.m.NR.Phi <- CalculatePhi(C.m.NR.v.Males_ase, C.m.NR.v.Fem_ase)
#######


# Figure plotting functions
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

# histogram plotting only the differnce in phi between two sample types
splicing.MFdiff.plot <- function(RedData, NRData, plotCol, color){
  tmp <- data.frame(FlyBaseID = RedData$FlyBaseID,
                    diff = RedData[[plotCol]] - NRData[[plotCol]])
  
  splice.plot <- ggplot() + 
    geom_histogram(data= tmp, 
                   aes(x= diff), color = color, fill = color, alpha = 0.5, binwidth = 0.08) +
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
          plot.margin = margin(6,6,6,6),
          axis.line.x = element_blank(),
          axis.line.y = element_blank())
  
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
  dat <- dat[dat[[x]] != dat[[y]],]
  
  count <- dat %>%
    # Count how many with each combination of X and Y being positive
    dplyr::count(right = .[[x]] > 0, # on the right side of plot?
                 top = .[[y]] > 0, # on the top side of plot?
                 UP = (!top & !right & .[[x]] < .[[y]]) | # quadrant III x < y
                   (top & right & .[[x]] > .[[y]])) %>%  # quadrant I x > y
    dplyr::mutate(perc = n/sum(n, na.rm = T)) %>% # calculate percentage of points relative to total number of points
    
    # this is another strange one for setting up the coordinates
    dplyr::mutate(conc = right & top | (!right & !top), # quadrant I and quadrant III (the concordant changes)
                  dir_UP = conc & UP, # concordant changes where x > y
                  dir_DOWN = conc & !UP) %>% # concordant changes where x < y
    
    # TRUE = 1, FALSE = 0
    # specificy coordinates for texts on plot
    # specificy coordinates for texts on plot
    dplyr::mutate(!!x := ifelse(!top, lim/2*(2*(right - 0.5)+(UP - 0.5)+((conc-0.001)*0.5)-(dir_UP*1.5)+(dir_DOWN*0.5)), 
                                ifelse(UP, lim/2*(2*(right-0.5)+(UP+0.5)+((conc-0.001)*0.5)-(dir_UP*1.5)+(dir_DOWN*0.5)), 
                                       lim/2*(2*(.05)+(-0.5)+((conc-0.001)*0.5)-(dir_UP*1.5)+(dir_DOWN*0.5)))),
                  
                  !!y := ifelse(!top, lim/2*(2*(top - 0.5)+(UP-0.5)), 
                                ifelse(UP, lim/2*(2*(top - 0.5)+(-0.5)), 
                                       lim/2*(2*(top - 0.5)+(0.5))) )) %>%
    # 
    dplyr::mutate(!!x := ifelse((!right & top & !UP), .[[x]]*2,
                                ifelse(right & !top & !UP, .[[x]]*2 , .[[x]])),
                  !!y := ifelse((!right & top & !UP), .[[y]]/1.5,
                                ifelse(right & !top & !UP, .[[y]]/1.5, .[[y]])))
  # 
  #   # adjustments for males
  #   dplyr::mutate(!!x := ifelse((!right & !top & UP), .[[x]]/1.5,
  #                               ifelse((!right & !top & !UP), .[[x]]/1.2, .[[x]] )),
  #                 !!y := ifelse((!right & !top & !UP), .[[y]]/1.3,
  #                               ifelse((!right & !top & UP), .[[y]]/1.2, .[[y]] )))
  
  
  print(count)
  
  return(count)
}

colour_quadrant <-  function(dat, x, y, colx, coly, colNonCon){
  col <- dat %>%
    # logical columns to define where the point is locates
    dplyr::mutate(right = .[[x]] > 0, # on the right part of plot?
                  top = .[[y]] > 0, # on the top part of plot?
                  # for the concordant quadrants (I & III)... 
                  DOWN = !top & abs(.[[x]]) > abs(.[[y]]) | # quadrant III where x > y
                    top & abs(.[[x]]) > abs(.[[y]])) %>% # quadrant I where x > y
    # add the colour
    dplyr:: mutate(quadrant = ifelse(right & top | (!right & !top), 
                                     ifelse(DOWN, colx, coly), 
                                     colNonCon))
  return(col)
}


# dot plot (lines can be commented in to add number/% of points for each octant)
splicing.MF.diff.dot.plot <- function(RedData, NRData, plotCol, color){
  tmp <- as_tibble(data.frame(FlyBaseID = RedData$FlyBaseID,
                              Red = RedData[[plotCol]],
                              NR = NRData[[plotCol]]))
  
  tmp <- na.omit(tmp)
  
  lim = max(abs(min(tmp$Red, na.rm = T)), 
            abs(min(tmp$NR, na.rm = T)), 
            max(tmp$Red, na.rm = T), 
            max(tmp$NR, na.rm = T), na.rm = T)
  lim = lim + 0.3
  
  # count the percentages
  quad_n <- quad_count(dat=tmp, x="Red", y="NR", lim=lim)
  # manage the colour of points
  # quad_col <- colour_quadrant(tmp, "Red", "NR", "grey2", "grey2", "grey20")
  
  splice.plot <- ggplot(tmp, aes(x = Red, y = NR)) +
    geom_point(size = 2, shape = 16, alpha = 0.5, color = color) + # quad_col$quadrant,
    # #CC3399, #0072B2, #666666
    
    # add lines to separate quadrants
    geom_abline(intercept = 0, slope = 0,  size = 0.5, linetype="solid", color = "black") +
    geom_hline(yintercept = 0,  size = 0.5, linetype="solid", color = "black") +
    geom_vline(xintercept = 0,  size = 0.5, linetype="solid", color = "black") +
    geom_abline(intercept = 0, slope = 1,  size = 0.5, linetype="dashed", color = "black") +
    # geom_abline(intercept = 0, slope = -1,  size = 0.5, linetype="dashed", color = "black") +
    
    # add percentages
    # geom_text(aes(label = paste(round(perc*100,digits=0),"%",sep="")), data = quad_n, size = 8.5) +
    coord_cartesian(xlim=c(-lim, lim), ylim = c(-lim,lim)) +
    labs(x = "Red", y = "Non-Red") +
    scale_x_continuous(breaks = c(-1, 0, 1)) +
    scale_y_continuous(breaks = c(-1, 0, 1)) +
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
          plot.margin = margin(6,6,6,6),
          axis.line.x = element_blank(),
          axis.line.y = element_blank()
    )
  return(splice.plot)
}




# figures and t-tests comparing (M-F)/(M+F) metric (PHI)
# generates Suppl. Table S6
######
sampleTypes <- c("A.m", "A.f", "C.m")
Phi.test.table <- data.frame(sampleType = sampleTypes) %>% 
  mutate(Red.MvF.compare = paste0(sampleType,".Red.Phi"), 
         NR.MvF.compare = paste0(sampleType,".NR.Phi")) 

for(i in 1:dim(Phi.test.table)[1]){
  RedData <- get(paste0(Phi.test.table[i,2]))
  # subset appropriately
  RedData <- RedData[ 
    RedData$FlyBaseID %in% subset.sss,]
  
  # subset appropriately
  NRData <- get(paste0(Phi.test.table[i, 3]))
  NRData <- NRData[
    NRData$FlyBaseID %in% subset.sss,]
  
  Phi.test.table$N <- length(!is.na(RedData$phi))
  Phi.test.table$avg.Red[i] <- mean(RedData$phi, na.rm = T)
  
  Phi.test.table$avg.NR[i] <- mean(NRData$phi, na.rm = T)
  
  test <- merge(RedData, NRData, by = "FlyBaseID")
  Phi.test.table$pval[i] <- t.test(test$phi.x, 
                                     test$phi.y, paired = T)$p.value
  Phi.test.table$diff[i] <- t.test(test$phi.x, 
                                     test$phi.y, paired = T)$estimate
  
  rm(test)
}

Phi.test.table

# write.table(Phi.test.table, file = "Results/tmp.csv", sep = ",", quote = FALSE, row.names = F)


#######


# Plotting Fig.5
########
Fig5A.A.m <- splicing.MF.diff.dot.plot(A.m.Red.Phi[A.m.Red.Phi$FlyBaseID %in% subset.sss,],
                                       A.m.NR.Phi[A.m.NR.Phi$FlyBaseID %in% subset.sss,], 
                                                  plotCol = "phi", color = "#0072B2") +
  coord_cartesian(xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1))

Fig5A.C.m <- splicing.MF.diff.dot.plot(C.m.Red.Phi[C.m.Red.Phi$FlyBaseID %in% subset.sss,],
                                       C.m.NR.Phi[C.m.NR.Phi$FlyBaseID %in% subset.sss,],
                                                  plotCol = "phi", color = "#666666") +
  coord_cartesian(xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1))

Fig5A.A.f <- splicing.MF.diff.dot.plot(A.f.Red.Phi[A.f.Red.Phi$FlyBaseID %in% subset.sss,],
                                       A.f.NR.Phi[A.f.NR.Phi$FlyBaseID %in% subset.sss,], 
                                                  plotCol = "phi", color = "#D55E00") +
  coord_cartesian(xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1))



Fig_5A <- ggarrange(NA,
                    Fig5A.A.f + ggtitle(expression(bold("Exp. Females"))) +
                         theme(axis.title.x = element_blank(), 
                               axis.title.y = element_blank(),
                               plot.title = element_text(hjust = 0.5, size = 30, vjust = 1.5, color = "#D55E00"),
                               panel.border = element_rect(colour = "#D55E00", fill=NA, size=3)),
                       NA,  
                    Fig5A.A.m + ggtitle(expression(bold("Exp. Males"))) +
                         theme(axis.title.x = element_blank(), 
                               axis.title.y = element_blank(),
                               plot.title = element_text(hjust = 0.5, size= 30, vjust = 1.5, color = "#0072B2"),
                               panel.border = element_rect(colour = "#0072B2", fill=NA, size=3)),
                       NA, 
                    Fig5A.C.m + ggtitle(expression(bold("Ctrl. Males"))) +
                         theme(axis.title.x = element_blank(), 
                               axis.title.y = element_blank(),
                               plot.title = element_text(hjust = 0.5, size = 30, vjust = 1.5, color = "#666666"),
                               panel.border = element_rect(colour = "#666666", fill=NA, size=3)),
                       widths = c(0.025, 1, 0.05, 1, 0.05, 1),
                       ncol = 6)


Fig_5A <- annotate_figure(Fig_5A, left = text_grob(expression(bold(italic(phi)["NonRed"])), 
                                                      rot = 90, size = 40),
                          bottom = text_grob(expression(bold(italic(phi)["Red"])), 
                                             size = 40))


Fig5B.A.m <- splicing.MFdiff.plot(A.m.Red.Phi[A.m.Red.Phi$FlyBaseID %in% subset.sss,],
                                  A.m.NR.Phi[A.m.NR.Phi$FlyBaseID %in% subset.sss,], 
                                          plotCol = "phi", color = "#0072B2") +
  annotate("label", label = expression(atop(bar(x)*" = 0.049", italic("P")*" < 10"^-5*"***")), 
           x =  1.05, y = 315, size = 8.5, label.padding=unit(1, "lines"))

Fig5B.A.f <- splicing.MFdiff.plot(A.f.Red.Phi[A.f.Red.Phi$FlyBaseID %in% subset.sss,],
                                  A.f.NR.Phi[A.f.NR.Phi$FlyBaseID %in% subset.sss,], 
                                          plotCol = "phi", color = "#D55E00") +
  annotate("label", label = expression(atop(bar(x)*" = -0.008", italic("P")*" = 0.014"*"*")),
           x =  -1.05, y = 625, size = 8.5,  label.padding=unit(1, "lines"))

Fig5B.C.m <- splicing.MFdiff.plot(C.m.Red.Phi[C.m.Red.Phi$FlyBaseID %in% subset.sss,],
                                  C.m.NR.Phi[C.m.NR.Phi$FlyBaseID %in% subset.sss,], 
                                          plotCol = "phi", color = "#666666")+
  annotate("label", label = expression(atop(bar(x)*" = -0.015", italic("P")*" < 10"^-5*"***")), 
           x =  -1.05, y = 405, size = 8.5,  label.padding=unit(1, "lines"))


Fig_5B <- ggarrange(NA,
                       Fig5B.A.f + coord_cartesian(xlim=c(-2,2)) +
                         theme(axis.title.x = element_blank(), 
                               axis.title.y = element_blank(),
                               plot.title = element_text(hjust = 0.5, size = 30, vjust = 1.5, color = "#D55E00"),
                               panel.border = element_rect(colour = "#D55E00", fill=NA, size=3)),
                       NA,  
                       Fig5B.A.m  + coord_cartesian(xlim=c(-2,2)) +
                         theme(axis.title.x = element_blank(), 
                               axis.title.y = element_blank(),
                               plot.title = element_text(hjust = 0.5, size= 30, vjust = 1.5, color = "#0072B2"),
                               panel.border = element_rect(colour = "#0072B2", fill=NA, size=3)),
                       NA, 
                       Fig5B.C.m  + coord_cartesian(xlim=c(-2,2)) +
                         theme(axis.title.x = element_blank(), 
                               axis.title.y = element_blank(),
                               plot.title = element_text(hjust = 0.5, size = 30, vjust = 1.5, color = "#666666"),
                               panel.border = element_rect(colour = "#666666", fill=NA, size=3)),
                       widths = c(0.025, 1, 0.05, 1, 0.05, 1),
                       ncol = 6)

Fig_5B <- annotate_figure(Fig_5B, left = text_grob("Count", 
                                                      rot = 90, size = 40),
                          bottom = text_grob(expression(bold(italic(Delta))*bold(italic(phi))), 
                                             size = 40))


# comment in to save plot
# pdf(file = "~/Desktop/UofT/SSAV_RNA/Plots/final_2/Fig5_main.pdf",   # The directory you want to save the file in
#     width = 21, # 12, 24, 20 The width of the plot in inches
#     height = 14) # 10, 20, 13 The height of the plot in inches
# 
# ggarrange(Fig_5A, 
#           NA, 
#           Fig_5B, 
#           nrow = 3, heights = c(1, 0.07, 1),
#           labels = c("A)", NA, "B)"),
#           font.label = list(size = 40))
# 
# dev.off()

########


# Plotting Fig.6B
########
# compare SSAV vs Controls
A.m.C.m.Red <- compareSplicingProfiles(A.m.Red.norm.exp, C.m.Red.norm.exp)
A.m.C.m.NR <- compareSplicingProfiles(A.m.NR.norm.exp, C.m.NR.norm.exp)
test <- merge(A.m.C.m.Red, A.m.C.m.NR, by = "FlyBaseID")
test <- na.omit(test)
# test <- test[test$FlyBaseID %in% subset.sss,]

# use the file "Figure6.R" to set up plotting functions
source("Fig6_effectFuns.R")

Figure_6B <- plot_corr(test, x="percent.dissim.x", y="percent.dissim.y", 
                       colx = "red3", coly = "grey9", colNonCon = "grey",
                       xlab = "Red", ylab = "NonRed", title = "", lim = 1.2) +
  labs(x = expression(italic("d"["Red,C"])),
       y = expression(italic("d"["NonRed,C"]))) +
  theme(axis.title.x = element_text(vjust = 1.5, margin = margin(25,0,20,0)), 
        axis.title.y = element_text(vjust = 1, margin = margin(0,20,0,20))) +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1))
Figure_6B



# Permute phi for A.m(Red-NR) vs C.m(Red-NR)
A.m.Red.vs.NR <- data.frame(FlyBaseID = A.m.Red.Phi$FlyBaseID,
                            Phi.Red = A.m.Red.Phi$phi,
                            Phi.NR = A.m.NR.Phi$phi)
A.m.Red.vs.NR$diff.Am <- A.m.Red.vs.NR$Phi.Red - A.m.Red.vs.NR$Phi.NR
head(A.m.Red.vs.NR)

C.m.Red.vs.NR <- data.frame(FlyBaseID = C.m.Red.Phi$FlyBaseID,
                            Phi.Red = C.m.Red.Phi$phi,
                            Phi.NR = C.m.NR.Phi$phi)
C.m.Red.vs.NR$diff.Cm <- C.m.Red.vs.NR$Phi.Red - C.m.Red.vs.NR$Phi.NR
head(C.m.Red.vs.NR)

test <- merge(A.m.Red.vs.NR, C.m.Red.vs.NR, by = "FlyBaseID")
test <- na.omit(test)
test <- test[test$FlyBaseID %in% subset.sss,]

# run the R script "boot_permute.R" for permutation functions
source("boot_permute.R")
PairedTwoPerm(test, "diff.Am", "diff.Cm")
########

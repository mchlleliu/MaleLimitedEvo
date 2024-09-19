###################################
#
#                             Grieshop et al. 2024
#                             Author: Michelle Liu
#             DsRed experimental evolution - transcriptomics analysis
#           Analysis of Red vs. NonRed Differentially Spliced (DS) genes
# 
# 
###################################

rm(list=ls())
setwd("~/Desktop/UofT/SSAV_RNA/")

# Packages
#########
library(readxl)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggblend)
library(ggpubr)
library(ggstatsplot)
library(stringr)
library(DESeq2)
#########

# Get chromosome locations  
source("Chromosome_df.R")

# set up external SSS data from Mishra et al. 2022
########
jseq.Mishra = read.table("JunctionSeq/SDIU_ase/JSresults/SDIU_ASEallGenes.results.txt",
                         sep = "\t", header = TRUE)
colnames(jseq.Mishra)[2]="FlyBaseID" # change column name to reflect the rest of the dataset

# ------ First for the Mishra reference 
# remove genes that were not assayed in males and females
jseq.Mishra = jseq.Mishra[!is.na(jseq.Mishra$expr_F) & !is.na(jseq.Mishra$expr_M),]
dim(jseq.Mishra)

# remove novel counting bins
#   for the purpose of comparing Mishra with our populations, we want to make sure that the counting bins
#   used are the same.
jseq.Mishra = jseq.Mishra %>% filter(!str_detect(countbinID, "N"))
dim(jseq.Mishra) # check how many novel splice sites were removed

# assign genes with significant sex-specific splicing
jseq.Mishra <- jseq.Mishra %>% mutate(SSS = ifelse(geneWisePadj < 0.01, TRUE, FALSE))

# list of genes significant and non-significant
Mishra.sig.SSS <- unique(jseq.Mishra$FlyBaseID[jseq.Mishra$geneWisePadj < 0.01])
Mishra.nonsig.SSS <-  unique(jseq.Mishra$FlyBaseID[jseq.Mishra$geneWisePadj >= 0.01])
length(Mishra.sig.SSS)
length(Mishra.nonsig.SSS)

# save list of SSS genes from Mishra population
# write_delim(data.frame(Mishra.sig.SSS), file = "JunctionSeq/SDIU_ase/JSresults/Mishra.sig.SSS_genes.txt", delim = ",")


# Sanity check
# fisher's test for genes identified as SSS in Osada et al. and Mishra et al. 2022 population
# they should overlap significantly.

# get the Osada data analysed in Singh & Agrawal 2023
Osada <- read.csv(file="~/Desktop/UofT/SSAV_RNA/Data/SBGEandSSSdataForMBE.csv", sep=",", header=TRUE)
colnames(Osada)[2] <- "FlyBaseID"

Osada.sig.sss <- Osada$FlyBaseID[!is.na(Osada$SDIU.body.sig) & Osada$SDIU.body.sig]
Osada.nonsig.sss <- Osada$FlyBaseID[!is.na(Osada$SDIU.body.sig) & !Osada$SDIU.body.sig]
length(Osada.sig.sss)
length(Osada.nonsig.sss)

# make a column for the Mishra SSS status
test <- Osada %>% mutate(Mishra.SSS = ifelse(FlyBaseID %in% Mishra.sig.SSS, TRUE, FALSE))
fisher.test(test$SDIU.body.sig, test$Mishra.SSS) # check results

rm(test, Osada) # clear from ENV

########


# load JunctionSeq results from our samples
#######
# SSAV samples
jseq.A.f.geno = read.table("Results/A.f.geno.DS.results.tsv", 
                           sep = "\t", header = TRUE)
jseq.A.m.geno = read.table("Results/A.m.geno.DS.results.tsv", 
                           sep = "\t", header = TRUE)
jseq.C.m.geno = read.table("Results/C.m.geno.DS.results.tsv", 
                           sep = "\t", header = TRUE)
colnames(jseq.A.f.geno)[2]="FlyBaseID" 
colnames(jseq.A.m.geno)[2]="FlyBaseID" 
colnames(jseq.C.m.geno)[2]="FlyBaseID" 

# remove genes close to the DsRed marker
DsRed_genes <- read.delim(file="Data/dmel_2R_DsRed_ids.tsv", header=FALSE)
jseq.A.f.geno <- jseq.A.f.geno[!(jseq.A.f.geno$FlyBaseID %in% DsRed_genes$V1),]
jseq.A.m.geno <- jseq.A.m.geno[!(jseq.A.m.geno$FlyBaseID %in% DsRed_genes$V1),]
jseq.C.m.geno <- jseq.C.m.geno[!(jseq.C.m.geno$FlyBaseID %in% DsRed_genes$V1),]

# remove genes on the Y and on Chr 4
jseq.A.f.geno <- merge(jseq.A.f.geno, Chrs, by = "FlyBaseID")
jseq.A.m.geno <- merge(jseq.A.m.geno, Chrs, by = "FlyBaseID")
jseq.C.m.geno <- merge(jseq.C.m.geno, Chrs, by = "FlyBaseID")

jseq.A.f.geno <- jseq.A.f.geno[jseq.A.f.geno$Chr!= "Y" & jseq.A.f.geno$Chr!= "4",]
jseq.A.m.geno <- jseq.A.m.geno[jseq.A.m.geno$Chr!= "Y" & jseq.A.m.geno$Chr!= "4",]
jseq.C.m.geno <- jseq.C.m.geno[jseq.C.m.geno$Chr!= "Y" & jseq.C.m.geno$Chr!= "4",]

# make a list of samples. This is just so we can automate 
# it easier in subsequent steps bcs i don't want to write everything thrice T-T.
SSAV.sample.types <- data.frame(sampleType = c("A.m", "A.f", "C.m")) %>% 
  mutate(raw.JS.table = paste0("jseq.",sampleType,".geno"))
#######


# SSAV sample types data table:
SSAV.sample.types # check the format




# prepare data: filtering un-assayed genes, assign significance
#######

# ------ For SSAV samples
# set FDRThreshold to a less stringent one for the SSAV data
FDRThreshold = 0.1
mappedReadsThreshold = 50

for(i in SSAV.sample.types$raw.JS.table){
  tmp.JS.table <- get(paste0(i))
  # remove untested genes
  tmp.JS.table = tmp.JS.table[!is.na(tmp.JS.table$expr_Red) & !is.na(tmp.JS.table$expr_NR),]
  tmp.JS.table = tmp.JS.table[!is.na(tmp.JS.table$geneWisePadj),]
  
  print(dim(tmp.JS.table)) # check number of counting bins
  tmp.JS.table = tmp.JS.table %>% filter(!str_detect(countbinID, "N")) # remove novel splice sites
  print(dim(tmp.JS.table)) # check how many novel splice sites were removed
  
  # assign gene-wise significance 
  tmp.JS.table <- tmp.JS.table %>% mutate(sig.hit = ifelse(geneWisePadj <= FDRThreshold, TRUE, FALSE))
  
  assign(i, tmp.JS.table) 
  rm(tmp.JS.table)
}


# do not assign as significantly spliced if the controls also differ
jseq.A.f.geno <- jseq.A.f.geno %>% mutate(sig.hit = ifelse(!(FlyBaseID %in% jseq.C.m.geno[jseq.C.m.geno$sig.hit,]$FlyBaseID),
                                                           sig.hit, FALSE))
length(unique(jseq.A.f.geno[!is.na(jseq.A.f.geno$sig.hit) & jseq.A.f.geno$sig.hit,]$FlyBaseID))

jseq.A.m.geno <- jseq.A.m.geno %>% mutate(sig.hit = ifelse(!(FlyBaseID %in% jseq.C.m.geno[jseq.C.m.geno$sig.hit,]$FlyBaseID),
                                                           sig.hit, FALSE))
length(unique(jseq.A.m.geno[!is.na(jseq.A.m.geno$sig.hit) & jseq.A.m.geno$sig.hit,]$FlyBaseID))


# combine male and female results
jseq.All.geno <- merge(jseq.A.m.geno[,c("FlyBaseID", "geneWisePadj", "sig.hit")], 
                       jseq.A.f.geno[,c("FlyBaseID", "geneWisePadj", "sig.hit")], 
                       by = "FlyBaseID", all = T) 
# if significant in males and/or females, set as TRUE
jseq.All.geno <- jseq.All.geno[!duplicated(jseq.All.geno$FlyBaseID),] %>%
  dplyr::mutate(sig.hit = ifelse(!is.na(sig.hit.x) & sig.hit.x, TRUE, 
                                 ifelse(!is.na(sig.hit.y) & sig.hit.y, TRUE, 
                                        ifelse(is.na(sig.hit.y) & is.na(sig.hit.x), NA, FALSE))))


jseq.All.geno <-  jseq.All.geno[,c("FlyBaseID", "geneWisePadj.x", "sig.hit.x", 
                                   "geneWisePadj.y", "sig.hit.y", "sig.hit")]
# jseq.All.geno <- unique(na.omit(jseq.All.geno))
colnames(jseq.All.geno) = c("FlyBaseID", "geneWisePadj.M", "sig.hit.M", 
                            "geneWisePadj.F", "sig.hit.F", "sig.hit")
# save table
# write.table(jseq.All.geno, file="Results/jseq.All.geno.txt", quote=F, sep = "\t")

# # list of all candidate genes
# all.candidates <- unique(c(jseq.All.geno$FlyBaseID[jseq.All.geno$sig.hit], SSAV.geno$FlyBaseID[SSAV.geno$Sig]))
# # length(all.candidates)
# # write.table(all.candidates, file="Results/DE_AS_candidate.list.txt", quote=F, sep = "\t")

# add all sample dataset to the list of sample types
SSAV.sample.types <- rbind(SSAV.sample.types, c("A.all", "jseq.All.geno"))

#######



###### differential splicing analysis for SSAV samples
# check how many genes were spliced differently between Red and NR, 
# do some Fisher's exact tests for association with sex-specific splicing in Mishra data:

# load DE data to check for overlap in DE and DS genes ######
A.f.geno <- read.delim("Results/A.f.geno_candidates.tsv")
A.m.geno <- read.delim("Results/A.m.geno_candidates.tsv")
C.m.geno <- read.delim("Results/C.m.geno_candidates.tsv")
SSAV.geno <- read.delim("Results/All.geno_candidates.tsv")

# add chromosome position info
A.f.geno <- merge(A.f.geno, Chrs, by = "FlyBaseID")
A.m.geno <- merge(A.m.geno, Chrs, by = "FlyBaseID")
C.m.geno <- merge(C.m.geno, Chrs, by = "FlyBaseID")
SSAV.geno <- merge(SSAV.geno, Chrs, by = "FlyBaseID")

SSAV.sample.types$DEdata <- c("A.m.geno", "A.f.geno", "C.m.geno", "SSAV.geno")
#######


# generate data for Table S1 and Table S2
# initialize data frame object to store results (wow that variable name do be long)
fishers.test.RedNR.splice.results <- data.frame(sampleType = SSAV.sample.types$sampleType)

for(i in 1:dim(SSAV.sample.types)[1]){
  tmp.JS.table <- get(paste0(SSAV.sample.types[i, 2])) # get junctionSeq table for the current sample type
  tmp.JS.table <- tmp.JS.table[!is.na(tmp.JS.table$sig.hit),] # only get genes which could be assayed
  
  DE_data <- get(paste0(SSAV.sample.types[i, 3])) # get the DE data
  
  # if processing non-combined data, cut down the table so it is faster
  if(i < 4) {
    tmp.JS.table <- tmp.JS.table[,c(1,27)] # take only the FlyBaseID column and the significant DS/not column info
    tmp.JS.table = tmp.JS.table[!duplicated(tmp.JS.table$FlyBaseID),]
  }
  
  # number of all genes assayed
  fishers.test.RedNR.splice.results$N.all[i] = length(tmp.JS.table$FlyBaseID)
  # number of significantly DS genes
  fishers.test.RedNR.splice.results$N.sig[i] = length(tmp.JS.table$FlyBaseID[!is.na(tmp.JS.table$sig.hit) 
                                                                             & tmp.JS.table$sig.hit])
  
  # assign overlap with SSS genes as defined with the Mishra et al. 2022 data set
  tmp.JS.table = tmp.JS.table %>% 
    dplyr::mutate(Mishra.SSS = ifelse(FlyBaseID %in% jseq.Mishra[!is.na(jseq.Mishra$SSS) & 
                                                                   jseq.Mishra$SSS,]$FlyBaseID, 
                                      TRUE, ifelse(is.na(jseq.Mishra$SSS), NA, FALSE)))
  
  
  # any overlap with differentially expressed genes?
  tmp.JS.table = tmp.JS.table %>% 
    dplyr::mutate(DE.Sig = ifelse(FlyBaseID %in% DE_data[!is.na(DE_data$Sig) 
                                                         & DE_data$Sig,]$FlyBaseID, 
                                  TRUE, ifelse(is.na(DE_data$Sig), NA, FALSE)))
  
  # number of genes assayed with sig SSS status defined from Mishra et al. 2022 datset
  fishers.test.RedNR.splice.results$N.sig.Mishra[i] <- length(tmp.JS.table$FlyBaseID[!is.na(tmp.JS.table$Mishra.SSS) &
                                                                                       tmp.JS.table$Mishra.SSS])
  # number of genes assayed with non sig SSS status defined from Mishra et al. 2022 datset
  fishers.test.RedNR.splice.results$N.nonsig.Mishra[i] <- length(tmp.JS.table$FlyBaseID[!is.na(tmp.JS.table$Mishra.SSS) &
                                                                                          !tmp.JS.table$Mishra.SSS])
  
  # number of genes assayed that also has sig SSS status defined from Mishra et al. 2022 datset
  fishers.test.RedNR.splice.results$N.DS.nonsig.Mishra[i] <- length(tmp.JS.table$FlyBaseID[!is.na(tmp.JS.table$Mishra.SSS) &
                                                                                             !tmp.JS.table$Mishra.SSS &
                                                                                             !is.na(tmp.JS.table$sig.hit) &
                                                                                             tmp.JS.table$sig.hit ])
  # do fisher test for overlap between SSS and DS genes
  if(sum(tmp.JS.table$sig.hit) > 0 & sum(tmp.JS.table$Mishra.SSS) > 0){
    fishers.test.RedNR.splice.results$SSS.DS.pval[i] <- fisher.test(tmp.JS.table$sig.hit, tmp.JS.table$Mishra.SSS)$p.val
    
  }
  
  
  # number of genes assayed with significant DE status
  fishers.test.RedNR.splice.results$N.sig.DE[i] <- length(tmp.JS.table$FlyBaseID[!is.na(tmp.JS.table$DE.Sig) &
                                                                                   tmp.JS.table$DE.Sig])
  # number of genes assayed with significant DE AND DS status
  fishers.test.RedNR.splice.results$N.DS.sig.DE[i] <- length(tmp.JS.table$FlyBaseID[!is.na(tmp.JS.table$DE.Sig) &
                                                                                      tmp.JS.table$DE.Sig &
                                                                                      !is.na(tmp.JS.table$sig.hit) &
                                                                                      tmp.JS.table$sig.hit ])
  # number of genes assayed with non significant DE status
  fishers.test.RedNR.splice.results$N.nonsig.DE[i] <- length(tmp.JS.table$FlyBaseID[!is.na(tmp.JS.table$DE.Sig) &
                                                                                      !tmp.JS.table$DE.Sig])
  # number of genes assayed with non significant DE status but is significantly DS
  fishers.test.RedNR.splice.results$N.DS.nonsig.DE[i] <- length(tmp.JS.table$FlyBaseID[!is.na(tmp.JS.table$DE.Sig) &
                                                                                         !tmp.JS.table$DE.Sig &
                                                                                         !is.na(tmp.JS.table$sig.hit) &
                                                                                         tmp.JS.table$sig.hit ])
  # do fisher test for overlap between DE and DS genes
  if(sum(tmp.JS.table$sig.hit) > 0 & sum(tmp.JS.table$DE.Sig) > 0){
    fishers.test.RedNR.splice.results$DE.DS.pval[i] <- fisher.test(tmp.JS.table$sig.hit, tmp.JS.table$DE.Sig)$p.val
  }
  
  assign(paste0(SSAV.sample.types[i, 1],".tmp.fisher"), tmp.JS.table)
}

fishers.test.RedNR.splice.results # print table of results

# # save table
# write.table(fishers.test.RedNR.splice.results,
#             file = "Results/splicing.Red.NR.fisher.tests_filtered.csv", 
#             sep = ",", quote = FALSE, row.names = F)
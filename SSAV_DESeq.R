###################################
#
#                             Grieshop et al. 2023
#             DsRed experimental evolution - transcriptomics analysis
#             Data were prepared using file DsRed_transcriptomics.txt
# 
# 
###################################

rm(list=ls()) # Clears the environment

# set path to data files
Data_path <- "~/Desktop/UofT/PRJ1/ReadCount_KG_renamed/"
setwd(Data_path)

# Packages
##########
library(DESeq2)
library(apeglm)
library(ashr)
library(vsn)
library(gridExtra)
library(plyr)
library(dplyr)
library(car)
library(Rmisc)
library(Hmisc)
library(mltools)
library(boot)
library(clipr) 
library(broom)
library(tidyr)

##########


# make dataframe for DESeq2
##########

# list the files by category
# the ordering is based on treatment(A/C)_rep(1:6)_sex(F/M)
files <- list.files(Data_path)

# chop off the ".tsv" from those file names
samplename <- gsub('.{4}$', '', files)

# set up factor for condition "sex" (must match order that files are listed; see "samplename" and "files")
sex <- factor(c(rep(c("Female", "Male"), each = 2, times = 6), # Af, Am
                rep("Male", times = 12)))  # Cm

# set up factor for each genotype (Red or NR) "geno" (must match order that files are listed; see "samplename" and "files")
geno <- factor(c(rep(c("NR", "Red"), each = 1, times = 18)))

# set up factor for each treatment (A or C) "trt" (must match order that files are listed; see "samplename" and "files")
trt <- factor(c(rep("A", 24), 
                rep("C", 12)))

# set up factor for each treatment (Af, Am, or Cm) "trt2" 
trt2 <- factor(c(rep(c("Af", "Am"), each = 2, times = 6), 
                 rep("Cm", 12)))

# set up factor for batch "rep" (must match order that files are listed; see "samplename" and "files")
rep <- factor(c(as.character(c(rep(1:6, each = 4, times = 1), # 6 replicates for each of the 4 sex-geno in A
                               rep(1:6, each = 2, times = 1)))))  # 6 replicates for each of the male-geno in C


# put those into a dataframe 
sampleTable <- data.frame(sampleName = samplename, 
                          fileName = files, 
                          sex = sex,
                          geno = geno,
                          trt = trt,
                          trt2 = trt2,
                          rep = rep)

str(sampleTable) # Check that factors are factors, else use e.g.: sampleTable$rep <- factor(sampleTable$rep)

##########


# Split data frame up into different subsets 
##########

# trt A data (list of sample names, and files)
A <- sampleTable[(sampleTable$trt == "A"),] 
A$trt <- droplevels(A$trt)

# A.f data
A.f <- A[(A$sex == "Female"),]
A.f$sex <- droplevels(A.f$sex)

# A.m data
A.m <- A[(A$sex == "Male"),]
A.m$sex <- droplevels(A.m$sex)

# C.m
C.m <- sampleTable[(sampleTable$trt == "C"),] # There are only Male samples for C
C.m$trt <- droplevels(C.m$trt)

# Red.m
Red.m <- Males[(Males$geno == "Red"),]
Red.m$geno <- droplevels(Red.m$geno)

# NR.m
NR.m <- Males[(Males$geno == "NR"),]
NR.m$geno <- droplevels(NR.m$geno)

# A.f.C.m (contrast A females and C males)
A.f.C.m <- rbind(A.f, C.m) # For contrasting A.f to C.m
A.f.C.m_Red <- A.f.C.m[(A.f.C.m$geno == "Red"),] # within Red
A.f.C.m_Red$geno <- droplevels(A.f.C.m_Red$geno)
A.f.C.m_NR <- A.f.C.m[(A.f.C.m$geno == "NR"),] # within NR
A.f.C.m_NR$geno <- droplevels(A.f.C.m_NR$geno)

# A.f.nr_A.m.r (contrast A nonRed females and A red males)
A.f.nr <- A.f[(A.f$geno == "NR"),]
A.m.r <- A.m[(A.m$geno == "Red"),]
A.f.nr_A.m.r <- rbind(A.f.nr, A.m.r)

# A.Red
A.Red <- A[(A$geno == "Red"),]
A.Red$geno <- droplevels(A.Red$geno)

# A.NR
A.NR <- A[(A$geno == "NR"),]
A.NR$geno <- droplevels(A.NR$geno)

# Just Males (A and C)
Males <- sampleTable[(sampleTable$sex == "Male"),] 
Males$sex <- droplevels(Males$sex)

# Just Females (A)
Females <- sampleTable[(sampleTable$sex == "Female"),] 
Females$sex <- droplevels(Females$sex)

# A.Red.m
A.Red.m <- Males[(Males$geno == "Red") & Males$trt == "A",]
A.Red.m$geno <- droplevels(A.Red.m$geno)

# A.NR.m
A.NR.m <- Males[(Males$geno == "NR") & Males$trt == "A",]
A.NR.m$geno <- droplevels(A.NR.m$geno)

# A.Red.f
A.Red.f <- Females[(Females$geno == "Red"),]
A.Red.f$geno <- droplevels(A.Red.f$geno)

# A.NR.f
A.NR.f <- Females[(Females$geno == "NR"),]
A.NR.f$geno <- droplevels(A.NR.f$geno)

# A.Red.m.NR.f
A.Red.m.NR.f <- rbind(A.Red.m, A.NR.f)
A.Red.m.NR.f$geno <- droplevels(A.Red.m.NR.f$geno)

# A.NR.m.Red.f
A.NR.m.Red.f <- rbind(A.NR.m, A.Red.f)
A.NR.m.Red.f$geno <- droplevels(A.NR.m.Red.f$geno)

##########


# Set up the different contrast designs
# (sets up data and run DESeq2 with all the different comparisons that you want)
##########

# All data - for first PCA look (A.REDvNR)
dds.all <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                      design = ~ rep + trt2 + geno)

# For geno differences in A.f
dds.A.f.geno <- DESeqDataSetFromHTSeqCount(sampleTable = A.f, 
                                           design = ~ rep + geno) 

dds.A.f.geno_rep <- DESeqDataSetFromHTSeqCount(sampleTable = A.f, 
                                           design = ~ geno + rep) 

# For geno differences in A.m
dds.A.m.geno <- DESeqDataSetFromHTSeqCount(sampleTable = A.m, 
                                           design = ~ rep + geno) 

# For geno differences in C.m
dds.C.m.geno <- DESeqDataSetFromHTSeqCount(sampleTable = C.m, 
                                           design = ~ rep + geno) 

# For trt differences in NR.m
dds.NR.m.trt <- DESeqDataSetFromHTSeqCount(sampleTable = NR.m,
                                           design = ~ rep + trt)

# For trt differences in Red.m
dds.Red.m.trt <- DESeqDataSetFromHTSeqCount(sampleTable = Red.m,
                                            design = ~ rep + trt)

# Males
dds.Males <- DESeqDataSetFromHTSeqCount(sampleTable = Males,
                                        design = ~ rep + geno + trt) 

# A.f versus C.m within Red
dds.A.f.C.m_Red <- DESeqDataSetFromHTSeqCount(sampleTable = A.f.C.m_Red,
                                              design = ~ rep + sex) 

# A.f versus C.m within Red
dds.A.f.nr_A.m.r <- DESeqDataSetFromHTSeqCount(sampleTable = A.f.nr_A.m.r,
                                               design = ~ rep + sex) 

# A.f versus C.m within NR
dds.A.f.C.m_NR <- DESeqDataSetFromHTSeqCount(sampleTable = A.f.C.m_NR,
                                             design = ~ rep + sex) 

# Interaction
dds.Int <- DESeqDataSetFromHTSeqCount(sampleTable = Males,
                                      design = ~ rep + geno + trt + geno:trt) 

# For sex differences in A.Red
dds.A.Red.sex <- DESeqDataSetFromHTSeqCount(sampleTable = A.Red,
                                            design = ~ rep + sex)

dds.sex <- DESeqDataSetFromHTSeqCount(sampleTable = A,
                                      design = ~ rep + sex)

# For sex differences in A.NR
dds.A.NR.sex <- DESeqDataSetFromHTSeqCount(sampleTable = A.NR,
                                           design = ~ rep + sex)

# For SBGE in natural pair (Red males and NR females)
dds.A.Red.m.NR.f.sex <- DESeqDataSetFromHTSeqCount(sampleTable = A.Red.m.NR.f,
                                                   design = ~ rep + sex)

# For SBGE in unnatural pair (Red females and NR males)
dds.A.NR.m.Red.f.sex <- DESeqDataSetFromHTSeqCount(sampleTable = A.NR.m.Red.f,
                                                   design = ~ rep + sex)

##########


# \\||// #
# Set parameters and details for the focal contrast (modify/rerun for each desired contrast)
##########

# Adjust filtering criteria
minCountPerSample = 1 #min read count per sample
minAvgPerCat = 10 #min average read per category

# Specify the contrast
focal.contrast <-  dds.sex # Change accordingly

# Specify the samples for each category of the focal contrasts (some tricky ones commented here)
# this is based on the order of the sample names

# A.REDvNR
# numerator <- samplename[seq(2, 24, by = 2)] #A-Red
# denominator <- samplename[seq(1, 24, by = 2)] #A-NR

# dds.A.m.geno
# numerator <- samplename[seq(4, 24, by = 4)] # A.red.males
# denominator <- samplename[seq(3, 24, by = 4)] # A.nr.males

# dds.A.f.geno
# numerator <- samplename[seq(2, 24, by = 4)] # A.red.females
# denominator <- samplename[seq(1, 24, by = 4)] # A.nr.females

# dds.C.m.geno
# numerator <- samplename[seq(26, 36, by = 2)] # C.red.males
# denominator <- samplename[seq(25, 36, by =2)] # C.nr.males

# dds.Red.m.trt 
# numerator <- samplename[seq(4, 24, by = 4)] # A.Red.m
# denominator <- samplename[seq(26, 36, by = 2)] # C.Red.m

# dds.NR.m.trt
# numerator <- samplename[seq(3, 24, by = 4)] # A.NR.m
# denominator <- samplename[seq(25, 36, by = 2)] # C.NR.m

# dds.A.Red.sex
# numerator <- samplename[seq(4, 24, by = 4)]
# denominator <-samplename[seq(2, 24, by = 4)]

# dds.A.NR.sex
# numerator <- samplename[seq(3, 24, by = 4)]
# denominator <-samplename[seq(1, 24, by = 4)]

# # dds.A.Red.m.NR.f.sex
# numerator <- samplename[seq(4, 24, by = 4)]
# denominator <- samplename[seq(1, 24, by =4)]

# dds.A.NR.m.Red.f.sex
# numerator <- samplename[seq(3, 24, by = 4)]
# denominator <- samplename[seq(2, 24, by =4)]

# dds.A.sex
numerator <- samplename[c(seq(3, 24, by = 4), seq(4, 24, by = 4))]
denominator <-samplename[c(seq(1, 24, by = 4), seq(2, 24, by =4))]

# Analysis details
# print "str(sampleTable)" to see your options here
factor.numerator.denominator = c("sex", "Male", "Female") # that is c("factor", "numerator", "denominator"); flip if desired
alpha.threshold = 0.05 # this sets alpha for the adjusted pvalue; default is 0.1.
##########


# Run the focal contrast
##########

# 1. get read counts 
countdf = DESeq2::counts(focal.contrast) #read counts of all the genes for each sample

# 2a. Create logical vector for whether there's a low count in any sample of either the numerator or denominator
# prod() returns the product of all the values present in its arguments
# TRUE if low count and FALSE if the count is above threshold
lowCountAnySample.numerator = sapply(1:(dim(countdf)[1]), function (x) prod(countdf[x, numerator]) < minCountPerSample)  
lowCountAnySample.denominator = sapply(1:(dim(countdf)[1]), function (x) prod(countdf[x, denominator]) < minCountPerSample)  
lowCountAnySample.Either =  lowCountAnySample.numerator | lowCountAnySample.denominator

# 2b. Create logical vector for whether there's a good average count in both numerator and denominator
# calculate the mean reads of each gene across the samples
avgCounts.numerator = rowMeans(countdf[, numerator]) 
avgCounts.denominator = rowMeans(countdf[, denominator]) 
# create the logical vectors
goodAvgCount.numerator = avgCounts.numerator > minAvgPerCat
goodAvgCount.denominator = avgCounts.denominator > minAvgPerCat
goodAvgCount.Both = goodAvgCount.numerator & goodAvgCount.denominator

# 2c. Create logical vector indicating fulfillment of goodAvgCount.Both and NOT having lowCountAnySample.Either
# remove any genes that do not have avg. read count > threshold or any genes with low read count for any sample
keep.these = goodAvgCount.Both & (!lowCountAnySample.Either) 

# 2.d Inspect what you're about to discard
sum(lowCountAnySample.numerator) # number of genes with low read count for the numerator samples
sum(lowCountAnySample.denominator) # number of genes with low read count for the denominator samples
sum(lowCountAnySample.numerator & lowCountAnySample.denominator) # number that are the same genes
mean(avgCounts.numerator) # mean before filtering
mean(avgCounts.denominator) # mean before filtering
mean(avgCounts.numerator)/mean(avgCounts.denominator) # ratio of those means
# many other things you could do here to know what you're tossing out 

# 3. Filter according to the above
focal.contrast.filtered = focal.contrast[keep.these]

# 4. Do the analysis for that focal contrast using those filtered data
DESeq.Analysis = DESeq(focal.contrast.filtered) # Change to "focal.contrast" to bypass all filtering.

# 5a. Get the results of that analysis
DESeq.Results = results(DESeq.Analysis,
                        contrast = factor.numerator.denominator, # that is c("factor", "numerator", "denominator"); 
                        # by defult, it will take the last term in the design (e.g. the interaction term, if that's last)
                        alpha = alpha.threshold, independentFiltering=T) 

# 5b. Remove Y genes
# read in .tsv files of X and Y geneIDs (see end of DsRed "transcriptomics.txt")
Ychr <- read.delim("~/Desktop/UofT/PRJ1/data/Y.chromosome.genes.tsv", sep = '\t', header = TRUE)
# Flag genes near DsRed marker
DsRed_genes <- read.delim(file="~/Desktop/UofT/SSAV_RNA/Data/dmel_2R_DsRed_ids.tsv", header=FALSE)
# Make the FlyBaseIDs into rownames and index by this new column
DESeq.Results$FlyBaseID = rownames(DESeq.Results)
DESeq.Results <- DESeq.Results[!(DESeq.Results$FlyBaseID %in% Ychr$geneID),]
DESeq.Results <- DESeq.Results[!DESeq.Results$FlyBaseID %in% DsRed_genes$V1,]

# 6. Look at the metadata for DEseq's independent filtering 
plot(metadata(DESeq.Results)$filterNumRej, 
     type="b", ylab="number of rejections",
     xlab="quantiles of filter")
lines(metadata(DESeq.Results)$lo.fit, col="red")
abline(v=metadata(DESeq.Results)$filterTheta)

# 7. Have a basic look 
summary(DESeq.Results)
DESeq2::plotMA(DESeq.Results, ylim=c(-5,5), colSig = "red", colNonSig = "lightgray")


##########


# PCA - modify/use this to explore global differences among samples (before and after filtering)
##########

# check STATS notes to understand VST
# Swap out different dds.* designs from above or use focal.contrast or focal.contrast.filtered
vsd <- vst(focal.contrast, blind=FALSE) # Use focal.contrast and focal.contrast.filtered to see the effect of filtering on the data
# head(assay(vsd), 3)

# Then swap out different 'intgroup' variables as needed...
pcaData <- plotPCA(vsd, intgroup=c("geno", "rep"), returnData=TRUE) # the argument 'intgroup' should specify columns of colData(dds.*)
percentVar <- round(100 * attr(pcaData, "percentVar"))
# ...and 'color'/'shape' variables as needed...
ggplot(pcaData, aes(PC1, PC2, color=rep, shape=geno)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()

# pcaData[order(pcaData$PC1),]
##########


# \\||// #
# change se_trt/geno and exp_trt/geno accordingly
# Organize results into data frame (may need to adjust some column names)
##########

# Make the rownames a column in the df for indexing later
DESeq.Results$FlyBaseID = rownames(DESeq.Results) # fix accordingly (e.g. FlyBaseID)


# Grab the desired variables from the DESeq contrasts above  
Results.df <- data.frame(cbind(DESeq.Results$log2FoldChange, 
                               DESeq.Results$lfcSE, 
                               DESeq.Results$padj,
                               DESeq.Results$FlyBaseID # fix accordingly (e.g. FlyBaseID)
))
colnames(Results.df) <- c("exp_sex", "se_sex", "padj", "FlyBaseID") # fix accordingly (e.g. maybe c()[1] == "exp_trt", and/or c()[4] == "FlyBaseID")
str(Results.df)
Results.df$exp_sex <- as.numeric(Results.df$exp_sex) # fix accordingly
Results.df$se_sex <- as.numeric(Results.df$se_sex) # fix accordingly
Results.df$padj <- as.numeric(Results.df$padj)
str(Results.df)

## change this to the correct file name
write.table(Results.df, file = "~/Desktop/UofT/SSAV_RNA/Results/A.m.geno_raw.tsv", sep = "\t", # Fix file name accordingly
            row.names = FALSE, col.names = TRUE)
##########


# check correlation between log2FC in SSAV (A) and in ASE
test <- merge(Results.df, ASE, by = "FlyBaseID")
cor.test(test$exp_sex, test$exp_SBGE_ase)

check_test <- merge(Results.df,SDIU, by = "FlyBaseID")
ggplot(check_test) + 
  geom_density(aes(exp_geno), fill = "grey", color = "grey") +
  labs(x=expression("Log2FC sex-bias"),y="Number of loci") +
  theme_classic() +
  theme(axis.text.x = element_text(size=18, margin = margin(5,0,5,0), color = "black"),
        axis.text.y = element_text(size=12, margin = margin(0,5,0,0), color = "black"),
        axis.title.x = element_text(size=24, margin = margin(10,0,10,0), color = "black"),
        axis.title.y = element_text(size=15, margin = margin(0,10,0,0), color = "black"),
        plot.margin = margin(6,6,6,6)) 
cor.test(check_test$exp_SBGE_ase, check_test$Whole.SBGE.Osada)
plot_corr(check_test,"exp_SBGE_ase", "Whole.SBGE.Osada", colx = "purple", coly = "blue2", 
          colNonCon = "grey", xlab = "ASE", ylab = "Osada", lim = 15, title = "")

# \\||// #
# Assign significant genes
##########
# If not already in environment, read in df to specify the results df to be used 
# C.m.geno <- read.delim("Results/C.m.geno_candidates.tsv")
# Results.df <- read.delim(file = "~/Desktop/UofT/SSAV_RNA/Results/NR.m.trt_raw.tsv", header = TRUE)

# Function to sort significant and non-significant genes by adding logical column
assign_sig <- function(contrast_df){
  contrast_df <- na.omit(contrast_df)
  contrast_df$Sig = FALSE
  for(i in 1:nrow(contrast_df)){
    if(contrast_df$padj[i] < alpha.threshold & 
       !contrast_df$FlyBaseID[i] %in% C.m.geno$FlyBaseID[C.m.geno$Sig]){
      contrast_df$Sig[i] = TRUE
    }
  }
  print(dim(contrast_df[contrast_df$Sig == TRUE, ]))
  return(contrast_df)
}

# note: need to rm() observations with padj == NA
# change df name as needed
Results.df <- assign_sig(Results.df)
dim(Results.df[Results.df$Sig == TRUE, ]) # right number?


# this is only for A.m.geno
# Results.df <- Results.df[order(Results.df$padj),]
# Results.df <- Results.df[order(abs(Results.df$exp_geno), decreasing = T),]
# colnames(Results.df)[colnames(Results.df) == "Sig"] = "Top.Sig"
# Results.df$Sig <- FALSE
# Results.df[1:350,]$Sig <- TRUE
# droplevels(Results.df)

# get candidates for males, but exclude mtDNA genes and Chr4
Results.df <- Results.df[order(abs(Results.df$exp_geno), decreasing = T),]
colnames(Results.df)[colnames(Results.df) == "Sig"] = "Top.Sig"
Results.df$Sig <- FALSE
i = 1
while(dim(Results.df[Results.df$Sig,])[1] < 350){
  if(!Results.df$FlyBaseID[i] %in% Chrs$FlyBaseID[Chrs$Chr == "Y"] & 
     Results.df$FlyBaseID[i] %in% Chrs$FlyBaseID){
  Results.df$Sig[i] <- TRUE
  }
  i = i + 1
}
droplevels(Results.df)
dim(Results.df[Results.df$Sig,])

## change this to the correct file name
write.table(Results.df, file = "~/Desktop/UofT/SSAV_RNA/Results/A.m.geno_candidates.tsv", sep = "\t", # Fix file name accordingly
            row.names = FALSE, col.names = TRUE)
##########


# Combine SSAV males and females
########
setwd("~/Desktop/UofT/SSAV_RNA/")
A.f.geno <- read.delim("Results/A.f.geno_candidates.tsv")
A.m.geno <- read.delim("Results/A.m.geno_candidates.tsv")
C.m.geno <- read.delim("Results/C.m.geno_candidates.tsv")

SSAV.geno <- merge(A.m.geno, A.f.geno, by = "FlyBaseID", all = TRUE)
colnames(SSAV.geno) <- c("FlyBaseID", "A.m.exp_geno", "A.m.se_geno", "A.m.padj", "A.m.TopSig", "A.m.Sig",
                          "A.f.exp_geno", "A.f.se_geno", "A.f.padj", "A.f.Sig")
str(SSAV.geno)
# column denotes genes that are candidates in males or females
SSAV.geno <- SSAV.geno %>% mutate(Sig = ifelse(!is.na(A.m.Sig) & A.m.Sig, TRUE,
                                               ifelse(!is.na(A.f.Sig) & A.f.Sig, TRUE, 
                                                      ifelse(is.na(A.m.Sig) & is.na(A.f.Sig), NA, FALSE)))) 
SSAV.geno <- SSAV.geno[!is.na(SSAV.geno$Sig),]
dim(SSAV.geno[SSAV.geno$Sig,])

# only keep concordant changes
SSAV.geno.con <- na.omit(SSAV.geno[(SSAV.geno$A.f.exp_geno > 0 & SSAV.geno$A.m.exp_geno > 0) |
                             (SSAV.geno$A.f.exp_geno < 0 & SSAV.geno$A.m.exp_geno < 0),])

## save to file
write.table(SSAV.geno, file = "~/Desktop/UofT/SSAV_RNA/Results/All.geno_candidates.tsv", sep = "\t", # Fix file name accordingly
            row.names = FALSE, col.names = TRUE)
########


########
# candidates around DsRed marker excluded?
dim(SSAV.geno[SSAV.geno$FlyBaseID %in% DsRed_genes$V1 & SSAV.geno$Sig,])
dim(A.m.geno[A.m.geno$FlyBaseID %in% DsRed_genes$V1 & A.m.geno$Sig,])
dim(A.m.geno[A.m.geno$FlyBaseID %in% DsRed_genes$V1 & A.m.geno$Top.Sig,])
dim(A.f.geno[A.f.geno$FlyBaseID %in% DsRed_genes$V1 & A.f.geno$Sig,])
# total = 18/592
# male candidates = 6/200
# male significant genes = 1/16
# female candidates = 13/403

########

check_test <- merge(Results.df, ASE, by = "FlyBaseID")
hist(check_test$exp_SBGE_ase, breaks = 100)
cor.test(check_test$exp_geno.x, check_test$exp_geno.y)
View(check_test[check_test$Sig.x != check_test$Sig.y,])

A.f.countfiles.test <- data.frame(countdf[]) %>% 
  tibble::rownames_to_column(var="FlyBaseID") %>%
  dplyr::mutate(count.All = rowMeans(select_if(., is.numeric), na.rm = T))

red.tmp <- data.frame(countdf[,str_detect(colnames(countdf), "_Red")]) %>% 
  tibble::rownames_to_column(var="FlyBaseID") %>%
  dplyr::mutate(count.All = rowMeans(select_if(., is.numeric), na.rm = T))

nr.tmp <- data.frame(countdf[,str_detect(colnames(countdf), "_NR")]) %>% 
  tibble::rownames_to_column(var="FlyBaseID") %>%
  dplyr::mutate(count.All = rowMeans(select_if(., is.numeric), na.rm = T))

tmp <- data.frame(cbind(countdf, lowCountAnySample.Either)) %>%
  tibble::rownames_to_column(var="FlyBaseID")
A.f.countfiles.test <- ASE[ASE$FlyBaseID %in% tmp$FlyBaseID[tmp$lowCountAnySample.Either==0],]
hist(x = A.f.countfiles.test$exp_SBGE_ase)

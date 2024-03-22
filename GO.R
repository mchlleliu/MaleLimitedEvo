###################################
#
#                             Grieshop et al. 2023
#             DsRed experimental evolution - transcriptomics analysis
#                           GO enrichement of Candidate genes
# 
# 
###################################

library(clipr)
library(dplyr)
library(readxl)

# Use gProfiler
# Or use GOrilla to set background genes
# Load results
######
A.f.geno <- read.delim("Results/A.f.geno_candidates.tsv")
A.m.geno <- read.delim("Results/A.m.geno_candidates.tsv")
SSAV.geno <- read.delim("Results/All.geno_candidates.tsv")


# list of candidates including the DE and AS analyses
all.candidates <- read.delim("Results/DE_AS_candidate.list.txt", header=T)
colnames(all.candidates) = "FlyBaseID"
str(all.candidates)

# list of all genes considered in the DE and AS analyses
all.genes <- unique(c(SSAV.geno$FlyBaseID, jseq.All.geno$FlyBaseID))
dim(SSAV.geno)
dim(jseq.All.geno)
length(all.genes)

# only keep concordant changes
SSAV.geno.con <- SSAV.geno[(SSAV.geno$A.f.exp_geno > 0 & SSAV.geno$A.m.exp_geno > 0) |
                            (SSAV.geno$A.f.exp_geno < 0 & SSAV.geno$A.m.exp_geno < 0) |
                             (is.na(SSAV.geno$A.f.exp_geno) | is.na(SSAV.geno$A.m.exp_geno)) &
                             !(is.na(SSAV.geno$A.f.exp_geno) & is.na(SSAV.geno$A.m.exp_geno)),]

# AS genes
target <- unique(c(jseq.A.f.geno$FlyBaseID[!is.na(jseq.A.f.geno$sig.hit) & jseq.A.f.geno$sig.hit]),
                 jseq.A.m.geno$FlyBaseID[!is.na(jseq.A.m.geno$sig.hit) & jseq.A.m.geno$sig.hit])
background <- unique(c(jseq.A.f.geno$FlyBaseID, jseq.A.m.geno$FlyBaseID))
background <- background[!background %in% target]
write_clip(target)
length(target)
write_clip(background)
length(background)

target <- unique(c(jseq.A.m.geno$FlyBaseID[!is.na(jseq.A.m.geno$sig.hit) & jseq.A.m.geno$sig.hit]))
background <- unique(c(jseq.A.m.geno$FlyBaseID))
background <- background[!background %in% target]
write_clip(target)
length(target)
write_clip(background)
length(background)

#######

# list of all genes vs background genes
write_clip(all.candidates$FlyBaseID)
write_clip(all.genes[!all.genes %in% all.candidates$FlyBaseID])


# copy genes to clipboard
#######
# candidates UPregulated in males and females
target <- na.omit(SSAV.geno.con[SSAV.geno.con$Sig & 
                                   (SSAV.geno.con$A.m.exp_geno > 0 & SSAV.geno.con$A.f.exp_geno > 0),]$FlyBaseID)
write_clip(target)
length(target)
background <- na.omit(SSAV.geno$FlyBaseID[!SSAV.geno$FlyBaseID %in% target])
write_clip(background)
length(background)



# candidates DOWNregulated in males and/or females
target <- na.omit(SSAV.geno.con[SSAV.geno.con$Sig & 
                                  (SSAV.geno.con$A.m.exp_geno < 0 & SSAV.geno.con$A.f.exp_geno < 0),]$FlyBaseID)
write_clip(target)
length(target)
background <- na.omit(SSAV.geno$FlyBaseID[!SSAV.geno$FlyBaseID %in% target])
write_clip(background)
length(background)



# candidates UPregulated in females
target <- A.f.geno[A.f.geno$Sig & A.f.geno$exp_geno > 0,]$FlyBaseID
background <- A.f.geno[(!A.f.geno$Sig & A.f.geno$exp_geno > 0) |
                  A.f.geno$exp_geno < 0,]$FlyBaseID
write_clip(target)
length(target)
write_clip(background)
length(background)



# DOWN
target <- A.f.geno[A.f.geno$Sig & A.f.geno$exp_geno < 0,]$FlyBaseID
background <- A.f.geno[(!A.f.geno$Sig & A.f.geno$exp_geno < 0) |
                      A.f.geno$exp_geno > 0,]$FlyBaseID
write_clip(target)
length(target)
write_clip(background)
length(background)



# candidates UPregulated in males
taregt<- A.m.geno[A.m.geno$Sig & A.m.geno$exp_geno > 0,]$FlyBaseID
background <- A.m.geno[(!A.m.geno$Sig & A.m.geno$exp_geno > 0) |
                      A.m.geno$exp_geno < 0,]$FlyBaseID
write_clip(target)
length(target)
write_clip(background)
length(background)


# DOWN
target <- (A.m.geno[A.m.geno$Sig & A.m.geno$exp_geno < 0,]$FlyBaseID)
background <- (A.m.geno[(!A.m.geno$Sig & A.m.geno$exp_geno < 0) |
                      A.m.geno$exp_geno > 0,]$FlyBaseID)
write_clip(target)
length(target)
write_clip(background)
length(background)
#######



# randomize genes
#######
randomIDs <- sample(A.m.geno$FlyBaseID, 
                    length(target), 
                    replace = TRUE)
write_clip(randomIDs)
write_clip(A.m.geno[!A.m.geno$FlyBaseID %in% randomIDs,]$FlyBaseID)
dim(A.m.geno)

#######


# Load GO analysis results
# downloaded from gProfiler 
#######
All.geno.UP_GO <- read_excel("Results/GO/GO_results.xlsx", sheet = 1)
All.geno.DOWN_GO <- read_excel("Results/GO/GO_results.xlsx", sheet = 2)
A.m.geno.UP_GO <- read_excel("Results/GO/GO_results.xlsx", sheet = 3)
A.m.geno.DOWN_GO <- read_excel("Results/GO/GO_results.xlsx", sheet = 4)
A.f.geno.UP_GO <- read_excel("Results/GO/GO_results.xlsx", sheet = 5)
A.f.geno.DOWN_GO <- read_excel("Results/GO/GO_results.xlsx", sheet = 6)

#######



pdf(file = "~/Desktop/UofT/SSAV_RNA/Plots/Enrichment_Tests/GO_DOWN_all.pdf",  # The directory you want to save the file in
    width = 5, # The width of the plot in inches
    height = 15) # The height of the plot in inches
GO_DOWN 
dev.off()


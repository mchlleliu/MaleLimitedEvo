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

# only keep concordant changes
SSAV.geno.con <- na.omit(SSAV.geno[(SSAV.geno$A.f.exp_geno > 0 & SSAV.geno$A.m.exp_geno > 0) |
                                     (SSAV.geno$A.f.exp_geno < 0 & SSAV.geno$A.m.exp_geno < 0),])

#######

# copy genes to clipboard
#######
# candidates UPregulated in males and females
write_clip(SSAV.geno.con[SSAV.geno.con$Sig & SSAV.geno.con$A.m.exp_geno > 0,]$FlyBaseID)
write_clip(SSAV.geno.con[(!SSAV.geno.con$Sig & SSAV.geno.con$A.m.exp_geno > 0) | 
                           SSAV.geno.con$A.m.exp_geno < 0,]$FlyBaseID) # Background set

length(SSAV.geno.con[SSAV.geno.con$Sig & SSAV.geno.con$A.m.exp_geno > 0,]$FlyBaseID)
length(SSAV.geno.con[(!SSAV.geno.con$Sig & SSAV.geno.con$A.m.exp_geno > 0) | 
                           SSAV.geno.con$A.m.exp_geno < 0,]$FlyBaseID) # Background set

# candidates DOWNregulated in males and females
write_clip(SSAV.geno.con[SSAV.geno.con$Sig & SSAV.geno.con$A.m.exp_geno < 0,]$FlyBaseID)
write_clip(SSAV.geno.con[(!SSAV.geno.con$Sig & SSAV.geno.con$A.m.exp_geno < 0) | 
                           SSAV.geno.con$A.m.exp_geno > 0,]$FlyBaseID) 

length(SSAV.geno.con[SSAV.geno.con$Sig & SSAV.geno.con$A.m.exp_geno < 0,]$FlyBaseID)
length(SSAV.geno.con[(!SSAV.geno.con$Sig & SSAV.geno.con$A.m.exp_geno < 0) | 
                       SSAV.geno.con$A.m.exp_geno > 0,]$FlyBaseID)

# candidates UPregulated in females
write_clip(A.f.geno[A.f.geno$Sig & A.f.geno$exp_geno > 0,]$FlyBaseID)
write_clip(A.f.geno[(!A.f.geno$Sig & A.f.geno$exp_geno > 0) |
                  A.f.geno$exp_geno < 0,]$FlyBaseID)
dim(A.f.geno)
length(A.f.geno[A.f.geno$Sig & A.f.geno$exp_geno > 0,]$FlyBaseID)
length(A.f.geno[(!A.f.geno$Sig & A.f.geno$exp_geno > 0) |
                            A.f.geno$exp_geno < 0,]$FlyBaseID)


# DOWN
write_clip(A.f.geno[A.f.geno$Sig & A.f.geno$exp_geno < 0,]$FlyBaseID)
write_clip(A.f.geno[(!A.f.geno$Sig & A.f.geno$exp_geno < 0) |
                      A.f.geno$exp_geno > 0,]$FlyBaseID)


# candidates UPregulated in males
write_clip(A.m.geno[A.m.geno$Sig & A.m.geno$exp_geno > 0,]$FlyBaseID)
write_clip(A.m.geno[(!A.m.geno$Sig & A.m.geno$exp_geno > 0) |
                      A.m.geno$exp_geno < 0,]$FlyBaseID)
# DOWN
write_clip(A.m.geno[A.m.geno$Sig & A.m.geno$exp_geno < 0,]$FlyBaseID)
length(A.m.geno[(!A.m.geno$Sig & A.m.geno$exp_geno < 0) |
                      A.m.geno$exp_geno > 0,]$FlyBaseID)
#######



# randomize genes
#######
randomIDs <- sample(SSAV.geno.con$FlyBaseID, 
                    length(SSAV.geno.con[SSAV.geno.con$Sig & SSAV.geno.con$A.m.exp_geno < 0,]$FlyBaseID), 
                    replace = TRUE)
write_clip(randomIDs)
write_clip(SSAV.geno.con[!SSAV.geno.con$FlyBaseID %in% randomIDs,]$FlyBaseID)
dim(A.m.geno)

#######


# Load GO analysis results
# downloaded from gProfiler 
#######
All.geno.UP_GO <- read_excel("Results/GO/GO_all.xlsx", sheet = 1)
All.geno.UP_GO$`P-value` <- as.numeric(All.geno.UP_GO$`P-value`)
All.geno.UP_GO$`FDR q-value` <- as.numeric(All.geno.UP_GO$`FDR q-value`)

All.geno.DOWN_GO <- read_excel("Results/GO/GO_all.xlsx", sheet = 2)
All.geno.DOWN_GO$`P-value` <- as.numeric(All.geno.DOWN_GO$`P-value`)
All.geno.DOWN_GO$`FDR q-value` <- as.numeric(All.geno.DOWN_GO$`FDR q-value`)

A.f.geno.UP_GO <- read_excel("Results/GO/GO_all.xlsx", sheet = 3)
A.f.geno.UP_GO$`P-value` <- as.numeric(A.f.geno.UP_GO$`P-value`)
A.f.geno.UP_GO$`FDR q-value` <- as.numeric(A.f.geno.UP_GO$`FDR q-value`)

A.f.geno.DOWN_GO <- read_excel("Results/GO/GO_all.xlsx", sheet = 4)
A.f.geno.DOWN_GO$`P-value` <- as.numeric(A.f.geno.DOWN_GO$`P-value`)
A.f.geno.DOWN_GO$`FDR q-value` <- as.numeric(A.f.geno.DOWN_GO$`FDR q-value`)

A.m.geno.UP_GO <- read_excel("Results/GO/GO_all.xlsx", sheet = 5)
A.m.geno.UP_GO$`P-value` <- as.numeric(A.m.geno.UP_GO$`P-value`)
A.m.geno.UP_GO$`FDR q-value` <- as.numeric(A.m.geno.UP_GO$`FDR q-value`)

A.m.geno.DOWN_GO <- read_excel("Results/GO/GO_all.xlsx", sheet = 6)
A.m.geno.DOWN_GO$`P-value` <- as.numeric(A.m.geno.DOWN_GO$`P-value`)
A.m.geno.DOWN_GO$`FDR q-value` <- as.numeric(A.m.geno.DOWN_GO$`FDR q-value`)


Fem_GO_DOWN <- ggplot(A.f.geno.DOWN_GO[A.f.geno.DOWN_GO$`P-value` < 0.00005,]) + 
  geom_bar(aes(x=`GO Term`, y=b, fill=`P-value`), stat = "identity") + 
  scale_fill_continuous(low="blue", high="red", name = "P-val") +
  theme_classic() +
  theme(legend.text = element_text(size = 20, color = "black"),
        plot.margin = margin(6,6,6,6), 
        axis.text.y = element_text(size = 10, color = "black")) +
  coord_flip() +
  labs(y="number of genes",
       x="GO Term") 
dim(All.geno.UP_GO)
#######



pdf(file = "~/Desktop/UofT/SSAV_RNA/Plots/Enrichment_Tests/GO_DOWN_all.pdf",  # The directory you want to save the file in
    width = 5, # The width of the plot in inches
    height = 15) # The height of the plot in inches
GO_DOWN 
dev.off()


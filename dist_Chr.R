###################################
#
#                             Grieshop et al. 2023
#                             Author: Michelle Liu
#             DsRed experimental evolution - transcriptomics analysis
#                  Distribution of SA candidates by chromosomes
# 
# 
###################################

rm(list=ls()) # Clears the environment
setwd("~/Desktop/UofT/SSAV_RNA/")


# packages required
#########
library(tidyr)
library(plyr)
library(dplyr)
library(broom)
library(ggplot2)
######



# Get chromosome locations  
##########

all.genes <- read.delim(file="~/Desktop/UofT/SSAV_RNA/Data/all.genes.tsv", sep="\t", header=FALSE)
colnames(all.genes) = c("FlyBaseID")

Xchr <- read.delim(file="~/Desktop/UofT/SSAV_RNA/Data/X.chromosome.genes.tsv", sep="\t", header=TRUE)
colnames(Xchr) = c("FlyBaseID")
Xchr$Chr <- rep("X", dim(Xchr)[1])

Ychr <- read.delim(file="~/Desktop/UofT/SSAV_RNA/Data/Y.chromosome.genes.tsv", sep="\t", header=TRUE)
colnames(Ychr) = c("FlyBaseID")
Ychr$Chr <- rep("Y", dim(Ychr)[1])


chr2L <- read.delim(file="~/Desktop/UofT/SSAV_RNA/Data/2L.chromosome.genes.tsv", sep="\t", header=FALSE)
colnames(chr2L) = c("FlyBaseID")
chr2L$Chr <- rep("2L", dim(chr2L)[1])
chr2R <- read.delim(file="~/Desktop/UofT/SSAV_RNA/Data/2R.chromosome.genes.tsv", sep="\t", header=FALSE)
colnames(chr2R) = c("FlyBaseID")
chr2R$Chr <- rep("2R", dim(chr2R)[1])
# 
chr2 <- rbind(chr2L, chr2R)
chr2$Chr <- rep("2", dim(chr2)[1])

chr3L <- read.delim(file="~/Desktop/UofT/SSAV_RNA/Data/3L.chromosome.genes.tsv", sep="\t", header=FALSE)
colnames(chr3L) = c("FlyBaseID")
chr3L$Chr <- rep("3L", dim(chr3L)[1])
chr3R <- read.delim(file="~/Desktop/UofT/SSAV_RNA/Data/3R.chromosome.genes.tsv", sep="\t", header=FALSE)
colnames(chr3R) = c("FlyBaseID")
chr3R$Chr <- rep("3R", dim(chr3R)[1])
#
chr3 <- rbind(chr3L, chr3R)
chr3$Chr <- rep("3", dim(chr3)[1])

Chrs <- rbind(Xchr, Ychr, chr2, chr3) # Not all genes; just X, Y, 2, and 3.
Chrs$Chr <- as.factor(Chrs$Chr)

Chrs_All <- rbind(Xchr, Ychr, chr2L, chr2R, chr3L, chr3R)
Chrs_All$Chr <- as.factor(Chrs_All$Chr)
##########

# Prepare plotting dataset
########


# load results if not loaded in env.
A.f.geno <- read.delim("Results/A.f.geno_candidates.tsv")
A.m.geno <- read.delim("Results/A.m.geno_candidates.tsv")
SSAV.geno <- read.delim("Results/All.geno_candidates.tsv")
# jseq.All.geno <- read.delim("Results/jseq.All.geno.txt")

# include Chr 
A.m.geno_Chr <- merge(A.m.geno, Chrs, by = "FlyBaseID", all = TRUE)
A.m.geno_Chr <- A.m.geno_Chr[!is.na(A.m.geno_Chr$Sig) & !is.na(A.m.geno_Chr$Chr),]
A.f.geno_Chr <- merge(A.f.geno, Chrs, by = "FlyBaseID", all = TRUE)
A.f.geno_Chr <- A.f.geno_Chr[!is.na(A.f.geno_Chr$Sig) & !is.na(A.f.geno_Chr$Chr),]

jseq.All.geno_Chr <- merge(jseq.All.geno[,c("FlyBaseID", "geneWisePadj.x", "sig.hit.x", "geneWisePadj.y", "sig.hit.y", "sig.hit")],
                           Chrs, by = "FlyBaseID", all = T)
jseq.All.geno_Chr <- unique(na.omit(jseq.All.geno_Chr))
jseq.All.geno_Chr <- jseq.All.geno_Chr[!jseq.All.geno_Chr$FlyBaseID %in% DsRed_genes$V1,]
colnames(jseq.All.geno_Chr)[6] = "Sig"

# column denotes genes that are candidates in males or females
SSAV.geno_Chr <- merge(SSAV.geno, Chrs, by = "FlyBaseID", all = TRUE)
SSAV.geno_Chr <- SSAV.geno_Chr[!is.na(SSAV.geno_Chr$Sig) & !is.na(SSAV.geno_Chr$Chr),]

sampleTable <- data.frame(sampleType = c("A.m", "A.f", "A.all", "AS"),
                          df = c("A.m.geno_Chr", "A.f.geno_Chr", "SSAV.geno_Chr", "jseq.All.geno_Chr"))

#########



## binned proportion of candidate genes
propChr <- function(dat){
  total_All <- dim(dat)[1]
  total_Sig <- dim(dat[dat$Sig,])[1]
  total_NS <- total_All - total_Sig
  
  total_Chr <- dat %>% 
    dplyr::group_by(Chr) %>%
    dplyr::summarise(Chr_count = n()) 
  
  fract <- dat %>% 
    dplyr::group_by(Chr, Sig) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::mutate(frac_Sig = ifelse(Sig, count/total_Sig, count/total_NS))
  
  fract <- merge(fract, total_Chr, by = "Chr") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(frac_by_Chr = count/Chr_count,
           lower = list(binom.test(count, Chr_count)),
           upper = lower$conf.int[2],
           lower = lower$conf.int[1])

  return(fract)
}
plotChrprop <- function(dat, xlab){
  plot_dat <- propChr(dat)
  
  prop_all <- data.frame(dim(dat[dat$Sig,])[1]/dim(dat)[1])
  prop_all$lower <- binom.test(dim(dat[dat$Sig,])[1], dim(dat)[1])$conf.int[1]
  prop_all$upper <- binom.test(dim(dat[dat$Sig,])[1], dim(dat)[1])$conf.int[2]
  
  plot_vec <- ggplot(plot_dat[plot_dat$Sig,], aes_string(x = "Chr", y = "frac_by_Chr")) + 
    geom_errorbar(ymin = plot_dat[plot_dat$Sig,]$lower, ymax = plot_dat[plot_dat$Sig,]$upper,
                  width = 0.5, size = 0.75) +
    geom_point(fill = "forestgreen", color = "forestgreen", size = 7) + 
    labs(x = print(xlab), # "omegaA_MK" = expression(italic("\u03c9A")[MK]); "alpha_MK" = expression(italic("\u03b1")[MK])
         y = "fraction of candidate genes") +
    scale_y_continuous(limits = c(0, 0.20)) +
    geom_hline(yintercept = prop_all[,1],
               linetype = "dashed", show.legend = TRUE, size = 0.6) +
    annotate('ribbon', x = c(-Inf, Inf), ymin = prop_all$lower, ymax = prop_all$upper, 
             alpha = 0.15, fill = 'forestgreen') +
    geom_text(data = plot_dat[plot_dat$Sig,], aes_string(x = "Chr", y = "upper" , label = "count" ), 
              size = 7.5, vjust = -0.7) +
    theme_classic() +
    theme(plot.title.position = c("panel"),
          legend.title = element_blank(),
          legend.position = c("bottom"),
          #legend.justification = c("right", "bottom"),
          #legend.box.just = "left",
          #legend.box.background = element_rect(),
          legend.box.background = element_rect(),
          #legend.box.margin = margin(4, 6, 6, 6),
          legend.text = element_text(size = 20, color = "black"),
          axis.text.x = element_text(size=20, margin = margin(5,0,0,0), color = "black"),
          axis.text.y = element_text(size=20, margin = margin(0,5,0,0), color = "black"),
          axis.title.x = element_text(size=30, margin = margin(10,0,0,0), color = "black"),
          axis.title.y = element_text(size=30, margin = margin(0,10,0,0), color = "black"),
          plot.title = element_text(size=40, margin = margin(0,0,0,0), color = "black"),
          plot.margin = margin(6,6,6,6))
  return(plot_vec)
}

# propChr(jseq.All.geno_Chr)
# plotChrprop(jseq.All.geno_Chr, "Chr") + coord_cartesian(ylim=c(-0.02, 0.05))


propChr(SSAV.geno_Chr[SSAV.geno_Chr$Chr != "Y",])
propChr(A.f.geno_Chr[!is.na(A.f.geno_Chr$Sig) &
                       !is.na(A.f.geno_Chr$Chr),])
propChr(A.m.geno_Chr[!is.na(A.m.geno_Chr$Sig) &
                        !is.na(A.m.geno_Chr$Chr) & A.m.geno_Chr$Chr != "Y",])



# Fisher's exact tests for diff. in proportions
######
fisher.test.results <- sampleTable
all.Chr <- NULL
Chr2.v.Chr3 <- NULL
Chr2.v.X <- NULL
Chr3.v.X <- NULL

for(i in 1:dim(sampleTable)[1]){
  dat <- get(paste0(sampleTable[i,2]))
  
  print(paste0(sampleTable[i,1]," SA candidates distribution across chromosomes:"))
  test <- fisher.test(dat$Sig, dat$Chr)
  print(test)
  all.Chr <- c(all.Chr, test$p.value)
  
  print(paste0(sampleTable[i,1]," SA candidates distribution Chr 2 vs Chr 3:"))
  test <- fisher.test(dat[dat$Chr == "2" | dat$Chr == "3",]$Sig, 
                    dat[dat$Chr == "2" | dat$Chr == "3",]$Chr)
  print(test)
  Chr2.v.Chr3 <- c(Chr2.v.Chr3, test$p.value)
  
  print(paste0(sampleTable[i,1]," SA candidates distribution Chr 2 vs X Chr:"))
  test <- fisher.test(dat[dat$Chr == "2" | dat$Chr == "X",]$Sig, 
                    dat[dat$Chr == "2" | dat$Chr == "X",]$Chr)
  print(test)
  Chr2.v.X <- c(Chr2.v.X, test$p.value)
  
  
  print(paste0(sampleTable[i,1]," SA candidates distribution Chr 3 vs X Chr:"))
  test <- fisher.test(dat[dat$Chr == "3" | dat$Chr == "X",]$Sig, 
                    dat[dat$Chr == "3" | dat$Chr == "X",]$Chr)
  Chr3.v.X <- c(Chr3.v.X, test$p.value)

  Sys.sleep(1.5)
}

fisher.test.results$all <- all.Chr
fisher.test.results$Chr2.v.Chr3 <- Chr2.v.Chr3
fisher.test.results$Chr2.v.X <- Chr2.v.X
fisher.test.results$Chr3.v.X <- Chr3.v.X
rm(all.Chr, Chr2.v.X, Chr2.v.Chr3, Chr3.v.X)

fisher.test.results[,-2]
write.table(fisher.test.results[,-2], file = "Results/Chr.dist.fisher.results.txt", quote = F, sep = "\t")
######



All_geno <- plotChrprop(SSAV.geno_Chr[SSAV.geno_Chr$Chr != "Y",], "Chromosome")

A.f_Sig <- plotChrprop(A.f.geno_Chr[!is.na(A.f.geno_Chr$Sig) &
                                       !is.na(A.f.geno_Chr$Chr),], "Chromosome")
A.m_Sig <- plotChrprop(A.m.geno_Chr[!is.na(A.m.geno_Chr$Sig) &
                           !is.na(A.m.geno_Chr$Chr) & A.m.geno_Chr$Chr != "Y",], 
            "Chromosome") 




pdf(file = "~/Desktop/UofT/SSAV_RNA/Plots/Chromosomal_dist_ALL.pdf",   # The directory you want to save the file in
    width = 10, # 12, 24 The width of the plot in inches
    height = 8) # 10, 20 The height of the plot in inches
# ggarrange(A.m_Sig, NA, A.f_Sig,
#           labels = c("A)", NA, "B)"),
#           widths = c(1, 0.05, 1),
#           ncol = 3,
#           font.label = list(size = 30))
All_geno
dev.off()





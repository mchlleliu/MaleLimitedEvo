# ASE Splicing analysis

library(DESeq2)
library(tidyverse)
library(GenomicFeatures)

# load data from JunctionSeq analysis
jseq.ASE = read.table("JunctionSeq/SDIU_ase/JSresults/SDIU_ASEallGenes.results.txt",
                      sep = "\t", header = TRUE)


# Set global variables
FDRThreshold = 0.1
mappedReadsThreshold = 50


# remove genes that were not assayed
jseq.ASE = jseq.ASE[!is.na(jseq.ASE$expr_F) & 
                      !is.na(jseq.ASE$expr_M),]

# remove exons with less than 50 reads mapping 
jseq.ASE = jseq.ASE[jseq.ASE$expr_F > mappedReadsThreshold & 
                      jseq.ASE$expr_M > mappedReadsThreshold,]
colnames(jseq.ASE)[2]="FlyBaseID"
str(jseq.ASE)

# asign genes with significant SSS
jseq.ASE <- jseq.ASE %>% mutate(SSS = ifelse(geneWisePadj <= FDRThreshold, TRUE, FALSE))


# list of genes with significant and non-significant SSS
ASE.sig.SSS <- unique(jseq.ASE$FlyBaseID[jseq.ASE$geneWisePadj < FDRThreshold])
ASE.nonsig.SSS <-  unique(jseq.ASE$FlyBaseID[jseq.ASE$geneWisePadj > FDRThreshold])
length(ASE.sig.SSS)
length(ASE.nonsig.SSS)



# get genes with low expression
decoder.ASE <- read.table("JunctionSeq/SDIU_ase/QoRTs.decoder.file.for.JunctionSeq.txt", header=T, stringsAsFactors = F)
colnames(decoder.ASE)[1]="unique.ID"
decoder.ASE$rep = rep(seq(1:3),4)

# split decoder into males and females
M.decoder.ASE <- decoder.ASE[decoder.ASE$sex=="M",]
F.decoder.ASE <- decoder.ASE[decoder.ASE$sex=="F",]

countFiles.ASE <- paste0("JunctionSeq/SDIU_ase/count.files/",decoder.ASE$unique.ID,"/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")
countFiles.m.ASE <- paste0("JunctionSeq/SDIU_ase/count.files/",M.decoder.ASE$unique.ID,"/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")
countFiles.f.ASE <- paste0("JunctionSeq/SDIU_ase/count.files/",F.decoder.ASE$unique.ID,"/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")


# use this function to read and combine all sample counts and create a counts matrix
# data.files is a list containing the paths to sample count files
makeCountMatrix <- function(data.files, decoder.file){
  count.file <- NULL
  for(i in data.files){
    tmp <- read.delim(i, header=F, sep = "\t")
    tmp <- tmp %>% separate(V1, into = c("FlyBaseID", "countBin"), sep = ":") %>%
      group_by(FlyBaseID) %>%
      summarise(count = sum(V2))
    if(!is.null(count.file)){
      count.file <- merge(count.file, tmp, by = "FlyBaseID", all = T)}
    else
      count.file <- tmp
  }

  colnames(count.file) <- c("FlyBaseID", decoder.file$unique.ID)

  return(count.file)
}


ASE.count.matrix <- count_file_reader_DESeqMatrix(countFiles.ASE, decoder.ASE)
ASE.count.matrix <- ASE.count.matrix %>% filter(!str_detect(FlyBaseID, "\\+"))
ASE.count.matrix <- ASE.count.matrix %>% remove_rownames %>% column_to_rownames(var="FlyBaseID")
ASE.colData <- decoder.ASE %>% remove_rownames %>% column_to_rownames(var="unique.ID")
dds.ASE <- DESeqDataSetFromMatrix(countData = ASE.count.matrix, 
                                  colData = ASE.colData, 
                                  design = ~ rep + sex)


load("JunctionSeq/SDIU_ase/Drosophila_melanogaster.BDGP6.28.102.exonLengths.RData", verbose = T)

exonic <- GRangesList(exonic)
exonic <- exonic[names(exonic) %in% rownames(ASE.count.matrix)]

rowRanges(dds.ASE) <- exonic

FPKM.ASE <- fpkm(dds.ASE)

filter.low.exp.genes.FEM <- data.frame(FPKM.ASE) %>% select(., contains("F_")) %>%
  mutate(totalCounts = rowSums(select_if(., is.numeric)))
filter.low.exp.genes.MALE <- data.frame(FPKM.ASE) %>% select(., contains("M_")) %>%
  mutate(totalCounts = rowSums(select_if(., is.numeric)))

filter.low.exp.genes.FEM.list <- rownames(filter.low.exp.genes.FEM[!is.na(filter.low.exp.genes.FEM$totalCounts) &
                                                                     filter.low.exp.genes.FEM$totalCounts  > quantile(filter.low.exp.genes.FEM$totalCounts, 0.25),])
filter.low.exp.genes.MALE.list <- rownames(filter.low.exp.genes.MALE[!is.na(filter.low.exp.genes.MALE$totalCounts) &
                                                                       filter.low.exp.genes.MALE$totalCounts  > quantile(filter.low.exp.genes.MALE$totalCounts, 0.25),])

filter.low.exp.genes.list <- unique(c(filter.low.exp.genes.MALE.list, filter.low.exp.genes.FEM.list))
filter.low.exp.genes.list <- filter.low.exp.genes.list[!(filter.low.exp.genes.list %in% Chrs[Chrs$Chr == "Y",]$FlyBaseID)]


filter.low.exp.genes.FEM.list.10perc <- rownames(filter.low.exp.genes.FEM[!is.na(filter.low.exp.genes.FEM$totalCounts) &
                                                                     filter.low.exp.genes.FEM$totalCounts  > quantile(filter.low.exp.genes.FEM$totalCounts, 0.10),])
filter.low.exp.genes.MALE.list.10perc <- rownames(filter.low.exp.genes.MALE[!is.na(filter.low.exp.genes.MALE$totalCounts) &
                                                                       filter.low.exp.genes.MALE$totalCounts  > quantile(filter.low.exp.genes.MALE$totalCounts, 0.10),])

filter.low.exp.genes.list.10perc <- unique(c(filter.low.exp.genes.MALE.list.10perc, filter.low.exp.genes.FEM.list.10perc))
filter.low.exp.genes.list.10perc <- filter.low.exp.genes.list.10perc[!(filter.low.exp.genes.list.10perc %in% Chrs[Chrs$Chr == "Y",]$FlyBaseID)]


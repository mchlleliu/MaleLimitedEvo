###################################
#
#                             Grieshop et al. 2023
#                        Original Author: Amardeep Singh
#                       Michelle's notes + editted pipeline
#             DsRed experimental evolution - transcriptomics analysis
#                 Running JunctionSeq to quantify isoform usage
# 
# 
###################################

require(QoRTs)

# Read QoRTs outputs into R to generate a sizeFactor by which to normalize read counts
# Tell QoRTs which directory to look into for directory that contain outputs. Each directory should be named to correspond to a sample (which should be decoded by the decoder file) and give it the location of your decoder file
qorts.results <- read.qc.results.data("/plas1/michelle.liu/SSAV/QoRTsCheck/", 
                                      decoder.files="/plas1/michelle.liu/Dmel_BDGP6.28/JunctionSeq.files/QoRTs.decoder.file.txt", 
                                      calc.DESeq2=TRUE)

# Save size factors for each sample in a text file
get.size.factors(qorts.results, outfile ="/plas1/michelle.liu/Dmel_BDGP6.28/JunctionSeq.files/RAL.size.Factors.GEO.txt")


### go back to counting through the command line
# ........


# after the Linux pipeline is done... run JunctionSeq using R
rm(list=ls())
require(DESeq2)
require(JunctionSeq)
require(BiocParallel) # This is a package used for paralellization of jobs

# Set global variables
numCores = 10
FDRThreshold = 0.01
mappedReadsThreshold = 50 # mean normalized read-pair counts?


# Load in decoder files and add fields for conditions

# Run this only once to edit the decoder file
# decoder.file =  read.delim("/plas1/michelle.liu/Dmel_BDGP6.28/JunctionSeq.files/QoRTs.decoder.file.txt", header = TRUE, stringsAsFactors=FALSE)
# sex.tmp = gsub("[A,C]..","", decoder.file$unique.ID)
# decoder.file$sex =  gsub("\\_.*","", sex.tmp)
# decoder.file$geno = gsub("[A,C]....","", decoder.file$unique.ID)
# decoder.file$pop = gsub("[[:digit:]]\\_.*","", decoder.file$unique.ID)
# write.table(decoder.file, "/plas1/michelle.liu/Dmel_BDGP6.28/JunctionSeq.files/QoRTs.decoder.file.for.JunctionSeq.txt", quote = F, row.names = FALSE, col.names = TRUE, sep = "\t")

# Load flat GFF generated earlier
flat_GFF <- "/plas1/michelle.liu/Dmel_BDGP6.28/JunctionSeq.files/count.files/withNovel.forJunctionSeq.gff.gz"

# Loading in decoder file
decoder <- read.table("/plas1/michelle.liu/Dmel_BDGP6.28/JunctionSeq.files/QoRTs.decoder.file.for.JunctionSeq.txt", 
                      header = TRUE, stringsAsFactors = FALSE)


# separate sample groups
A.f.decoder <- decoder[decoder$sex == "F",]
A.m.decoder <- decoder[decoder$sex == "M" & decoder$pop == "A",]
C.m.decoder <- decoder[decoder$sex == "M" & decoder$pop == "C",]

# set up pathways to the count files
#######
countFiles.A.f <- paste0("/plas1/michelle.liu/Dmel_BDGP6.28/JunctionSeq.files/count.files/", 
                         A.f.decoder$unique.ID, "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")

countFiles.A.m <- paste0("/plas1/michelle.liu/Dmel_BDGP6.28/JunctionSeq.files/count.files/", 
                         A.m.decoder$unique.ID, "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")

countFiles.C.m <- paste0("/plas1/michelle.liu/Dmel_BDGP6.28/JunctionSeq.files/count.files/", 
                         C.m.decoder$unique.ID, "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")
#######

###  Run the differential exon usage (DEU) analysis:
# Create a design dataframe which has a column for the condition you are testing
design.df.A.f <- data.frame(condition = factor(A.f.decoder$geno))
design.df.A.m <- data.frame(condition = factor(A.m.decoder$geno))
design.df.C.m <- data.frame(condition = factor(C.m.decoder$geno))

# change accordingly to which comparison you want to do.
focal.countFiles = countFiles.A.f
focal.design = design.df.A.f
focal.decoder = A.f.decoder

# Building the count set object that JunctionSeq will analyze and add to it all of the parameters of the analysis
count.set.object <- readJunctionSeqCounts(countfiles = focal.countFiles, 
                                          samplenames = focal.decoder$unique.ID, 
                                          design = focal.design, 
                                          flat.gff.file = flat_GFF, 
                                          nCores = numCores, verbose = TRUE)

# Generate size factors for normalization and load them into the count.set.object								
count.set.object <- estimateJunctionSeqSizeFactors(count.set.object)			

# Generate test specific dispersion estimates and load into count.set.object					  
count.set.object <- estimateJunctionSeqDispersions(count.set.object, nCores = numCores)	

# Fit the observed dispersions to a regression to create a fitted dispersion
count.set.object <- fitJunctionSeqDispersionFunction(count.set.object)						  

# Perform the hypothesis tests to test for differential splice junction/exon usage (DEU)
count.set.object <- testForDiffUsage(count.set.object, nCores = numCores)

# Calculate effect sizes and parameter estimates
count.set.object <- estimateEffectSizes(count.set.object, nCores = numCores)

# Save output to file
writeCompleteResults(count.set.object, 
                     outfile.prefix="/plas1/michelle.liu/SSAV/JunctionSeq/C.m.geno/C.m.geno.splicing", 
                     gzip.output = TRUE, 
                     FDR.threshold = FDRThreshold, 
                     save.allGenes = TRUE, 
                     save.sigGenes = TRUE, 
                     save.fit = FALSE, 
                     save.VST = FALSE, 
                     save.bedTracks = TRUE, bedtrack.format = c("BED", "GTF", "GFF3"), 
                     save.jscs = TRUE, verbose = TRUE)






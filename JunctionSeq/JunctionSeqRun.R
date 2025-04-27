###################################
#
#                             Grieshop et al. 2025
#                        Original Author: Amardeep Singh
#                       (Michelle's notes + editted pipeline)
#             DsRed experimental evolution - transcriptomics analysis
#           Running JunctionSeq to test for differential isoform usage
# 
###################################

# make sure to run "QoRTsRun.sh" and "GetSizeFactors.R" to generate the count files needed

rm(list=ls())
require(dplyr)
require(QoRTs)
require(DESeq2)
require(JunctionSeq)
require(stringr)
require(BiocParallel) # This is a package used for paralellization of jobs

# Set global variables
numCores = 10
FDRThreshold = 0.01
mappedReadsThreshold = 50 

# Load flat GFF generated earlier
flat_GFF <- "/plas1/michelle.liu/SSAV/QoRTs/sampleINFO.size.Factors.txt"

# Loading in sample info with factors (i.e., decoder file)
decoder <- read.table("/plas1/michelle.liu/SSAV/sample_INFO.tsv", header = FALSE)
decoder <- decoder[,-2] # don't need the second column (fastq file prefix)
colnames(decoder) <- c("unique.ID", "trt", "rep", "sex", "geno")

# separate sample groups
A.f.decoder <- decoder[str_detect(decoder[,1], "F"),]
A.m.decoder <- decoder[str_detect(decoder[,1], "M") & str_detect(decoder[,1], "A"),]
C.m.decoder <- decoder[str_detect(decoder[,1], "M") & str_detect(decoder[,1], "C"),]

# paths to the count files
#######
qorts_dir = "/plas1/michelle.liu/SSAV/QoRTs/"
countFiles.A.f <- paste0(qorts_dir, A.f.decoder$unique.ID, 
                         "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")
countFiles.A.m <- paste0(qorts_dir, A.m.decoder$unique.ID, 
                         "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")
countFiles.C.m <- paste0(qorts_dir, C.m.decoder$unique.ID,
                         "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")
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

# Save output to file (change prefix to fit which focal contrast)
writeCompleteResults(count.set.object, 
                     outfile.prefix="/plas1/michelle.liu/SSAV/JunctionSeq/A.f.geno", 
                     gzip.output = TRUE, 
                     FDR.threshold = FDRThreshold, 
                     save.allGenes = TRUE, 
                     save.sigGenes = TRUE, 
                     save.fit = FALSE, 
                     save.VST = FALSE, 
                     save.bedTracks = TRUE, bedtrack.format = c("BED", "GTF", "GFF3"), 
                     save.jscs = TRUE, verbose = TRUE)






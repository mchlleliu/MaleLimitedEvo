# JunctionSeq for ASE data



#!/bin/bash
# Script to run QoRTs
# run STAR aligner and sort as necessary before to generate BAM files.

# flags:
# -f: path to file (from home) containing sample information. The 2nd column of this file should contain the file prefix
# -b: directory containing Aligned BAM files
# -g: path to reference GFF/GTF file
# -o: output directory


# default path to sorted BAM files
BAM_PATH="."
OUT_DIR="."

# while do loop to read flags and set variables
while getopts "f:b:g:o:" flag
do
  # Look through the possible cases where $OPTARG is a variable made by getopts
  case "${flag}" in
    f) SAMPLE_INFO=$OPTARG
       echo "Sample file is: $SAMPLE_INFO"
       ;;
    b) BAM_DIR=$OPTARG
       echo "BAM directory is: $BAM_DIR" >&2
       ;;
    g) PATH_TO_GFF=$OPTARG
	echo "Reference GFF/GTF file is: $PATH_TO_GFF" >&2
	;;
	o) OUT_DIR=$OPTARG
	echo "Outputting to: $OUT_DIR" >&2
	;;
    # Exit the program if flag is unknown
   \?) echo "Invalid option: -$OPTARG" >&2
       exit 1
       ;;
  esac
done



# As Amardeep has noted, QoRTs seem to be memory hungry
# You may get some warnings like "NOTE: Unmatched Read Buffer Size > 100000 [Mem usage:[808MB / 3500MB]]"
# this is generally not a problem but to deal with this, you can either:
# increase allocated memory (modify the -XmxN tag), or
# sort your bam files by name (samtools sort -n -O bam -o file.sort.bam file.bam) and use the --nameSorted QoRTs parameter

cut -f 1 $SAMPLE_INFO | parallel -j 5 \
	"mkdir $OUT_DIR/{} && java -Xmx10G -jar /plas1/michelle.liu/bin/QoRTs-STABLE.jar QC \
		--stranded --maxReadLength 101 \
			$BAM_DIR/{}.name.sorted.bam $PATH_TO_GFF $OUT_DIR/{}/"


echo "Done running QoRTs :)"



nohup ./QoRTs_parallel.sh -f /plas1/michelle.liu/ASE_data/whole_bodies_sample.txt -b /plas1/michelle.liu/ASE_data/samtools.sortedBAM -g /plas1/michelle.liu/Dmel_BDGP6.28/Drosophila_melanogaster.BDGP6.28.102.gtf -o /plas1/michelle.liu/ASE_data/QoRTsCheck


cut -f 1 whole_bodies_sample.txt > QoRTs.decoder.ASE.txt
# a header is needed for this
cat "sample.ID"


# ---- R Code ---
require(QoRTs)

# Read QoRTs outputs into R to generate a sizeFactor by which to normalize read counts
# Tell QoRTs which directory to look into for directory that contain outputs. Each directory should be named to correspond to a sample (which should be decoded by the decoder file) and give it the location of your decoder file
qorts.results <- read.qc.results.data("/plas1/michelle.liu/ASE_data/QoRTsCheck/", 
decoder.files="/plas1/michelle.liu/ASE_data/QoRTsCheck/QoRTs.decoder.ASE.txt", 
calc.DESeq2=TRUE);

# Save size factors for each sample in a text file
get.size.factors(qorts.results, outfile ="/plas1/michelle.liu/ASE_data/JunctionSeq.files/RAL.size.Factors.GEO.txt");


# ---- Linux ----
# Create novel junction splice sites
java -Xmx20G -jar /plas1/michelle.liu/bin/QoRTs-STABLE.jar mergeNovelSplices --minCount 20 --stranded \
/plas1/michelle.liu/ASE_data/QoRTsCheck \
/plas1/michelle.liu/ASE_data/JunctionSeq.files/RAL.size.Factors.GEO.txt \
/plas1/michelle.liu/Dmel_BDGP6.28/Drosophila_melanogaster.BDGP6.28.102.gtf \
/plas1/michelle.liu/ASE_data/JunctionSeq.files/count.files

# ---- R Code ----

### JunctionSeq Script for Body Tissue

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
decoder.file =  read.delim("/plas1/michelle.liu/ASE_data/JunctionSeq.files/QoRTs.decoder.file.txt", header = TRUE, stringsAsFactors=FALSE)
sex.tmp=gsub("\\_.*","",decoder.file$sample.ID)
decoder.file$sex =  sex.tmp

write.table(decoder.file, "/plas1/michelle.liu/ASE_data/JunctionSeq.files/QoRTs.decoder.file.for.JunctionSeq.txt", quote=F, row.names=F, col.names=T, sep="\t")


# Loading in decoder file
decoder <- read.table("/plas1/michelle.liu/ASE_data/JunctionSeq.files/QoRTs.decoder.file.for.JunctionSeq.txt", header = TRUE, stringsAsFactors = FALSE)


countFiles <- paste0("/plas1/michelle.liu/ASE_data/JunctionSeq.files/count.files/",decoder$sample.ID,"/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")



###  Run the differential exon usage (DEU) analysis:
design.df <- data.frame(condition = factor(decoder$sex))



# Building the count set object that JunctionSeq will analyze and add to it all of the parameters of the analysis
count.set.object <- readJunctionSeqCounts(countfiles = countFiles, samplenames = decoder$sample.ID, design = design.df, flat.gff.file = "/plas1/michelle.liu/ASE_data/JunctionSeq.files/count.files/withNovel.forJunctionSeq.gff.gz", nCores = numCores, verbose = TRUE)

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
writeCompleteResults(count.set.object, outfile.prefix="/plas1/michelle.liu/ASE_data/JunctionSeq.files/SDIU_ASE", gzip.output = TRUE, FDR.threshold = FDRThreshold, save.allGenes = TRUE, save.sigGenes = TRUE, save.fit = FALSE, save.VST = FALSE, save.bedTracks = TRUE, bedtrack.format = c("BED", "GTF", "GFF3"), save.jscs = TRUE, verbose = TRUE)



# get gene lengths from GTF file
# -- on R
library(GenomicFeatures)
txdb <- makeTxDbFromGFF(gtfFile, format="gtf")
exonic <- exonsBy(txdb, by="gene")
save(exonic, file="/plas1/michelle.liu/Dmel_BDGP6.28/Drosophila_melanogaster.BDGP6.28.102.exonLengths.RData")
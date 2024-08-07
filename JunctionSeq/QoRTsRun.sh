#!/bin/bash
# Script to run QoRTs



# flags:
# -f: path to file (from home) containing sample information. The 2nd column of this file should contain the file prefix
# -b: directory containing Aligned BAM files
# -g: path to reference GFF file
# -o: output directory
# -h: help

# default path to sorted BAM files
BAM_PATH="."
OUT_DIR="."

# while do loop to read flags and set variables
while getopts "f:b:g:o:h" flag
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
	echo "Reference GFF file is: $PATH_TO_GFF" >&2
	;;
	o) OUT_DIR=$OPTARG
	echo "Outputting to: $OUT_DIR" >&2
	;;
    h) echo "flags:"
	echo "-f: path to file (from home) containing sample information. The 2nd column of this file should contain the file"
	echo "-b: directory containing Aligned BAM files"
	echo "-g: path to reference GFF file"
	echo "-o: output directory"
	echo "-h: help"
	exit 1
	;;
    # Exit the program if flag is unknown
   \?) echo "Invalid option: -$OPTARG" >&2
       exit 1
       ;;
  esac
done


# --------- make a flat gff -----------
# this is used to tell JunctionSeq to only test annotated splice junctions
# allocate 5G for java.
# --stranded option should be the same option used in the count-generation step
# input (path to reference gtf)
# output file name
# edit this to reflect the path to where you store QoRTs
java -Xmx5G -jar /plas1/michelle.liu/bin/QoRTs-STABLE.jar makeFlatGff --stranded $PATH_TO_GFF $OUT_DIR/JunctionSeq.flat.gff.gz



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

# ----------------------------------------------
# run this like so:
nohup ./QoRTs_run.sh -f /plas1/michelle.liu/SSAV/sample_INFO.tsv \
  -b /plas1/michelle.liu/SSAV/sortedBAM \
  -g /plas1/michelle.liu/Dmel_BDGP6.28/Drosophila_melanogaster.BDGP6.28.102.gtf \
  -o /plas1/michelle.liu/SSAV/QoRTsCheck




# For the next step, you would need to dip into R to get QoRTs to calculate the RAL size factors for the sample.
# Check the R script "JunctionSeqRun.R" for this step. below is a comment-out version

# ---- R Code ---
# require(QoRTs)

## Read QoRTs outputs into R to generate a sizeFactor by which to normalize read counts
## Tell QoRTs which directory to look into for directory that contain outputs. Each directory should be named to correspond to a sample (which should be decoded by the decoder file) and give it the location of your decoder file
# qorts.results <- read.qc.results.data("/plas1/michelle.liu/SSAV/QoRTsCheck/", 
#            decoder.files="/plas1/michelle.liu/Dmel_BDGP6.28/JunctionSeq.files/QoRTs.decoder.file.txt", 
#            calc.DESeq2=TRUE);

## Save size factors for each sample in a text file
# get.size.factors(qorts.results, outfile ="/plas1/michelle.liu/Dmel_BDGP6.28/JunctionSeq.files/RAL.size.Factors.GEO.txt");


# Back to command line
# ---- Linux ----
# Create novel junction splice sites and count them as well.
# (This isn't very relevant in the final analysis that we did, because we did not include the novel sites that could not be compared between sample groups)
# input: - directory where the counts per sample QoRTs output are stored
#        - path to GTF or GFF
#        - where you want the count files (now containing novel exons) to be stored 
java -Xmx20G -jar /plas1/michelle.liu/bin/QoRTs-STABLE.jar mergeNovelSplices --minCount 20 --stranded \
    /plas1/michelle.liu/SSAV/QoRTsCheck \
    /plas1/michelle.liu/Dmel_BDGP6.28/JunctionSeq.files/RAL.size.Factors.GEO.txt \
    /plas1/michelle.liu/Dmel_BDGP6.28/Drosophila_melanogaster.BDGP6.28.102.gtf \
    /plas1/michelle.liu/Dmel_BDGP6.28/JunctionSeq.files/count.files



# the rest of the analysis with JunctionSeq is done in R. Check "JunctionSeqRun.R"


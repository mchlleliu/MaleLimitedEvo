#!/bin/bash
# Script to run QoRTs to generate count file of exons and splice junctions
# These files eventually are used for the JunctionSeq analysis done on R (see "JunctionSeqRun.R")
# Requires "GetSizeFactors.R"

# usage:
# ./QoRTs_run.sh -f /path/to/sample_INFO/file -i /path/to/input/BAM/directory -g /path/to/GFF -o /path/to/output/directory

# arguments:
# -f | --filenames: path to file (from home) containing sample information. The 2nd column of this file should contain the file prefix
# -i | --input_dir: directory containing Aligned BAM files
# -g | --gff path to reference GFF file
# -o | --output_dir: output directory
# -h | --help: help

# default path to sorted BAM files
IN_DIR="."
OUT_DIR="."

QoRTs="/plas1/michelle.liu/bin/QoRTs-STABLE.jar"

VALID_ARGS=$(getopt -o f:i:g:o:h --long filenames:,input_dir:,gff:,output_dir:,help -- "$@")
# if argument invalid, exit program
if [[ $? -ne 0 ]]; then
    exit 1;
fi

# while do loop to read flags and set variables
while [ : ]; do
  # Look through the possible cases where $OPTARG is a variable made by getopts
  case "$1" in
    -f | --filenames) 
        SAMPLE_INFO=$2
        echo "Sample info file is: $SAMPLE_INFO"
        shift 2
        ;;
    -i | --input_dir) 
        IN_DIR=$2
        echo "input BAM directory is: $IN_DIR" >&2
        shift 2
        ;;
    -g | --gff) 
        PATH_TO_GFF=$2
	      echo "Reference GFF file is: $PATH_TO_GFF" >&2
	      shift 2
	      ;;
	  -o | --output_dir) 
	      OUT_DIR=$2
	      echo "Outputting to: $OUT_DIR" >&2
	      shift 2
	      ;;
	   # Help flag prints example for how to run script then exits program
    -h | --help) 
    echo "arguments & flags:"
	  echo "-f | --filenames: path to file (from home) containing sample information. The 1st column of this file should contain the sample ID/file prefix"
	  echo "-i | --input_dir: directory containing input BAM files"
  	echo "-g | --gff: path to reference GFF file"
	  echo "-o | --output_dir: output directory"
  	echo "-h | --help: print this help message"
  	echo ""
  	echo "usage:"
    echo "./QoRTs_run.sh -f /path/to/sample_INFO/file -i /path/to/input/BAM/directory -g /path/to/GFF -o /path/to/output/directory"
	  exit 1
	      ;;
    # Exit the program if flag is unknown
    -- ) shift; break ;;
    * ) break ;;
  esac
done



# actual pipeline code starts here

# --------- make a flat gff -----------
echo "Generating flattened annotation..."
# Overlapping genes are "flattened" so that each exon and splice junction is assigned a unique identifier.
#   this is generated for the statistical test of differential isoform usage performed by JunctionSeq, but ...
#   NOTE below that QoRTs use the unflattned GFF/GTF to count number of the reads overlapping each feature.
# Options:
#   allocate 5G for java.
#   --stranded flag should be the same option used in the count-generation step
#   input (path to reference gtf)
#   output file name
java -jar -Xmx5G $QoRTs makeFlatGff --stranded $PATH_TO_GFF $OUT_DIR/JunctionSeq.flat.gff.gz

echo "Flat annotation generated."



# --------- count reads overlapping each feature -----------
echo "Quantifying exons and splice junctions using QoRTs..."

# QoRTs seem to be memory hungry
# You may get some warnings like "NOTE: Unmatched Read Buffer Size > 100000 [Mem usage:[808MB / 3500MB]]"
# this is generally not a problem but to deal with this, you can either:
# increase allocated memory (modify the -XmxN tag), or
# sort your bam files by name (samtools sort -n -O bam -o file.sort.bam file.bam) and use the --nameSorted QoRTs parameter

cut -f 2 $SAMPLE_INFO | parallel -j 5 \
	"mkdir $OUT_DIR/{} && java -Xmx10G -jar $QoRTs QC \
		--stranded --maxReadLength 101 \
			$IN_DIR/{}_query.bam $PATH_TO_GFF $OUT_DIR/{}/"


echo "Done running QoRTs."




# ------ Obtaining count files again, but this time with novel junctions ------
# (Note: This isn't very relevant in the final analysis that we did, because we did not include the novel sites that could not be compared between sample groups)

# For this next step, you would need to dip into R to get QoRTs to calculate the RAL size factors for the sample.
echo "Obtaining sample size factors for read count normalization..."

# here I am calling the R script that has been modified to work with this current script
Rscript ./GetSizeFactors.R -q $OUT_DIR/ -f $SAMPLE_INFO -o $OUT_DIR/

# below is the original comment-out version (without adjustments for shell-calling)
# ---- R Code ---
# require(QoRTs)

## Read QoRTs outputs into R to generate a sizeFactor by which to normalize read counts
## Tell QoRTs which directory to look into for directory that contain outputs. Each directory should be named to correspond to a sample (which should be decoded by the decoder file) and give it the location of your decoder file
# qorts.results <- read.qc.results.data("/plas1/michelle.liu/SSAV/QoRTsCheck/", 
#            decoder.files="/plas1/michelle.liu/Dmel_BDGP6.28/JunctionSeq.files/QoRTs.decoder.file.txt", 
#            calc.DESeq2=TRUE);

## Save size factors for each sample in a text file
# get.size.factors(qorts.results, outfile ="/plas1/michelle.liu/Dmel_BDGP6.28/JunctionSeq.files/RAL.size.Factors.GEO.txt");

echo "Sample size factors obtained."


echo "Quantify novel splice junctions"
# Back to command line
# ---- Linux ----
# Create novel junction splice sites and generate count files with novel junctions
# input: - directory where the counts per sample QoRTs output are stored
#        - path to GTF or GFF (again, note that this is NOT the flat gff)
#        - where you want the count files (now containing novel exons) to be stored 
java -Xmx10G -jar $QoRTs mergeNovelSplices --minCount 20 --stranded \
    $OUT_DIR/ \
    $OUT_DIR/sampleINFO.size.Factors.txt \
    $PATH_TO_GFF \
    $OUT_DIR/

echo "Done!"


# the rest of the analysis with JunctionSeq is done in R. Check "JunctionSeqRun.R"
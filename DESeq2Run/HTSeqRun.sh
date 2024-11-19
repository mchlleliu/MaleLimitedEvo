#!/bin/bash
# running HTSeq using processed BAM files to generate count data for DESeq2


# flags:
# -i: path to input DIRECTORY.
# -g: path to reference GTF file
# -o: path to output DIRECTORY

# default paths
INPUT_DIR="."
OUTPUT_DIR="."

# while do loop to read flags and set variables
while getopts "i:g:o:" flag
do
  case "${flag}" in
    i) INPUT_DIR=$OPTARG
    echo "Path to input BAM directory is: $INPUT_DIR"
    ;;
    g) GTF_FILE=$OPTARG
    echo "GTF file is: $GTF_FILE" >&2
    ;;
    o) OUTPUT_DIR=$OPTARG
    echo "Output file is: $OUTPUT_DIR" >&2
    ;;
    # Exit the program if flag is unknown
    \?) echo "Invalid option: -$OPTARG" >&2
    exit 1
    ;;
  esac
done

cd $INPUT_DIR

SAMPLES=$(ls *_query.bam)

parallel -j 10 htseq-count -s no \
  --nonunique none --format bam \
  --secondary-alignments ignore --supplementary-alignments ignore \
  {} $GTF_FILE ">" $OUTPUT_DIR/{.}.tsv ::: $SAMPLES


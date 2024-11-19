#!/bin/bash
# Script to process FastQ reads to BAM files for HTseq

# flags:
# -f: path to FILE (from home) containing sample information. The 2nd column of this file should contain the file prefix
# -i: [optional] path to DIRECTORY containing input fastq files (from home)
# -d: [optional] path to reference DIRECTORY (from home)
# -g: path to annotation FILE
# -o: path to output DIRECTORY


# Set default path to data as the current directory
INPUT_DIR="."
REF_DIR="."
OUTPUT_DIR="."

# while do loop to read flags and set variables
while getopts "f:i:d:g:o:" flag
do
  # Look through the possible cases where $OPTARG is a variable made by getopts
  case "${flag}" in
    f) SAMPLE_INFO=$OPTARG
       echo "Sample file is: $SAMPLE_INFO"
       ;;
    i) INPUT_DIR=$OPTARG
       echo "Path to input Fastq files is: $INPUT_DIR" >&2
       ;;
    d) REF_DIR=$OPTARG
	echo "Reference directory is: $REF_DIR" >&2
	;;
    g) GTF_FILE=$OPTARG
	echo "GTF file is: $GTF_FILE" >&2
	;;
	o) OUTPUT_DIR=$OPTARG
	echo "Output directory is: $OUTPUT_DIR" >&2
	;;
    # Exit the program if flag is unknown
   \?) echo "Invalid option: -$OPTARG" >&2
       exit 1
       ;;
  esac
done


# ----- activate conda ENV -------
conda init
conda activate

mkdir $OUTPUT_DIR/Processed_BAM
mkdir $OUTPUT_DIR/tmp
REF=$(ls $REF_DIR/*.toplevel.fa)
SAMPLES=$(cut -f 2 $SAMPLE_INFO)

# ------------ Generate Genome ---------------- #
# run STAR to index the reference genome (Comment out if not done before)
# STAR --runThreadN 12 --runMode genomeGenerate --genomeDir $REF_DIR \
#    --genomeFastaFiles $REF_DIR/*.toplevel.fa --genomeSAindexNbases 12 \
#    --sjdbGTFfile $GTF_FILE --sjdbOverhang 100


# ------------------ Align reads -------------------- #
# parallel -j 5 STAR --runThreadN 12 --genomeDir $REF_DIR \
# 	--sjdbGTFfile $GTF_FILE --sjdbOverhang 100 \
# 	--readFilesIn $INPUT_DIR/{}_R1_subset.fastq.gz $INPUT_DIR/{}_R2_subset.fastq.gz \
# 	--readFilesCommand zcat --outFileNamePrefix $OUTPUT_DIR/{}_ \
# 	--outSAMtype BAM Unsorted --outSAMunmapped Within --twopassMode Basic ::: $SAMPLES


# ---- Sort STAR outputs by queryname by before MergeBam ---- #
parallel -j 5 java -jar /plas1/michelle.liu/bin/picard.jar SortSam \
    	I=$OUTPUT_DIR/{}_Aligned.out.bam \
    	O=$OUTPUT_DIR/tmp/{}_sort.bam \
    	SO=queryname ::: $SAMPLES


## ------ generate unmapped BAM -- needed for MergeBam
parallel -j 5 java -jar /plas1/michelle.liu/bin/picard.jar FastqToSam \
	  F1=$INPUT_DIR/{}_R1_subset.fastq.gz \
	  F2=$INPUT_DIR/{}_R2_subset.fastq.gz \
	  O=$OUTPUT_DIR/tmp/{}_unmapped.bam \
	  SM={} \
	  LB={} \
	  PL=Illumina \
	  RG={} ::: $SAMPLES

# --------- MergeBam -----------
parallel -j 5	java -jar /plas1/michelle.liu/bin/picard.jar MergeBamAlignment \
	  ALIGNED=$OUTPUT_DIR/tmp/{}_sort.bam \
	  UNMAPPED=$OUTPUT_DIR/tmp/{}_unmapped.bam \
	  O=$OUTPUT_DIR/tmp/{}_merged.bam \
	  R=$REF ::: $SAMPLES


rm $OUTPUT_DIR/tmp/*_unmapped.bam $OUTPUT_DIR/tmp/*_sort.bam


	# ---------------- MarkDuplicates ------------------
	# Too CPU-heavy, so had to use --java-options '-XX:ParallelGCThreads=1'
	# That worked, and keeps each process at ~ 100 %CPU, but they occasionally spike up to
	# ...200 %CPU. Also, each uses about 2.3 %MEM.
parallel -j 5	\
  java -jar -XX:+PrintGCDetails -XX:+UseParallelGC /plas1/michelle.liu/bin/picard.jar MarkDuplicates \
	  I=$OUTPUT_DIR/tmp/{}_merged.bam \
	  O=$OUTPUT_DIR/tmp/{}_dpl.bam \
	  M=$OUTPUT_DIR/{}_dpl.metrics.txt ::: $SAMPLES

rm $OUTPUT_DIR/tmp/*_merged.bam


## Sort BAM files AGAIN after the MergBam and MarkDuplicates steps BEFORE running HTSeq
parallel -j 5	java -jar /plas1/michelle.liu/bin/picard.jar SortSam \
	  I=$OUTPUT_DIR/tmp/{}_dpl.bam \
	  O=$OUTPUT_DIR/Processed_BAM/{}_query.bam \
	  SO=queryname ::: $SAMPLES


rm -r $OUTPUT_DIR/tmp
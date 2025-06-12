#!/bin/bash
#SBATCH --array=0-35
#SBATCH --time=5-0
#SBATCH --mem=80G
#SBATCH --cpus-per-task=12
#SBATCH --job-name=MaleLimitedEvo_Pipeline
#SBATCH -o /gpfs/home/rpr23sxu/scratch/Data/MaleLimitedEvo/ReadCounts/Output_Messages/MaleLimitedEvo_Pipeline-%a.out
#SBATCH -e /gpfs/home/rpr23sxu/scratch/Data/MaleLimitedEvo/ReadCounts/Error_Messages/MaleLimitedEvo_Pipeline-%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user= #<your.email@address>


# ----------------------------------------------------------------------
# Conda environment setup
# To recreate the RNA-Seq_env conda environment used by this script:
#   conda env create -f RNA-Seq_env.yml
# or, using the explicit package list:
#   conda create --name RNA-Seq_env --file RNA-Seq_env.txt
# Make sure RNA-Seq_env.yml and RNA-Seq_env.txt are present in this directory.
# ----------------------------------------------------------------------


# Load required modules
module load python/anaconda/2024.06  # Load Anaconda module to enable Conda commands

# Load Java 20+ for GATK compatibility
module unload java || true
module load OpenJDK/jdk-20.0.2

# Check Java version
JAVA_VERSION=$(java -version 2>&1 | awk -F[\".] '/version/ {print $2}')
if [ "$JAVA_VERSION" -lt 17 ]; then
  echo "Error: Java 17+ is required for GATK 4.2.6.1. Current version: $(java -version 2>&1 | head -n 1)"
  exit 1
fi

# Initialize Conda
source /gpfs/software/ada/python/anaconda/2024.06/etc/profile.d/conda.sh

# Activate RNA-Seq environment
conda activate RNA-Seq_env

# Test STAR installation
if ! STAR --version &>/dev/null; then
  echo "Error: STAR is not installed or not accessible in RNA-Seq_env."
  exit 1
fi

# Test HTSeq installation
if ! python -m HTSeq.scripts.count --help &>/dev/null; then
  echo "Error: HTSeq is not installed or not accessible in RNA-Seq_env."
  exit 1
fi

# Add GATK to PATH explicitly for this script
export PATH=/gpfs/home/rpr23sxu/gatk-4.2.6.1:$PATH

# Define input/output directories
GENOME_DIR="/gpfs/home/rpr23sxu/scratch/References/STAR"
GENOME_FASTA="/gpfs/home/rpr23sxu/scratch/References/Drosophila_melanogaster.BDGP6.28.dna.toplevel.fa"
GTF_FILE="/gpfs/home/rpr23sxu/scratch/References/Drosophila_melanogaster.BDGP6.28.102.gtf"
FASTQ_DIR="/gpfs/home/rpr23sxu/scratch/Data/MaleLimitedEvo/FastQ/raw_data"
OUTPUT_DIR="/gpfs/home/rpr23sxu/scratch/Data/MaleLimitedEvo/ReadCounts"

# Ensure OUTPUT_DIR exists
mkdir -p "$OUTPUT_DIR"

# Extract sample names from FASTQ files
FASTQ_FILES=($(ls $FASTQ_DIR/*_R1.fastq.gz | xargs -n 1 basename | rev | cut -c 13- | rev | uniq))

# Validate SLURM array index
if [ -z "${FASTQ_FILES[$SLURM_ARRAY_TASK_ID]}" ]; then
  echo "Error: SLURM_ARRAY_TASK_ID ($SLURM_ARRAY_TASK_ID) is out of range."
  exit 1
fi

SAMPLE=$(echo "${FASTQ_FILES[$SLURM_ARRAY_TASK_ID]}" | xargs)  # trim whitespace

# Validate input FASTQ files
if [ ! -f "${FASTQ_DIR}/${SAMPLE}_R1.fastq.gz" ] || [ ! -f "${FASTQ_DIR}/${SAMPLE}_R2.fastq.gz" ]; then
  echo "Error: Missing FASTQ files for sample ${SAMPLE}."
  exit 1
fi

# Step 1: Map RNA-seq reads with 2-pass procedure
STAR --runThreadN 12 --genomeDir $GENOME_DIR --sjdbGTFfile $GTF_FILE --sjdbOverhang 100 \
  --readFilesIn ${FASTQ_DIR}/${SAMPLE}_R1.fastq.gz ${FASTQ_DIR}/${SAMPLE}_R2.fastq.gz \
  --readFilesCommand zcat --outSAMtype BAM Unsorted --outFileNamePrefix ${OUTPUT_DIR}/${SAMPLE} \
  --twopassMode Basic 2>> ${OUTPUT_DIR}/${SAMPLE}_STAR_error.log

# Step 2: Generate unmapped BAM and merge with aligned BAM

# Generate unmapped BAM
gatk FastqToSam -F1 ${FASTQ_DIR}/${SAMPLE}_R1.fastq.gz -F2 ${FASTQ_DIR}/${SAMPLE}_R2.fastq.gz \
  -O ${OUTPUT_DIR}/${SAMPLE}_unmapped.bam -SM ${SAMPLE} -LB ${SAMPLE} -PL Illumina -RG ${SAMPLE} \
  2>> ${OUTPUT_DIR}/${SAMPLE}_FastqToSam_error.log

# Sort unmapped BAM by queryname (NEW: adheres to original pipeline step 4)
gatk SortSam -I ${OUTPUT_DIR}/${SAMPLE}_unmapped.bam -O ${OUTPUT_DIR}/${SAMPLE}_unmapped_query.bam -SO queryname

# Merge BAM files 
gatk MergeBamAlignment \
  -ALIGNED ${OUTPUT_DIR}/${SAMPLE}Aligned.out.bam \
  -UNMAPPED ${OUTPUT_DIR}/${SAMPLE}_unmapped_query.bam \
  -O ${OUTPUT_DIR}/${SAMPLE}_merged.bam \
  -R $GENOME_FASTA 2>> ${OUTPUT_DIR}/${SAMPLE}_MergeBam_error.log

# Step 3: Mark duplicates and sort
gatk MarkDuplicates -I ${OUTPUT_DIR}/${SAMPLE}_merged.bam -O ${OUTPUT_DIR}/${SAMPLE}_dedup.bam \
  -M ${OUTPUT_DIR}/${SAMPLE}_dpl.txt 2>> ${OUTPUT_DIR}/${SAMPLE}_MarkDuplicates_error.log

gatk SortSam -I ${OUTPUT_DIR}/${SAMPLE}_dedup.bam -O ${OUTPUT_DIR}/${SAMPLE}_sorted.bam \
  -SO queryname 2>> ${OUTPUT_DIR}/${SAMPLE}_SortSam_error.log

# Step 4: Run HTSeq count
python -m HTSeq.scripts.count -s no --nonunique none --format bam --secondary-alignments ignore --supplementary-alignments ignore \
  ${OUTPUT_DIR}/${SAMPLE}_sorted.bam $GTF_FILE > ${OUTPUT_DIR}/${SAMPLE}.tsv 2>> ${OUTPUT_DIR}/${SAMPLE}_HTSeq_error.log

# Deactivate Conda environment
conda deactivate

# Optional: Cleanup intermediate BAM files (keep only the final sorted BAM and .tsv)
rm -f ${OUTPUT_DIR}/${SAMPLE}Aligned.out.bam
rm -f ${OUTPUT_DIR}/${SAMPLE}_unmapped.bam
rm -f ${OUTPUT_DIR}/${SAMPLE}_unmapped_query.bam
rm -f ${OUTPUT_DIR}/${SAMPLE}_merged.bam
rm -f ${OUTPUT_DIR}/${SAMPLE}_dedup.bam

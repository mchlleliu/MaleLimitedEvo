#!/bin/bash
#SBATCH --array=0-1
#SBATCH --time=5-0
#SBATCH --mem=80G
#SBATCH --cpus-per-task=12
#SBATCH --job-name=MaleLimitedEvo_Pipeline
#SBATCH -o /gpfs/home/rpr23sxu/scratch/Data/MaleLimitedEvo/ReadCounts/Output_Messages/MaleLimitedEvo_Pipeline-%a.out
#SBATCH -e /gpfs/home/rpr23sxu/scratch/Data/MaleLimitedEvo/ReadCounts/Error_Messages/MaleLimitedEvo_Pipeline-%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rpr23sxu@uea.ac.uk

set -e  # Exit immediately if a command exits with a non-zero status
trap 'echo "Error occurred at line $LINENO"; exit 1' ERR  # Log error line
set -x  # Enable command tracing for debugging

# This script requires the following Conda environments to be set up:
# 1. GATK environment with Java 17:
#    conda create -n gatk_env -c conda-forge openjdk=17
#    Download GATK manually and add it to the PATH in the script.
# 2. HTSeq environment:
#    conda create -n htseq_env -c bioconda htseq

# Load required modules
module load STAR/2.7.10a
module load python/anaconda/2024.06  # Load Anaconda module to enable Conda commands

# Initialize Conda
source /gpfs/software/ada/python/anaconda/2024.06/etc/profile.d/conda.sh

# Define input/output directories
GENOME_DIR="/gpfs/home/rpr23sxu/scratch/References/STAR"
GENOME_FASTA="/gpfs/home/rpr23sxu/scratch/References/Drosophila_melanogaster.BDGP6.28.dna.toplevel.fa"
GTF_FILE="/gpfs/home/rpr23sxu/scratch/References/Drosophila_melanogaster.BDGP6.28.102.gtf"
FASTQ_DIR="/gpfs/home/rpr23sxu/scratch/Data/MaleLimitedEvo/test"
OUTPUT_DIR="/gpfs/home/rpr23sxu/scratch/Data/MaleLimitedEvo/ReadCounts"

# Extract sample names from FASTQ files
FASTQ_FILES=($(ls $FASTQ_DIR/*_R1.fastq.gz | xargs -n 1 basename | rev | cut -c 13- | rev | uniq))

# Validate SLURM array index
if [ -z "${FASTQ_FILES[$SLURM_ARRAY_TASK_ID]}" ]; then
  echo "Error: SLURM_ARRAY_TASK_ID ($SLURM_ARRAY_TASK_ID) is out of range. Check the number of samples in $FASTQ_DIR."
  exit 1
fi

SAMPLE=${FASTQ_FILES[$SLURM_ARRAY_TASK_ID]}

# Validate input FASTQ files
if [ ! -f "${FASTQ_DIR}/${SAMPLE}_R1.fastq.gz" ] || [ ! -f "${FASTQ_DIR}/${SAMPLE}_R2.fastq.gz" ]; then
  echo "Error: Missing FASTQ files for sample ${SAMPLE} in $FASTQ_DIR."
  exit 1
fi

# Step 1: Map RNA-seq reads with 2-pass procedure
STAR --runThreadN 12 --genomeDir $GENOME_DIR --sjdbGTFfile $GTF_FILE --sjdbOverhang 100 \
  --readFilesIn ${FASTQ_DIR}/${SAMPLE}_R1.fastq.gz ${FASTQ_DIR}/${SAMPLE}_R2.fastq.gz \
  --readFilesCommand zcat --outSAMtype BAM Unsorted --outFileNamePrefix ${OUTPUT_DIR}/${SAMPLE} \
  --twopassMode Basic 2>> ${OUTPUT_DIR}/${SAMPLE}_STAR_error.log

# Check if STAR mapping was successful
if [ ! -f ${OUTPUT_DIR}/${SAMPLE}Aligned.out.bam ]; then
  echo "Error: STAR mapping failed for ${SAMPLE}. Check ${OUTPUT_DIR}/${SAMPLE}_STAR_error.log for details."
  exit 1
fi

# Step 2: Generate unmapped BAM and merge with aligned BAM
# Activate Conda environment for GATK
conda activate gatk_env

# Set JAVA_HOME to ensure the correct Java version is used
export JAVA_HOME=$(dirname $(dirname $(which java)))
export PATH=$JAVA_HOME/bin:$PATH

# Add GATK to PATH
export PATH=/gpfs/home/rpr23sxu/gatk-4.6.0.0:$PATH

gatk FastqToSam -F1 ${FASTQ_DIR}/${SAMPLE}_R1.fastq.gz -F2 ${FASTQ_DIR}/${SAMPLE}_R2.fastq.gz -O /dev/stdout -SM ${SAMPLE} -LB ${SAMPLE} -PL Illumina -RG ${SAMPLE} 2>> ${OUTPUT_DIR}/${SAMPLE}_FastqToSam_error.log | \
gatk MergeBamAlignment -ALIGNED ${OUTPUT_DIR}/${SAMPLE}Aligned.out.bam -UNMAPPED /dev/stdin -O /dev/stdout -R $GENOME_FASTA 2>> ${OUTPUT_DIR}/${SAMPLE}_MergeBam_error.log | \
gatk MarkDuplicates -I /dev/stdin -O /dev/stdout -M ${OUTPUT_DIR}/${SAMPLE}_dpl.txt 2>> ${OUTPUT_DIR}/${SAMPLE}_MarkDuplicates_error.log | \
gatk SortSam -I /dev/stdin -O ${OUTPUT_DIR}/${SAMPLE}_sorted.bam -SO queryname 2>> ${OUTPUT_DIR}/${SAMPLE}_SortSam_error.log

# Deactivate Conda environment for GATK
conda deactivate

# Check if final sorted BAM file is created
if [ ! -f ${OUTPUT_DIR}/${SAMPLE}_sorted.bam ]; then
  echo "Error: Final sorted BAM file not created for ${SAMPLE}. Check logs for details."
  exit 1
fi

# Step 3: Run HTSeq count
# Activate Conda environment for HTSeq
conda activate htseq_env

python -m HTSeq.scripts.count -s no --nonunique none --format bam --secondary-alignments ignore --supplementary-alignments ignore \
  ${OUTPUT_DIR}/${SAMPLE}_sorted.bam $GTF_FILE > ${OUTPUT_DIR}/${SAMPLE}.tsv 2>> ${OUTPUT_DIR}/${SAMPLE}_HTSeq_error.log

# Deactivate Conda environment for HTSeq
conda deactivate

# Check if HTSeq count was successful
if [ ! -f ${OUTPUT_DIR}/${SAMPLE}.tsv ]; then
  echo "Error: HTSeq count failed for ${SAMPLE}. Check ${OUTPUT_DIR}/${SAMPLE}_HTSeq_error.log for details."
  exit 1
fi
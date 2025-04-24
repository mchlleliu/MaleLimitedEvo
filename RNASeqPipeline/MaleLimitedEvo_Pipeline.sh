#!/bin/bash
#SBATCH --array=0-35
#SBATCH --time=5-0
#SBATCH --mem=80G
#SBATCH --cpus-per-task=12
#SBATCH --job-name=MaleLimitedEvo_Pipeline
#SBATCH -o /gpfs/data/GrieshopLab/Karl/Transcriptomics/Output_Messages/MaleLimitedEvo_Pipeline-%a.out
#SBATCH -e /gpfs/data/GrieshopLab/Karl/Transcriptomics/Error_Messages/MaleLimitedEvo_Pipeline-%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rpr23sxu@uea.ac.uk

# Load all required modules
module load STAR/2.7.9a
module load gatk/4.2.0.0
module load python/3.8.5
module load HTSeq/0.11.2

# Define input/output directories
GENOME_DIR="/gpfs/data/GrieshopLab/Karl/Transcriptomics/References/STAR"
GENOME_FASTA="/gpfs/data/GrieshopLab/Karl/Transcriptomics/References/Drosophila_melanogaster.BDGP6.28.dna.toplevel.fa"
GTF_FILE="/gpfs/data/GrieshopLab/Karl/Transcriptomics/References/Drosophila_melanogaster.BDGP6.28.102.gtf"
FASTQ_DIR="/gpfs/data/GrieshopLab/Karl/Transcriptomics/Data"
OUTPUT_DIR="/gpfs/data/GrieshopLab/Karl/Transcriptomics/Results"
FASTQ_FILES=($(ls $FASTQ_DIR/*.fastq.gz | rev | cut -c 13- | rev | uniq))
SAMPLE=${FASTQ_FILES[$SLURM_ARRAY_TASK_ID]}

# Step 1: Map RNA-seq reads with 2-pass procedure
STAR --runThreadN 12 --genomeDir $GENOME_DIR --sjdbGTFfile $GTF_FILE --sjdbOverhang 100 --readFilesIn 
${FASTQ_DIR}/${SAMPLE}_R1.fastq.gz ${FASTQ_DIR}/${SAMPLE}_R2.fastq.gz --readFilesCommand zcat --outSAMtype BAM 
Unsorted --outFileNamePrefix ${OUTPUT_DIR}/${SAMPLE} --twopassMode Basic

# Check if STAR mapping was successful
if [ ! -f ${OUTPUT_DIR}/${SAMPLE}Aligned.out.bam ]; then
  echo "Error: STAR mapping failed for ${SAMPLE}"
  exit 1
fi

# Step 2: Generate unmapped BAM and merge with aligned BAM
gatk FastqToSam -F1 ${FASTQ_DIR}/${SAMPLE}_R1.fastq.gz -F2 ${FASTQ_DIR}/${SAMPLE}_R2.fastq.gz -O /dev/stdout 
-SM ${SAMPLE} -LB ${SAMPLE} -PL Illumina -RG ${SAMPLE} | \
gatk MergeBamAlignment -ALIGNED ${OUTPUT_DIR}/${SAMPLE}Aligned.out.bam -UNMAPPED /dev/stdin -O /dev/stdout -R 
$GENOME_FASTA | \
gatk MarkDuplicates -I /dev/stdin -O /dev/stdout -M ${OUTPUT_DIR}/${SAMPLE}_dpl.txt | \
gatk SortSam -I /dev/stdin -O ${OUTPUT_DIR}/${SAMPLE}_sorted.bam -SO queryname

# Check if final sorted BAM file is created
if [ ! -f ${OUTPUT_DIR}/${SAMPLE}_sorted.bam ]; then
  echo "Error: Final sorted BAM file not created for ${SAMPLE}"
  exit 1
fi

# Step 3: Run HTSeq count
python -m HTSeq.scripts.count -s no --nonunique none --format bam --secondary-alignments ignore 
--supplementary-alignments ignore ${OUTPUT_DIR}/${SAMPLE}_sorted.bam $GTF_FILE > ${OUTPUT_DIR}/${SAMPLE}.tsv

# Check if HTSeq count was successful
if [ ! -f ${OUTPUT_DIR}/${SAMPLE}.tsv ]; then
  echo "Error: HTSeq count failed for ${SAMPLE}"
  exit 1
fi

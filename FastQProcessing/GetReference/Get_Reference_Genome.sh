#!/bin/bash

#SBATCH --time=2-0
#SBATCH --mem=20G
#SBATCH --job-name=Get_Reference_Genome
#SBATCH -o /gpfs/home/rpr23sxu/scratch/References/Output_Messages/Get_Reference_Genome.out
#SBATCH -e /gpfs/home/rpr23sxu/scratch/References/Error_Messages/Get_Reference_Genome.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rpr23sxu@uea.ac.uk

# Required Output Directories:
# - /Output_Messages
# - /Error_Messages
# - /References

# Required Input Files:
# - None (the script downloads the necessary files)

# Load all required modules
module load samtools/1.16.1
module load bwa/0.7.17
module load STAR/2.7.10a

# Define output directory
OUTPUT_DIR="/gpfs/home/rpr23sxu/scratch/References/"
REF_FILE="$OUTPUT_DIR/Drosophila_melanogaster.BDGP6.28.dna.toplevel.fa"
REF_FILE_GZ="$REF_FILE.gz"
GTF_FILE="$OUTPUT_DIR/Drosophila_melanogaster.BDGP6.28.102.gtf"
GENOME_DIR="/gpfs/home/rpr23sxu/scratch/References/STAR"

# Check if the reference genome file already exists
if [ ! -f "$REF_FILE" ]; then
  # Download latest reference genome from Ensembl
  wget https://ftp.ensembl.org/pub/release-102/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.28.dna.toplevel.fa.gz -P "$OUTPUT_DIR"

  # Unzip the .fa.gz reference file
  gunzip "$REF_FILE_GZ"
else
  echo "Reference genome file already exists. Skipping download."
fi

# Check if the GTF file already exists
if [ ! -f "$GTF_FILE" ]; then
  # Download the GTF file
  wget https://ftp.ensembl.org/pub/release-102/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.28.102.gtf.gz -P "$OUTPUT_DIR"
  
  # Unzip the .gtf.gz file
  gunzip "$GTF_FILE.gz"
else
  echo "GTF file already exists. Skipping download."
fi

# Check if BWA index files already exist
if [ ! -f "${REF_FILE}.bwt" ]; then
  # Index the downloaded reference genome with bwa
  bwa index -a bwtsw "$REF_FILE"
else
  echo "BWA index files already exist. Skipping BWA indexing."
fi

# Check if samtools index file already exists
if [ ! -f "${REF_FILE}.fai" ]; then
  # Index the reference genome with samtools
  samtools faidx "$REF_FILE"
else
  echo "Samtools index file already exists. Skipping samtools indexing."
fi

# Check if the GATK sequence dictionary file already exists
if [ ! -f "${REF_FILE%.fa}.dict" ] && [ ! -f "${REF_FILE%.fa}" ]; then
  echo "GATK sequence dictionary file not found. Generating it..."
  gatk CreateSequenceDictionary -R "$REF_FILE"
else
  echo "GATK sequence dictionary file already exists. Skipping dictionary generation."
fi

# Check if STAR genome indices already exist
if [ ! -d "$GENOME_DIR" ]; then
  mkdir -p "$GENOME_DIR"
  STAR --runThreadN 12 --runMode genomeGenerate --genomeSAindexNbases 12 --genomeDir "$GENOME_DIR" --genomeFastaFiles "$REF_FILE" --sjdbGTFfile "$GTF_FILE"
else
  echo "STAR genome indices already exist. Skipping STAR indexing."
fi


#!/bin/bash

#SBATCH --time=1-0
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=Rename_Files
#SBATCH -o /gpfs/home/rpr23sxu/Teaching/MaleLimitedEvo/FastQProcessing/RenamingFiles/Output_Messages/Rename_Files.out
#SBATCH -e /gpfs/home/rpr23sxu/Teaching/MaleLimitedEvo/FastQProcessing/RenamingFiles/Error_Messages/Rename_Files.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rpr23sxu@uea.ac.uk

set -e
# Exit immediately if a command exits with a non-zero status

set -u
# Treat unset variables as an error and exit immediately

set -o pipefail
# Return the exit status of the last command in the pipeline that failed

# Define the directory containing the .tsv files
TSV_DIR="/gpfs/home/rpr23sxu/Teaching/MaleLimitedEvo/FastQProcessing/RenamingFiles"

cd $TSV_DIR

## Rename ReadCount folders (for easier sorting in R)

for i in $(ls *.tsv)
# Loop through each file in the current directory with a .tsv extension

do
	# Extract the sample name components
	component1=$(echo "$i" | cut -d '.' -f5 | cut -d '_' -f 1)
	# First component: Extract the fifth field when split by '.' and take the first field when split by '_'
	
	component2=$(echo "$i" | cut -d '.' -f5 | cut -d '_' -f 2)
	# Second component: Extract the fifth field when split by '.' and take the second field when split by '_'
	
	component3=$(echo "$i" | cut -d '.' -f5 | cut -d '_' -f 3)
	# Third component: Extract the fifth field when split by '.' and take the third field when split by '_'
	
	# Reorder the components as desired
	sampleName="$component2"_"$component3"_"$component1"
	
	fileName="$sampleName"".tsv"
	mv $i $fileName
	# Rename the file to the new filename
done
# End of the loop

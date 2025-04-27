# Transcriptomics File Renaming

This directory contains files and scripts related to the renaming of transcriptomics data files for easier sorting and analysis in R.

## Files

- `Get_Reference_Genome.sh`: A bash script used to rename `.tsv` files based on specific components extracted from their original filenames.
- `readcount_filename_list.txt`: A list of the original filenames before renaming.
- `newfilenames.txt`: A list of the new filenames after renaming.

## Usage

To rename the files, run the `Get_Reference_Genome.sh` script in the directory containing the `.tsv` files:
```bash
sbatch Get_Reference_Genome.sh
```

## Description of `Get_Reference_Genome.sh`

The script performs the following steps:

1. Loops through each `.tsv` file in the directory.
2. Extracts specific components from the filename.
3. Reorders the components to form a new filename.
4. Renames the file to the new filename.

## Example

Original filename: `NS.1701.004.NEBNext_dual_i7_100---NEBNext_dual_i5_100.A2_F_NR_merged_dpl_query.tsv`

New filename: `F_NR_A2.tsv`

## Notes

- Ensure that the filenames follow the expected format for the script to work correctly.
- The script uses `set -e`, `set -u`, and `set -o pipefail` to handle errors and unset variables.


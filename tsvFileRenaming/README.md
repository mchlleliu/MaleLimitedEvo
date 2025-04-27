# Transcriptomics File Renaming

This directory contains files and scripts related to the renaming of transcriptomics data files for easier sorting and analysis in R.

## Files

- `renameTSVfiles.sh`: A bash script used to rename `.tsv` files based on specific components extracted from their original filenames. **Note**: This script is an example and does not directly produce the correct `.tsv` filenames required for the R script in `MaleLimitedEvo/DESeq2Run/DESeq2Run.R`.
- `readcount_filename_list.txt`: A list of the original filenames before renaming.
- `newfilenames.txt`: A list of the new filenames after renaming.

## Usage

To rename the files, run the `renameTSVfiles.sh` script in the directory containing the `.tsv` files:
```bash
sbatch renameTSVfiles.sh
```

## Description of `renameTSVfiles.sh`

The script performs the following steps:
1. Loops through each `.tsv` file in the directory.
2. Extracts specific components from the filename.
3. Reorders the components to form a new filename.
4. Renames the file to the new filename.

**Important**: The script is intended as an example of how filenames can be renamed. The output filenames may need further adjustment to match the exact format required by downstream scripts, such as `MaleLimitedEvo/DESeq2Run/DESeq2Run.R`.

## Example

The following example demonstrates how the renaming process works. Note that this example is illustrative and does not reflect the actual `.tsv` files used in the R script located in `MaleLimitedEvo/DESeq2Run/DESeq2Run.R`.

Original filename:  
`NS.1701.004.NEBNext_dual_i7_100---NEBNext_dual_i5_100.A2_F_NR_merged_dpl_query.tsv`

New filename:  
`F_NR_A2.tsv`

## Notes

- Ensure that the filenames follow the expected format for the script to work correctly.
- The script uses `set -e`, `set -u`, and `set -o pipefail` to handle errors and unset variables.
- The actual `.tsv` files used in the R script may have undergone additional processing or renaming steps not covered in this example.


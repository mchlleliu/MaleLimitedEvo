# FastQ File Renaming

This directory contains files and scripts related to the renaming of FastQ files for the Male-limited Experimental Evolution program in the Agrawal lab.

## Files

- `rename_files.sh`: A bash script used to rename `.fastq.gz` files based on a predefined mapping of original filenames to new filenames.
- `SampleName_key.csv`: A key mapping the original file numbers (1-36) to their corresponding sample names.
- `modified_filename_list.txt`: A list of the renamed filenames after running the script.

## Usage

To rename the files, run the `rename_files.sh` script in the directory containing the original `.fastq.gz` files:
```bash
bash rename_files.sh
```

## Description of `rename_files.sh`

The script performs the following steps:

1. Loops through each `.fastq.gz` file in the directory.
2. Renames the files based on a predefined mapping of numbers (1-36) to sample names.
3. Ensures that paired-end suffixes (`_R1` and `_R2`) and file extensions (`.fastq.gz`) are preserved.

## Example

Original filename: `36_R.fastq.gz`  
New filename: `C6_M_NR_R.fastq.gz`

## Notes

- The script requires the `rename` program. Run `rename --help` to verify its installation and usage.
- The renaming commands must be executed in the correct order to avoid conflicts (e.g., renaming `6_R` before `36_R` would cause errors).
- The renaming process was validated by comparing the renamed filenames in `modified_filename_list.txt` against the expected names derived from `SampleName_key.csv`.

## SampleName_key.csv

This file maps the original numbers (1-36) to their corresponding sample names:

| Number | Name  |
|--------|-------|
| 1      | A1RM  |
| 2      | A1NM  |
| 3      | A1RF  |
| 4      | A1NF  |
| ...    | ...   |
| 36     | C6NM  |

For the full mapping, refer to the `SampleName_key.csv` file.

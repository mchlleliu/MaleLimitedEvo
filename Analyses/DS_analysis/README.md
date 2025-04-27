# Differential Splicing Analysis

This folder contains scripts and resources for performing **Differential Splicing (DS) Analysis** as part of the study, **"Expression divergence in response to sex-biased selection"**. Each script generates specific figures or tables included in the study.

## Overview

The analyses in this folder investigate:
- Masculinization or feminization of splicing profiles in Red vs. NonRed populations
- Overlap between differentially spliced and differentially expressed genes
- Comparison of splicing profiles across experimental and control groups

## File Descriptions

1. **[DS_TableS3.S4.R](https://github.com/mchlleliu/MaleLimitedEvo/blob/main/Analyses/DS_analysis/DS_TableS3.S4.R)**  
   Generates **Tables S3 and S4**, analyzing differentially spliced genes between Red and NonRed gene pools. 
This script also examines overlap with sex-specific splicing (SSS) genes based on prior datasets.

2. **[SplicingPhi_Fig5.6B_TableS8.R](https://github.com/mchlleliu/MaleLimitedEvo/blob/main/Analyses/DS_analysis/SplicingPhi_Fig5.6B_TableS8.R)**  
   Produces **Figure 5**, **Figure 6B**, and **Table S8**, investigating masculinization and feminization of splicing profiles. The script calculates metrics like Phi (Φ) for comparing splicing profiles between experimental and control groups.

3. 
**[dimorphic.subset.list.txt](https://github.com/mchlleliu/MaleLimitedEvo/blob/main/Analyses/DS_analysis/dimorphic.subset.list.txt)** 
   A list of genes with consistent sexual dimorphism in splicing profile identified using several external 
datasets and used for testing the masculinization/feminization of splicing profiles in Red vs. NonRed samples.

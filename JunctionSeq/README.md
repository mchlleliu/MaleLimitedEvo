# JunctionSeq Analysis

This subdirectory contains scripts and resources for analyzing RNA-Seq data using **JunctionSeq**, a tool 
designed for detecting differential exon and splice junction usage in RNA-Seq datasets.

## Overview

The analyses in this folder focus on identifying and quantifying differential splicing events across 
Red and NonRed gene pools. 

## Files

1. **[QoRTsRun.sh](https://github.com/mchlleliu/MaleLimitedEvo/blob/main/JunctionSeq/QoRTsRun.sh)**  
   Script for running QoRTs to generate exon and splice junction count data for differential isoform 
usage analysis using Junctionseq. This script takes the sorted BAM files generated (see 
[FastQProcessing](https://github.com/mchlleliu/MaleLimitedEvo/blob/main/FastQProcessing) as input.

2. **[GetSizeFactors.R](https://github.com/mchlleliu/MaleLimitedEvo/blob/main/JunctionSeq/GetSizeFactors.R)**
   Script for calculating size factors, which are normalization factors used to account for library size   
differences in RNA-Seq data. Used within QoRTsRun.sh. This file should be placed within the same directory as 
QoRTaRun.sh.

3. **[JunctionSeqRun.R](https://github.com/mchlleliu/MaleLimitedEvo/blob/main/JunctionSeq/JunctionSeqRun.R)**
   The main R script for running the JunctionSeq analysis.

For questions or issues, contact [Michelle Liu](mailto:mchelle.liu@mail.utoronto.ca)

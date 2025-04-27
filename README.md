# MaleLimitedEvo

This repository serves as the documentation for the analyses and scripts used in the study, Grieshop, K., Liu, M. J., Frost, R. S., Lindsay, M. P., Bayoumi, M., Brengdahl, M. I., Molnar, R. I. & Agrawal, A. F. (2025). **"Expression divergence in response to sex-biased selection"**. 
RNA-Seq reads associated with this publication can be found on the SRA under BioProject accession PRJNA1184789

## Background

In six replicate populations of 1000 flies, a dominant marker (DsRed) was used to force a “Red” pool of genetically variable Chromosome 2 copies through exclusive father-to-son inheritance, while a complimentary pool of “NonRed” chromosomes was inherited primarily from mothers to daughters. After 100 generations, we demonstrated the effect of Red male-limited chromosomes in increasing male mating success. We analysed differentially expressed genes with and without Red chromosomes. 


## Contents

### Directories
1. **[Analyses](https://github.com/mchlleliu/MaleLimitedEvo/tree/main/Analyses)**  
   Contains scripts for performing differential expression (DE) and differential splicing (DS) analyses, generating figures, and statistical data exploration.

2. **[DESeq2Run](https://github.com/mchlleliu/MaleLimitedEvo/tree/main/DESeq2Run)**  
   Resources and scripts for running DESeq2.

3. **[JunctionSeq](https://github.com/mchlleliu/MaleLimitedEvo/tree/main/JunctionSeq)**  
   Scripts and resources for splicing analysis using JunctionSeq.

4. **[FastQProcessing](https://github.com/mchlleliu/MaleLimitedEvo/tree/main/FastQProcessing)**  
  Scripts and resources for processing FastQ files, including fetching reference genomes, sorting, and counting 
for downstream RNA-Seq analyses.

5. **[Results](https://github.com/mchlleliu/MaleLimitedEvo/tree/main/Results)**  
   A directory for output files from DESeq2 and JunctionSeq, used in the analyses.

### Files
1. **[README.md](https://github.com/mchlleliu/MaleLimitedEvo/blob/main/README.md)**  
   This README file, providing an overview of the repository contents.

2. **[sample_INFO.tsv](https://github.com/mchlleliu/MaleLimitedEvo/blob/main/sample_INFO.tsv)**  
   Metadata file containing detailed information on RNA-Seq samples used in the analyses.

---



## Contacts

Karl Grieshop<br />
Email: [k.grieshop@uea.ac.uk](mailto:K.Grieshop@uea.ac.uk)<br />
GitHub: [karlgrieshop](https://github.com/karlgrieshop)

Michelle J. Liu  
Email: [mchelle.liu@mail.utoronto.ca](mailto:mchelle.liu@mail.utoronto.ca)  
GitHub: [mchlleliu](https://github.com/mchlleliu)

Aneil F. Agrawal<br />
Email: [a.agrawal@utoronto.ca](mailto:a.agrawal@utoronto.ca)

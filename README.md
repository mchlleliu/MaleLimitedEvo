# Expression divergence in response to sex-biased selection
Repository of Code and Data associated with Grieshop et al. 202X: (DOI)
RNA-Seq reads associated with this publication can be found on the SRA under BioProject accession PRJNA1184789

## Notes
### Abbreviations used:
**SSAV**: "Segregating Sexually Antagonistic Variants". The moniker given to the Experimental DsRed populations
**A.m/f**: Experimental males/females. (These populations were first initialized as treatment "A", hence the name)
**C.m**: Control males


### sample_INFO.tsv file
[sample_INFO.tsv](https://github.com/mchlleliu/SSAV_RNA/blob/main/sample_INFO.tsv) is a tab-separated file (with no header) that contains columns denoting sampleID, prefix, and condition. This file is used in bash and R scripts to generate exon count tables (See [JunctionSeq](https://github.com/mchlleliu/SSAV_RNA/tree/main/JunctionSeq)). 
Column names:
unique.ID[\t]fastq.Prefix[\t]trt[\t]rep[\t]sex[\t]geno
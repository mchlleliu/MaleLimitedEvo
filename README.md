# (Insert paper title here)
Repository of Code and Data associated with Grieshop et al. 202X: (DOI)



## Notes
### "SSAV"
"SSAV" stand for "Segregating Sexually Antagonistic Variants" and is the moniker given to the Experimental DsRed populations

### sample_INFO.tsv file
[sample_INFO.tsv](https://github.com/mchlleliu/SSAV_RNA/blob/main/sample_INFO.tsv) is a tab-separated file (with no header) that contains columns denoting sampleID, prefix, and condition. This file is used in bash and R scripts to generate exon count tables (See [JunctionSeq](https://github.com/mchlleliu/SSAV_RNA/tree/main/JunctionSeq)). 
Column names:
unique.ID[\t]fastq.Prefix[\t]trt[\t]rep[\t]sex[\t]geno
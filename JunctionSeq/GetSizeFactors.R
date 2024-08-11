#!/usr/bin/env Rscript

# Usage:
#     Rscript GetSizeFactors.R -q /QoRTs/output/directory -f /path/to/sample/info/file.txt -o /path/to/output/file.txt

options(show.error.locations = TRUE)

if(!require("optparse")) install.packages("optparse")
require(optparse)

# list of options:
option_list <- list(make_option(c("-q","--qorts_dir"), default = getwd(), metavar="character", type = "character", help="Directory of QoRTs count files"),
                    make_option(c("-f","--filenames"), default=NA, metavar="character", type = "character", help=".txt file containing sample IDs/file prefix"),
                    make_option(c("-o","--out_dir"), default = getwd(), metavar="character", type = "character", help="Output directory location")
)

# get arguments
opt <- parse_args(OptionParser(option_list=option_list), positional_arguments = TRUE)

decoder <- opt$options$filenames
decoder <- read.table(decoder, header=FALSE)
decoder <- decoder[,1]

# generate size factors to normalize each sample read counts
if(!require("QoRTs")) install.packages("QoRTs")
require(QoRTs)

# Read QoRTs outputs into R to generate a sizeFactor by which to normalize read counts
# Tell QoRTs which directory to look into for directory that contain outputs. 
# Each directory should be named to correspond to a sample (which should be decoded by the decoder file) and give it the location of your decoder file
qorts.results <- read.qc.results.data(opt$options$qorts_dir, 
                                      decoder = decoder, 
                                      calc.DESeq2 = TRUE)

# Save size factors for each sample in a text file
get.size.factors(qorts.results, outfile = paste0(opt$options$out_dir,"/RAL.size.Factors.GEO.txt"))


### go back to counting through the command line
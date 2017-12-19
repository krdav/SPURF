# These modules are necessary for sumrep.

print('Installing modules for sumrep')

install.packages(c("alakazam", "ape", "data.table", "dplyr", "HDMD", "jsonlite", "magrittr", "pegas", "Peptides", "RecordLinkage", "shazam", "seqinr", "stringdist", "textmineR", "yaml"))

source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")


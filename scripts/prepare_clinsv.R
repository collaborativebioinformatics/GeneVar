library(dplyr)

## input arguments
args = commandArgs(TRUE)
## args = c('nstd102.GRCh38.variant_call.tsv.gz', 'clinsv.tsv.gz')
clinsv.f = args[1]
out.f = args[2]


clin.df = read.table(clinsv.f, as.is=TRUE, skip=1, comment='', sep='\t', header=TRUE)
head(clin.df)

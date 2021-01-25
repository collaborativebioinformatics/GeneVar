library(dplyr)

## input arguments
args = commandArgs(TRUE)
## args = c('gnomad_v2.1_sv.sites.bed.gz', 'gnomad.af.tsv')
gnomad.f = args[1]
out.f = args[2]

message('read gnomad BED and select/rename a few columns...')
gnomad = read.table(gnomad.f, as.is=TRUE, header=TRUE, comment.char = '')
gnomad = gnomad %>% select(X.chrom, start, end, svtype, AF) %>%
  rename(chr=X.chrom, type=svtype, af=AF)
  
message('rename SV types into simpler types...')
table(gnomad$type)
gnomad = gnomad %>% mutate(svtype=ifelse(grepl('CN', type), 'CNV', type),
                           svtype=ifelse(grepl('INS', type), 'INS', svtype),
                           svtype=ifelse(grepl('DEL', type), 'DEL', svtype),
                           svtype=ifelse(grepl('DUP', type), 'DUP', svtype),
                           svtype=ifelse(grepl('INV', type), 'INV', svtype))
head(gnomad)

message('extract maximum allele frequency...')
gnomad$af.max = unlist(lapply(strsplit(gnomad$af, split=','), function(x) max(as.numeric(x))))
gnomad$af.max.no2 = unlist(lapply(strsplit(gnomad$af, split=','), function(x){
  if(length(x)>=3) x = x[-3]
  max(as.numeric(x))
}))
gnomad = gnomad %>% mutate(af=ifelse(svtype=='CNV', af.max.no2, af.max)) %>%
  select(chr, start, end, svtype, af)
head(gnomad)

message('write output TSV file...')
write.table(gnomad, file=out.f, sep='\t', row.names=FALSE, quote=FALSE)

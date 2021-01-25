library(dplyr)
library(VariantAnnotation)
library(GenomicRanges)

## input arguments
args = commandArgs(TRUE)
## args = c('dbvar38.2.tsv.gz', 'dbvar37.2.id.af.tsv.gz', 'clinsv.tsv.gz', 'clinvar.vcf.gz', 'test.tsv')
dbvar.f = args[1]
idaf.f = args[2]
clinsv.f = args[3]
clinvar.f = args[4]
out.f = args[5]

message('read dbvar...')
vars.df = read.table(dbvar.f, as.is=TRUE, header=TRUE)
head(vars.df)

message('read allele frequencies...')
idaf = read.table(idaf.f, as.is=TRUE, header=TRUE)

message('add frequencies...')
vars.df = merge(vars.df, idaf, all.x=TRUE) %>% dplyr::select(chr, start, end, variant_id, type, af) %>%
  mutate(af=ifelse(is.na(af), 0, af))
head(vars.df)

message('read clinical SVs...')
clinsv = read.table(clinsv.f, as.is=TRUE, sep='\t')
colnames(clinsv) = c('variant_id', 'clinsig')
clinsv = clinsv %>% mutate(clinsig=factor(clinsig, levels=c("Pathogenic", "Likely pathogenic", "Benign", "Benign/Likely benign", "Likely benign"))) %>%
  filter(!is.na(clinsig)) %>%
  arrange(clinsig) %>%
  group_by(variant_id) %>% summarize(clinsig=paste0(clinsig, collapse=';'))
vars.df = merge(vars.df, clinsv, all.x=TRUE) %>% dplyr::select(chr, start, end, everything())
head(vars.df)

message('read clinvar VCF...')
clinvar = readVcf(clinvar.f)
rowRanges(clinvar)$clinsig = info(clinvar)$CLNSIG
clinvar = rowRanges(clinvar)
clinvar$clinsig = unlist(lapply(clinvar$clinsig, '[', 1))

message('overlap with clinvar SNV...')
vars.gr = makeGRangesFromDataFrame(vars.df)
vars.clinvar = findOverlaps(vars.gr, clinvar) %>% as.data.frame %>%
  mutate(clinsig=clinvar$clinsig[subjectHits],
         clinsig=factor(clinsig, levels=c("Pathogenic", "Likely pathogenic", "Benign", "Benign/Likely benign", "Likely benign"))) %>%
  filter(!is.na(clinsig)) %>%
  group_by(queryHits, clinsig) %>% summarize(n=n()) %>%
  arrange(clinsig) %>%
  group_by(queryHits) %>% summarize(clinsig=paste0(clinsig, '(', n, ')', collapse=';'))

vars.df$clin.snv = 0
vars.df$clin.snv[vars.clinvar$queryHits] = vars.clinvar$clinsig
head(vars.df)

message('sort...')
vars.df = vars.df[order(vars.df$chr, vars.df$start, vars.df$end),]

write.table(vars.df, file=out.f, sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)

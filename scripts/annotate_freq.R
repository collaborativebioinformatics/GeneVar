library(dplyr)
library(GenomicRanges)

## input arguments
args = commandArgs(TRUE)
## args = c('dbvar37.2.tsv', 'gnomad.af.tsv', 'dbvar37.2.id.af.tsv.gz')
dbvar.f = args[1]
gnomad.f = args[2]
out.f = args[3]

message('read gnomadsv...')
gnomad = read.table(gnomad.f, as.is=TRUE, header=TRUE)
head(gnomad)

message('read dbvar...')
vars.df = read.table(dbvar.f, as.is=TRUE, header=TRUE)
head(vars.df)

message('match SVs per type...')
MIN.ROL = .9
vars.af = lapply(unique(vars.df$type), function(tt){
  vars.gr = vars.df %>% filter(type==tt) %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  gnomad.gr = gnomad %>% filter(svtype==tt) %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  findOverlaps(vars.gr, gnomad.gr) %>% as.data.frame %>%
    mutate(af=gnomad.gr$af[subjectHits],
           variant_id=vars.gr$variant_id[queryHits],
           qw=width(vars.gr)[queryHits], sw=width(gnomad.gr)[subjectHits],
           qsw=width(pintersect(vars.gr[queryHits], gnomad.gr[subjectHits]))) %>%
    filter(qsw>MIN.ROL * qw,
           qsw > MIN.ROL * sw) %>% 
    group_by(variant_id) %>% summarize(af=max(af))
}) %>% bind_rows()

head(vars.af)
summary(vars.af$af)

write.table(vars.af, file=out.f, sep='\t', row.names=FALSE, quote=FALSE)

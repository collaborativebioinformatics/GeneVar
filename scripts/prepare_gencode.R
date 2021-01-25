library(rtracklayer)
library(GenomicRanges)
library(dplyr)

## input arguments
args = commandArgs(TRUE)
## args = c('gencode.v35.annotation.gff3.gz', 'gencode.tsv')
genc.f = args[1]
out.f = args[2]

message('read gencode annotation...')
genc = import(genc.f)
genc

message('susbet to information of interest...')
genc = subset(genc, type %in% c('gene', 'transcript', 'CDS', 'five_prime_UTR', 'three_prime_UTR'))
mcols(genc) = mcols(genc)[,c('type', 'gene_name', 'gene_type', 'gene_id', 'transcript_id', 'exon_number')]
genc

genc = genc %>% as.data.frame %>% mutate(chr=gsub('chr', '', seqnames)) %>%
  select(chr, start, end, strand, type, gene_name, gene_id, gene_type, transcript_id, exon_number)
write.table(genc, file=out.f, row.names=FALSE, quote=FALSE, sep='\t')


## vars = read.table('all_variants_chr21.tsv', as.is=TRUE, sep='\t', header=TRUE)
## vars.gr = vars %>% mutate(chr=paste0('chr', chr)) %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE)
## vars.gr
## var.gene.df = findOverlaps(vars.gr, genc) %>% as.data.frame %>%
##     filter(genc$gene_type[subjectHits] == 'protein_coding') %>%
##     mutate(variant_id=vars.gr$variant_id[queryHits],
##            gene_id=genc$gene_id[subjectHits],
##            gene_name=genc$gene_name[subjectHits],
##            transcript_id=genc$transcript_id[subjectHits],
##            exon_number=genc$exon_number[subjectHits], type=genc$type[subjectHits]) %>% 
##     select(-queryHits, -subjectHits) %>% unique
## head(var.gene.df)
## nrow(var.gene.df)
## var.gene.df = var.gene.df %>% group_by(variant_id, gene_id, gene_name, transcript_id, type) %>% summarize(exon_number=paste(sort(unique(exon_number)), collapse=';'))
## head(var.gene.df)
## nrow(var.gene.df)

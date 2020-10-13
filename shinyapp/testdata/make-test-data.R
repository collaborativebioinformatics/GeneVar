library(dplyr)

## genes
genes = paste0('gene', 1:100)

## variants
N = 1000
vars = tibble(variant_id=paste0('sv', 1:N),
              type=sample(c('deletion', 'duplication'), N, TRUE),
              chr=sample(1:22, N, TRUE),
              start=round(runif(N, 1, 1e8))) %>%
  mutate(end=start + round(runif(N, 50, 1e6)))
head(vars)
write.table(vars, file='all_variants.tsv', quote=FALSE, row.names=FALSE, sep='\t')

## variant-gene pair including gene impact
gene.var = tibble(variant_id=sample(vars$variant_id, 500, TRUE),
                  gene_id=sample(genes, 500, TRUE),
                  elt_type=sample(c('exon', 'UTR'), 500, TRUE)) %>%
  mutate(elt_info=ifelse(elt_type=='exon',
                         paste('exon', sample(1:20, n(), TRUE)),
                         sample(c("5p UTR", "3p UTR"), n(), TRUE)))
head(gene.var)
write.table(gene.var, file='gene_variants.tsv', quote=FALSE, row.names=FALSE, sep='\t')

## overlap clinically important SV
clinsv = tibble(variant_id=sample(vars$variant_id, 100),
                pathogenic_clinvar_sv=TRUE)
write.table(clinsv, file='clinsv_variants.tsv', quote=FALSE, row.names=FALSE, sep='\t')

## overlap clinically important SNV/indels
clinsnv = tibble(variant_id=sample(vars$variant_id, 200),
                 pathogenic_clinvar_snv_indel=TRUE)
write.table(clinsnv, file='clinsnv_variants.tsv', quote=FALSE, row.names=FALSE, sep='\t')

## allele frequency
af = tibble(variant_id=sample(vars$variant_id, 300), af=runif(300))
write.table(af, file='af.tsv', quote=FALSE, row.names=FALSE, sep='\t')

## external resources
ext.urls = tibble(gene_id=genes, omim_url='https://www.youtube.com/watch?v=dQw4w9WgXcQ')
write.table(ext.urls, file='ext_urls.tsv', quote=FALSE, row.names=FALSE, sep='\t')

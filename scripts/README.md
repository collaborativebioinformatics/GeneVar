## Trace

```
snakemake -n dbvar38.ann.tsv.gz
```

## Subset dbVar variant to one chromosome

We used a subset of the data to develop the app: variants and genes in chr21.
The different dbVar versions (GRCh37 and GRCh38) and genes were subset as shown in the [extract-chr21-genes-variants.ipynb](extract-chr21-genes-variants.ipynb) notebook.

## Allele frequency

### Extract allele frequency from gnomAD-SV

The [AF_extract.py](AF_extract.py) script reads a gnomAd_SV vcf file and extracts the allele frequencies of all variants and outputs them into a CSV file

Usage:

`python AF_extract.py -vcf <path to the gnomAD_SV vcf file> -out <path to the output file>`

### Match allele frequencies with dbVar SVs

The [gnomad-variants-af.ipynb](gnomad-variants-af.ipynb) notebook shows how we matched variants in dbVar (GRCh37) and the gnomAD-SV. 
Briefly: 50% reciprocal overlap per SV type.

R/Bioconductor packages used include:
- GenomicRanges
- dplyr

## Variants overlapping genes and gene impact

### Script to extract information for one gene

We used [annotGeneSV.R](annotGeneSV.R) (at first [linkGeneWithSV.R](linkGeneWithSV.R)) to read variants from dbVar and Gencode annotation and extract relevant information.

### Overlapping all genes with all variants in chr21

The [variant-gene-overlap.ipynb](variant-gene-overlap.ipynb) notebook shows how we overlapped variants in chr21 with Gencode to extract gene impact information.

R/Bioconductor packages used include:
- rtracklayer
- GenomicRanges
- dplyr

## Overlap dbVar SVs with SNVs of known clinical significance

We overlapped SVs with ClinVar using bedtools. 
The overlaps were then summarize in a TSV file using code in the [make-misc-tsvs.ipynb](make-misc-tsvs.ipynb) notebook.

## Annotate dbVar SVs with clinical significance

We used study nstd102 that contain clinical SVs.
The code to make the TSV for this annotation is part of the [make-misc-tsvs.ipynb](make-misc-tsvs.ipynb) notebook.

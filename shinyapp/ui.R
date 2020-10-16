library(shiny)
library(DT)
library(shinydashboard)

data = list()
data$all_variants = read.table('./all_variants_chr21.tsv', as.is=TRUE, header=TRUE, sep='\t')
data$gene_variants = read.table('./variants.genes.chr21.tsv.gz', as.is=TRUE, header=TRUE, sep='\t')

## list all genes
genes = unique(c(data$gene_variants$gene_id, data$gene_variants$gene_name, data$gene_variants$transcript_id))
genes = unique(c('ENST00000400454.6', genes))

## merge general variant info
svtypes = sort(unique(data$all_variants$type))
svsize.max = max(data$all_variants$end - data$all_variants$start)

ui <- dashboardPage(
  dashboardHeader(title='GeneVar'),
  dashboardSidebar(
    selectizeInput('gene_search', 'Gene', genes),
    div(p(' Search by gene name, gene id (ENSG...),, or transcript ID (ENST...)'),
        p(' Examples:'),
        p('  DSCAM'),
        p('  ENSG00000171587.15'),
        p('  ENST00000400454.6')),
    checkboxGroupInput('svtypes', "SV type", svtypes, svtypes),
    numericInput('size.min', 'Minimum SV size (bp)', 0, 0),
    numericInput('size.max', 'Maximum SV size (bp)', svsize.max, svsize.max)
  ),
  dashboardBody(
    htmlOutput('title'),
    fluidRow(
      ## A static infoBox
      infoBoxOutput("sv_box"),
      infoBoxOutput("path_sv_box"),
      infoBoxOutput("path_snv_box")
    ),
    shiny::htmlOutput('omim_url', class='btn btn-default action-button shiny-bound-input'),
    shiny::htmlOutput('gtex_url', class='btn btn-default action-button shiny-bound-input'),
    shiny::htmlOutput('gnomad_url', class='btn btn-default action-button shiny-bound-input'),
    hr(),
    dataTableOutput('vars_table'),
    hr(),
    h2('Allele frequency distribution'),
    plotOutput('af_plot')
  )
)

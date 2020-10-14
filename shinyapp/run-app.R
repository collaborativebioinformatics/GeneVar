library(dplyr)
library(shiny)
library(DT)
library(ggplot2)
library(shinydashboard)

## read data (for now from the testdata folder) into a list
## data.files = c('af.tsv', 'all_variants.tsv', 'clinsnv_variants.tsv',
##                'clinsv_variants.tsv', 'ext_urls.tsv', 'gene_variants.tsv')
## data = lapply(data.files, function(ff){
##   read.table(paste0('testdata/', ff), as.is=TRUE, header=TRUE, sep='\t')
## })
## names(data) = gsub('.tsv', '', data.files)

## read real data
data = list()
data$all_variants = read.table('./all_variants_chr21.tsv', as.is=TRUE, header=TRUE, sep='\t')
data$gene_variants = read.table('./variants.genes.chr21.tsv.gz', as.is=TRUE, header=TRUE, sep='\t')
data$clinsv_variants = read.table('./sv.clinical.variants.chr21.nstd120.tsv', as.is=TRUE, header=TRUE, sep='\t')
data$clinsnv_variants = read.table('./final_clinvar_dbvar_results_summary.txt', as.is=TRUE, header=TRUE, sep='\t')
data$af = read.table('./gnomad-af-variants-chr21.tsv', as.is=TRUE, header=TRUE, sep='\t')

## save pathogenic variants and modify the main table to have a summary
snvs.path = data$clinsnv_variants %>% filter(effect %in% c('Pathogenic', 'Likely pathogenic')) %>% .$variant_id
data$clinsnv_variants = data$clinsnv_variants %>%
  mutate(effect=factor(effect, levels=c("Pathogenic", "Likely pathogenic", "Benign", "Benign/Likely benign", "Likely benign"))) %>%
  filter(!is.na(effect)) %>% 
  arrange(effect) %>% group_by(variant_id) %>% summarize(clinical_snv=paste0(effect, '(', n, ')', collapse=';'))

## same for SVs
svs.path = data$clinsv_variants %>% filter(sv_clinical_significance %in% c('Pathogenic', 'Likely pathogenic')) %>% .$variant_id
data$clinsv_variants = data$clinsv_variants %>%
  mutate(effect=factor(sv_clinical_significance,
                       levels=c("Pathogenic", "Likely pathogenic", "Benign", "Benign/Likely benign", "Likely benign"))) %>%
  filter(!is.na(effect)) %>%
  group_by(variant_id, effect) %>% summarize(n=n()) %>% 
  arrange(effect) %>% group_by(variant_id) %>% summarize(clinical_sv=paste0(effect, '(', n, ')', collapse=';'))

## list all genes
genes = unique(c(data$gene_variants$gene_id, data$gene_variants$gene_name, data$gene_variants$transcript_id))
genes = unique(c('ENST00000400454.6', genes))

## merge general variant info
vars.df = data$all_variants %>%
  merge(data$clinsv_variants, all.x=TRUE) %>%
  merge(data$clinsnv_variants, all.x=TRUE) %>%
  merge(data$af, all.x=TRUE)

## format columns for better DataTable experience
vars.df = vars.df %>% mutate(type=factor(type), coord=paste0(chr, ':', start, '-', end), size=end-start) %>%
  select(-chr, -start, -end) %>%
  select(variant_id, coord, type, size, af, everything())
svtypes = sort(unique(vars.df$type))
svsize.max = max(vars.df$size)

## function to add links to the variant table
dtify <- function(df){
  df %>% mutate(variant_id=paste0('<a href="https://www.ncbi.nlm.nih.gov/dbvar/variants/', variant_id, '" target="_blank">', variant_id, '</a>'))
}

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

## server side of the app
server <- function(input, output) {
  ## reactive conductor to apply the filtering only once for all elements that need it
  selVars <- reactive({
    message('Gene: ', input$gene_search)
    gene.var = data$gene_variants %>%
      filter(gene_id==input$gene_search | gene_name==input$gene_search | transcript_id==input$gene_search) %>%
      group_by(variant_id, type) %>%
      summarize(exon_number=ifelse(any(exon_number!=''), paste0(sort(unique(exon_number)), collapse='|'), '')) %>% 
      mutate(elt_type=ifelse(exon_number!='', paste0(type, '(', exon_number, ')'), type)) %>%
      select(-type)      
    vars.sel = vars.df %>% filter(variant_id %in% gene.var$variant_id,
                                  type %in% input$svtypes,
                                  size >= input$size.min,
                                  size <= input$size.max)
    vars.sel = merge(gene.var, vars.sel, all.y=TRUE) %>%
      group_by(variant_id, coord, type, size, af, clinical_sv, clinical_snv) %>%
      summarize(gene_impact=paste(sort(unique(elt_type)), collapse=';')) %>% ungroup
    vars.sel %>% arrange(desc(clinical_sv))
  })
  geneName <- reactive({
    gene_name = input$gene_search
    if(all(gene_name != data$gene_variants$gene_name)){
      gene.var = data$gene_variants %>%
        filter(gene_id==input$gene_search | gene_name==input$gene_search | transcript_id==input$gene_search)
      gene_name = gene.var$gene_name[1]
    }
    return(gene_name)
  })
  ## Text
  output$title = renderText({
    paste0('<h1>', input$gene_search, '</h1>')
  })
  output$omim_url = renderText({
    return(as.character(a('OMIM', href=paste0('https://www.genenames.org/tools/search/#!/genes?query=', geneName()), target='_blank')))
  })
  output$gtex_url = renderText({
    return(as.character(a('GTEx', href=paste0('https://gtexportal.org/home/gene/', geneName()), target='_blank')))
  })
  output$gnomad_url = renderText({
    return(as.character(a('gnomAD', href=paste0('https://gnomad.broadinstitute.org/gene/', geneName()), target='_blank')))
  })
  ## boxes
  output$sv_box <- renderInfoBox({
    infoBox("SVs", nrow(selVars()), icon=icon("dna"), color="blue")
  })
  output$path_sv_box <- renderInfoBox({
    infoBox("Clinical SVs", sum(selVars()$variant_id %in% svs.path),
            icon=icon("stethoscope"), color="red")
  })
  output$path_snv_box <- renderInfoBox({
    infoBox("Overlap Clinical SNVs", sum(selVars()$variant_id %in% snvs.path),
            icon=icon("stethoscope"), color="red")
  })
  ## dynamic tables
  output$vars_table <- renderDataTable(
    datatable(dtify(selVars()),
    filter='top',
    rownames=FALSE,
    escape=FALSE,
    options=list(pageLength=15, searching=FALSE)))
  ## Graph
  output$af_plot = renderPlot({
    ggplot(selVars(), aes(x=af)) + geom_histogram() + theme_bw() + xlab('allele frequency') + xlim(-.1,1.1)
  })
}

#### launch app locally
runApp(list(ui=ui, server=server))

#### launch app on UCSC server
## runApp(list(ui=ui, server=server), port=3457, host='0.0.0.0')

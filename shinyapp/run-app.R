library(dplyr)
library(shiny)
library(DT)
library(ggplot2)
library(shinydashboard)

## read data (for now from the testdata folder) into a list
data.files = c('af.tsv', 'all_variants.tsv', 'clinsnv_variants.tsv',
               'clinsv_variants.tsv', 'ext_urls.tsv', 'gene_variants.tsv')
data = lapply(data.files, function(ff){
  read.table(paste0('testdata/', ff), as.is=TRUE, header=TRUE, sep='\t')
})
names(data) = gsub('.tsv', '', data.files)

## list all genes
gene.df = tibble(gene_id=unique(data$ext_urls$gene_id, data$gene_variants$gene_id))
## add number of SVs
gene.df = data$gene_variants %>% group_by(gene_id) %>% summarize(nb.svs=n()) %>%
  merge(gene.df, all.y=TRUE) %>% mutate(nb.svs=ifelse(is.na(nb.svs), 0, nb.svs))

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
    selectizeInput('gene_search', 'Gene', gene.df$gene_id),
    checkboxGroupInput('svtypes', "SV type", svtypes, svtypes),
    numericInput('size.min', 'Minimum SV size (bp)', 0, 0),
    numericInput('size.max', 'Maximum SV size (bp)', svsize.max, svsize.max)
  ),
  dashboardBody(
    htmlOutput('title'),
    fluidRow(
      ## A static infoBox
      infoBoxOutput("sv_box"),
      infoBoxOutput("path_sv_box")
    ),
    shiny::htmlOutput('omim_url', class='btn btn-default action-button shiny-bound-input'),
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
    gene.var = data$gene_variants %>% filter(gene_id==input$gene_search)
    vars.sel = vars.df %>% filter(variant_id %in% gene.var$variant_id,
                                  type %in% input$svtypes,
                                  size >= input$size.min,
                                  size <= input$size.max)
    vars.sel
  })
  ## Text
  output$title = renderText({
    paste0('<h1>', input$gene_search, '</h1>')
  })
  output$omim_url = renderText({
    urls = data$ext_urls %>% filter(gene_id==input$gene_search)
    if(!is.null(urls$omim_url)){
      return(as.character(a('OMIM', href=urls$omim_url, target='_blank')))
    }
    return('')
  })
  ## boxes
  output$sv_box <- renderInfoBox({
    infoBox("SVs", nrow(selVars()), icon=icon("dna"), color="blue")
  })
  output$path_sv_box <- renderInfoBox({
    infoBox("Clinical SVs", sum(selVars()$pathogenic_clinvar_sv, na.rm=TRUE),
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
    ggplot(selVars(), aes(x=af)) + geom_histogram() + theme_bw() + xlab('allele frequency')
  })
}

#### launch app locally
runApp(list(ui=ui, server=server))

#### launch app on UCSC server
## runApp(list(ui=ui, server=server), port=3457, host='0.0.0.0')

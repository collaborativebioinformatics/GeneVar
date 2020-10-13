library(dplyr)
library(shiny)
library(DT)
library(ggplot2)

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
vars.df = vars.df %>% mutate(type=factor(type), chr=factor(chr))

## function to compute gene-level stats and return HTML code


## app interface
ui <- fluidPage(
  titlePanel("GeneVar"),
  sidebarLayout(
    sidebarPanel(width=3, dataTableOutput('gene_search_table')),
    mainPanel(
      width=9,
      fluidRow(htmlOutput('title'), br(), htmlOutput('gene_info')),
      hr(),
      fluidRow(dataTableOutput('vars_table')),
      hr(),
      fluidRow(htmlOutput('af_pre'), plotOutput('af_plot'))
    )
  )
)

## server side of the app
server <- function(input, output) {
  ## dynamic tables
  output$gene_search_table <- renderDataTable(
    gene.df,
    filter='top',
    rownames=FALSE,
    options=list(pageLength=15),
    selection='single')
  output$vars_table <- renderDataTable(
    datatable({
      ## if no gene was select, don't try to make a table
      vars.sel = tibble()
      if(!is.null(input$gene_search_table_rows_selected)){
        sel.ii = input$gene_search_table_rows_selected
        gene.sel = gene.df$gene_id[sel.ii]
        message('Gene: ', gene.sel)
        gene.var = data$gene_variants %>% filter(gene_id==gene.sel)
        vars.sel = vars.df %>% filter(variant_id %in% gene.var$variant_id)
      }
      vars.sel
    },
    filter='top',
    rownames=FALSE,
    options=list(pageLength=15, searching=FALSE),
    selection='single'))
  ## Text
  output$title = renderText({
    sel.ii = input$gene_search_table_rows_selected
    gene.sel = gene.df$gene_id[sel.ii]
    paste0('<h1>', gene.sel, '</h1>')
  })
  output$af_pre = renderText({
    if(!is.null(input$gene_search_table_rows_selected)){
      return(as.character(h1('Allele frequency distribution')))
    }
    return(NULL)
  })
  output$gene_info = renderText({
    if(!is.null(input$gene_search_table_rows_selected)){
      sel.ii = input$gene_search_table_rows_selected
      gene.sel = gene.df$gene_id[sel.ii]
      gene.var = data$gene_variants %>% filter(gene_id==gene.sel)
      return(as.character(p('Number of SVs: ', nrow(gene.var))))
    }
    return(NULL)
  })
  ## Graph
  output$af_plot = renderPlot({
    ggp = NULL
    if(!is.null(input$gene_search_table_rows_selected)){
      sel.ii = input$gene_search_table_rows_selected
      gene.sel = gene.df$gene_id[sel.ii]
      message('Gene: ', gene.sel)
      gene.var = data$gene_variants %>% filter(gene_id==gene.sel)
      vars.sel = vars.df %>% filter(variant_id %in% gene.var$variant_id)
      ggp = ggplot(vars.sel, aes(x=af)) + geom_histogram() + theme_bw() + xlab('allele frequency')
    }
    ggp
  })
}

#### launch app locally
runApp(list(ui=ui, server=server))

#### launch app on UCSC server
## runApp(list(ui=ui, server=server), port=3457, host='0.0.0.0')

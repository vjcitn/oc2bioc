validate_tsv = function(fn) {
  tmp = read.delim(fn, sep="\t", h=FALSE)
  nr = nrow(tmp)
  nc = grep("chr", tmp[,1])
  stopifnot(length(nc)==nr)
  stopifnot(ncol(tmp)==6)
  return(TRUE)
}

validate_tab = function(tab) {
  nr = nrow(tab)
  nc = grep("chr", tab[,1])
  stopifnot(length(nc)==nr)
  stopifnot(ncol(tab)==6)
  return(TRUE)
}

#' simple app
#' @import shiny
#' @import plotly
#' @import ggplot2
#' @export
ocapp = function() {
 ui = fluidPage(
  sidebarLayout(
   sidebarPanel(
    helpText("Variant Annotation with OpenCRAVAT"),
    checkboxGroupInput("build", "ref. build", choices=c("hg38", "hg19"), selected="hg38"),
    fileInput("basicfile", "TSV file")
    ),
    mainPanel(
     tabsetPanel(
      tabPanel("GeneOnt",
       plotly::plotlyOutput("plot1")
       ),
      tabPanel("SeqOnt",
       plotly::plotlyOutput("plot2")
       ),
      tabPanel("Variants", DT::dataTableOutput("vartab")),
      tabPanel("Genes", DT::dataTableOutput("genetab"))
      )
    )
   )
  )
 server = function(input, output) {
  get_data = reactive({
   file = input$basicfile
   req(file)
   ext = tools::file_ext(file$datapath)
   validate(need(ext == "tsv", "Use a tsv file"))
   gg <- paste("cravat", file$datapath, "-t text -l", input$build)
   dat <- read.delim(file$datapath, h=FALSE, sep="\t")
   showNotification("starting cravat...")
   system(gg)
   showNotification("done.")
   list(tab=dat, sqlite=paste0(file$datapath, ".sqlite", sep=""))
   })
  output$vartab = DT::renderDataTable({
    con = DBI::dbConnect(RSQLite::SQLite(), get_data()$sqlite)
    get_oc_tab(con)
   }) # defaults
  output$genetab = DT::renderDataTable({
    con = DBI::dbConnect(RSQLite::SQLite(), get_data()$sqlite)
    tmp = get_oc_tab(con, "gene")
    tmp[which(tmp[,3]>0),]
   }) # defaults
  output$plot1 = plotly::renderPlotly({
   con = DBI::dbConnect(RSQLite::SQLite(), get_data()$sqlite)
   plotly::ggplotly(barplot_gene_ontology(con))
   })
  output$plot2 = plotly::renderPlotly({
   con = DBI::dbConnect(RSQLite::SQLite(), get_data()$sqlite)
   plotly::ggplotly(barplot_sequence_ontology(con))
   })
 }
 runApp(list(ui=ui, server=server))
}

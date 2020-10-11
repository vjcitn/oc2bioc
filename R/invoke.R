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
#' @rawNamespace import("plotly", except=c("config", "last_plot", "select"))
#' @import ggplot2
#' @param cravat_cmd character(1) string to invoke cravat
#' @param sqlite_to_home logical(1) if TRUE will copy generated sqlite to home, with basename of input filename as prefix to .sqlite
#' @export
ocapp = function(cravat_cmd="cravat", sqlite_to_home=TRUE) {
 ui = fluidPage(
  sidebarLayout(
   sidebarPanel(
    helpText("Variant Annotation with", a(href="https://github.com/KarchinLab/open-cravat/wiki", "OpenCRAVAT")),
    checkboxGroupInput("build", "ref. build", choices=c("hg38", "hg19"), selected="hg38"),
    fileInput("basicfile", "TSV file"), width=2
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
      tabPanel("Genes", DT::dataTableOutput("genetab")),
      tabPanel("Versions", verbatimTextOutput("vertab"))
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
   gg <- paste(cravat_cmd, file$datapath, "-t text -l", input$build)
   dat <- read.delim(file$datapath, h=FALSE, sep="\t")
   showNotification("starting cravat...")
   system(gg)
   showNotification("done.")
   if (sqlite_to_home) file.copy(paste0(file$datapath, ".sqlite", sep=""), paste("~/", file$name, ".sqlite", sep=""))
   list(tab=dat, sqlite=paste0(file$datapath, ".sqlite", sep=""))
   })
  output$vartab = DT::renderDataTable({
    con = DBI::dbConnect(RSQLite::SQLite(), get_data()$sqlite)
    tmp = get_oc_tab(con)
    if ("Link" %in% colnames(tmp)) {
     ok = which(!is.na(tmp$Link)); # make active links for mupit plots
     tmp$Link[ok] = paste("<a href='", tmp$Link[ok], "' target='_blank'>mupit</a>", sep="")
     }
    tmp
   }, escape=FALSE) # defaults
  output$genetab = DT::renderDataTable({
    con = DBI::dbConnect(RSQLite::SQLite(), get_data()$sqlite)
    tmp = get_oc_tab(con, "gene")
    tmp[which(tmp[,3]>0),]
   }) # defaults
  output$vertab = renderPrint({
    con = DBI::dbConnect(RSQLite::SQLite(), get_data()$sqlite)
    DBI::dbReadTable(con, "info")
   })
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

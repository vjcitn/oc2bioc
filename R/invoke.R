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
    uiOutput("selector"),
    checkboxGroupInput("build", "ref. build", choices=c("hg38", "hg19"), selected="hg38"),
    actionButton("stopBtn", "Stop app")
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
  observeEvent(input$stopBtn, { 
    con = DBI::dbConnect(RSQLite::SQLite(), get_data()$sqlite)
    variants = get_oc_tab(con)
    genes = get_oc_tab(con, "gene")
    stopApp(returnValue=list(variants=variants, genes=genes))
    })
  output$selector = renderUI({
    tsvs = dir("~", patt="tsv$")
    tsvsf = dir("~", patt="tsv$", full=TRUE)
    # should validate what is found
    validate(need(length(tsvs)>0, "no tsv files found"))
    checkboxGroupInput("picked_tsv", "Choose a TSV for variant annotation", choiceNames  = tsvs, choiceValues=tsvsf)
    })
  get_data = reactive({
   tf = paste0(tempfile(), ".tsv")
   req(input$picked_tsv)
   file.copy(input$picked_tsv, tf)  # the cravat run occurs in temp area
   gg <- paste(cravat_cmd, tf, "-t text -l", input$build)
   showNotification("starting cravat...")
   system(gg)
   showNotification("done.")
   if (sqlite_to_home) file.copy(paste0(tf, ".sqlite", sep=""), 
        paste("~/", basename(input$picked_tsv), ".sqlite", sep=""))
   list(tab=NA, sqlite=paste0(tf, ".sqlite", sep=""))
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

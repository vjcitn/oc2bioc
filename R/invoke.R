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
#' @import shinytoastr
#' @rawNamespace import("plotly", except=c("config", "last_plot", "select"))
#' @import ggplot2
#' @param cravat_cmd character(1) string to invoke cravat
#' @param sqlite_to_home logical(1) if TRUE will copy generated sqlite to home, with basename of input filename as prefix to .sqlite
#' @export
ocapp = function(cravat_cmd="cravat", sqlite_to_home=TRUE) {
 ui = fluidPage(
  useToastr(),
  sidebarLayout(
   sidebarPanel(
    helpText("Variant Annotation with", a(href="https://github.com/KarchinLab/open-cravat/wiki", "OpenCRAVAT")),
    uiOutput("selector"),
    checkboxGroupInput("build", "ref. build", choices=c("hg38", "hg19"), selected="hg38"),
    actionButton("stopBtn", "Stop app")
    ),
    mainPanel(
     tabsetPanel(
      tabPanel("SeqOnt",
       plotly::plotlyOutput("plot2")
       ),
      tabPanel("GeneOnt",
       plotly::plotlyOutput("plot1")
       ),
      tabPanel("Variants", DT::dataTableOutput("vartab")),
      tabPanel("Genes", DT::dataTableOutput("genetab")),
      tabPanel("Versions", DT::dataTableOutput("vertab")),
      tabPanel("Log", htmlOutput("textstream_output"))
      )
    )
   )
  )
 server = function(input, output, session) {
  thelog <<- NULL
# https://stackoverflow.com/questions/50650616/stream-system-output-to-shiny-front-end-continuously
   rv <- reactiveValues(textstream = c(""),
                       timer = reactiveTimer(1000),
                       started = TRUE)
   observe({
     rv$timer()
     req(thelog)
     rv$textstream = paste(readLines(thelog), collapse="<br/>")
     })
   output$textstream_output <- renderUI({
     HTML(rv$textstream)
   })
  observeEvent(input$stopBtn, { 
    con = DBI::dbConnect(RSQLite::SQLite(), get_data()$sqlite)
    variants = get_oc_tab(con)
    genes = get_oc_tab(con, "gene")
    DBI::dbDisconnect(con)
    stopApp(returnValue=list(variants=variants, genes=genes))
    })
  output$selector = renderUI({
    tsvs = dir("~", patt="tsv$")
    tsvsf = dir("~", patt="tsv$", full=TRUE)
    # should validate what is found
    validate(need(length(tsvs)>0, "no tsv files found"))
    checkboxGroupInput("picked_tsv", "Choose a TSV for variant annotation", choiceNames  = tsvs, choiceValues=tsvsf)
    })
  output$GOcomment = renderPrint({
    req(input$picked_tsv)
    paste("OpenCRAVAT GO analysis of variants in ", input$picked_tsv)
    })
  get_data = reactive({
   tf = paste0(tempfile(), ".tsv")
   req(input$picked_tsv)
   validate(need(length(input$picked_tsv)==1, "check only one box"))
   file.copy(input$picked_tsv, tf)  # the cravat run occurs in temp area
   gg <- paste(cravat_cmd, tf, "-t text -l", input$build)
   showNotification("starting cravat...")
   toastr_info("running cravat...")
   system(gg)
   thelog <<- paste0(tf, ".log")
   showNotification("done.")
   if (sqlite_to_home) file.copy(paste0(tf, ".sqlite", sep=""), 
        paste("~/", basename(input$picked_tsv), ".sqlite", sep=""))
   list(tab=NA, sqlite=paste0(tf, ".sqlite", sep=""))
   })
  output$vartab = DT::renderDataTable({
    con = DBI::dbConnect(RSQLite::SQLite(), get_data()$sqlite)
    tmp = get_oc_tab(con)
    DBI::dbDisconnect(con)
    if ("Link" %in% colnames(tmp)) {
     ok = which(!is.na(tmp$Link)); # make active links for mupit plots
     tmp$Link[ok] = paste("<a href='", tmp$Link[ok], "' target='_blank'>mupit</a>", sep="")
     }
    tmp
   }, escape=FALSE) # defaults
  output$genetab = DT::renderDataTable({
    con = DBI::dbConnect(RSQLite::SQLite(), get_data()$sqlite)
    tmp = get_oc_tab(con, "gene")
    DBI::dbDisconnect(con)
    tmp[which(tmp[,3]>0),]
   }) # defaults
  output$vertab = DT::renderDataTable({
    con = DBI::dbConnect(RSQLite::SQLite(), get_data()$sqlite)
    tmp <- DBI::dbReadTable(con, "info")
    DBI::dbDisconnect(con)
    tmp
   })
  output$plot1 = plotly::renderPlotly({
   con = DBI::dbConnect(RSQLite::SQLite(), get_data()$sqlite)
   plotly::ggplotly(barplot_gene_ontology(con) + labs(title=paste("GO analysis of variants in",
      basename(input$picked_tsv)), subtitle="(15 most frequent categories among mutated genes)"))
   })
  output$plot2 = plotly::renderPlotly({
   con = DBI::dbConnect(RSQLite::SQLite(), get_data()$sqlite)
   plotly::ggplotly(barplot_sequence_ontology(con) + labs(title=paste("Sequence Ontology classes of variants in",
     basename(input$picked_tsv))))
   })
 }
 runApp(list(ui=ui, server=server))
}

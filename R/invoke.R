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
      tabPanel("About", helpText("Source code at github.com/vjcitn/oc2bioc"), 
                 helpText("Gene Ontology table is constructed using 15 most frequent
                 occurrence of GO tags in MF or BP subontologies
                 for which association in humans is annotated with evidence code among
                 TAS, IPI, IDA, IMP, EXP, IGI.  Mapping uses org.Hs.eg.db and GO.db.  Counts are truncated at 1000.")),
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
    info = get_oc_tab(con, "info", use_header=FALSE)
    DBI::dbDisconnect(con)
    stopApp(returnValue=list(variants=variants, genes=genes, info=info))
    })
  output$selector = renderUI({
    tsvs = dir("~", patt="tsv$|sqlite$")
    tsvsf = dir("~", patt="tsv$|sqlite$", full=TRUE)
    # should validate what is found
    validate(need(length(tsvs)>0, "no tsv files found"))
    checkboxGroupInput("picked_resource", "Choose a TSV for variant annotation, or SQLite for exploration", choiceNames  = tsvs, choiceValues=tsvsf)
    })
  output$GOcomment = renderPrint({
    req(input$picked_resource)
    paste("OpenCRAVAT GO analysis of variants in ", input$picked_resource)
    })
  get_data = reactive({
   tf = paste0(tempfile(), ".tsv")
   req(input$picked_resource)
   validate(need(length(input$picked_resource)==1, "check only one box"))
   if (length(grep("sqlite$", input$picked_resource))==0) {
    file.copy(input$picked_resource, tf)  # the cravat run occurs in temp area
    gg <- paste(cravat_cmd, tf, "-t text -l", input$build)
    showNotification("starting cravat...")
    toastr_info("running cravat...")
    system(gg)
    thelog <<- paste0(tf, ".log")
    showNotification("done.")
    tmp_sqlite_name = paste0(tf, ".sqlite", sep="")
    if (sqlite_to_home) file.copy(tmp_sqlite_name,
        paste("~/", basename(input$picked_resource), ".sqlite", sep=""))
    sqlite_name = tmp_sqlite_name
    } else sqlite_name = input$picked_resource
   list(tab=NA, sqlite=sqlite_name)
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
    tmp  # [which(tmp[,3]>0),]
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
      basename(input$picked_resource)), subtitle="(15 most frequent categories among mutated genes)"))
   })
  output$plot2 = plotly::renderPlotly({
   con = DBI::dbConnect(RSQLite::SQLite(), get_data()$sqlite)
   plotly::ggplotly(barplot_sequence_ontology(con) + labs(title=paste("Sequence Ontology classes of variants in",
     basename(input$picked_resource))))
   })
 }
 runApp(list(ui=ui, server=server))
}

#' var_in_tx: example output of queryOC -- variant in transcribed region
#' @docType data
#' @format httr GET result
#' @source demo run
#' @note The value component of 'cookies' was deleted.
#' @examples
#' var_in_tx
#' str(httr::content(var_in_tx))
"var_in_tx"

#' nonco_var: example output of queryOC -- variant in non-coding region
#' @docType data
#' @format httr GET result
#' @source demo run
#' @note The value component of 'cookies' was deleted.
#' @examples
#' names(httr::content(nonco_var))
"nonco_var"

#' a chain object for transforming hg19 to GRCh38 coordinates
#' @docType data
#' @format rtracklayer chain structure
#' @examples
#' ch19to38
"ch19to38"

#' mutations (SNV) from TCGA ACC
#' @docType data
#' @format GRanges
#' @examples
#' gr38[1:5,1:6]
"gr38"

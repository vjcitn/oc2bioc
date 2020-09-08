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
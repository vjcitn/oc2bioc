#' utility for working with multi-experiment results from ENCODE TFBS in opencravat
#' @param x a properly filtered data.frame with the relevant columns
#' @param tf character(1) symbol for TF as recorded by ENCODE
#' @note The issue is that many TF experiments are recorded in a single excel spreadsheet cell. 
#' A column gives the TFs studied with semicolon delimiter, and other columns give the
#' cell line, quality of result.   These columns should be coupled with relevant variant-level
#' metadata
#' @examples
#' \dontrun{
#' x = read.csv("IPF GWAS Variants - Variant.csv")
#' inds = grep("YY1", x$ENCODE.TFBS)
#' xyy = x[inds,]
#' dd = data.frame(t(apply(xyy[,19:23], 1, pullTF)))
#' hh = cbind(xyy, dd)
#' impactYY1 = hh[,c(1:10,43,54:58)]
#' write.csv(impactYY1, "encodeYY1.csv")
#' }
pullTF = function(x, tf="YY1") {
 y = lapply(x, strsplit, ";")
 ind = grep(tf, y[[1]][[1]])[1]
 sapply(y, function(x) x[[1]][ind])
 }


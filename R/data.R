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

#' annotation of mutations in OpenCRAVAT's crx format
#' @docType data
#' @format data.frame
#' @note 500 records retained
#' @examples
#' crx_demo[1:5,1:6]
"crx_demo"

#' Sequence Ontology abbreviations and tags used with OpenCRAVAT
#' @docType data
#' @format data.frame
#' @examples
#' head(SO_map)
"SO_map"

#' Sequence Ontology ingested from OBO
#' @docType data
#' @format ontologyIndex instance
#' @note Ingested from so.obo at https://github.com/The-Sequence-Ontology/SO-Ontologies/Ontology_Files
#' @examples
#' head(names(SO_onto))
#' head(SO_onto$name)
"SO_onto"

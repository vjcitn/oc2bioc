#' convert a GRanges to a data.frame that can be written out as a POSTable file for OpenCRAVAT
#' @importFrom rtracklayer liftOver
#' @importFrom GenomicRanges width start
#' @import GenomeInfoDb
#' @param gr GRanges instance
#' @param chain rtracklayer chain for liftOver
#' @param genome character(1) genome tag
#' @examples
#' tf = tempfile()
#' mtb = make_oc_POSTable(gr38)
#' write.table(mtb, file=tf, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
#' read.delim(tf, sep="\t", nrow=5, h=FALSE)
#' @export
make_oc_POSTable = function(gr, chain=NULL, genome=NULL) {
 stopifnot(inherits(gr, "GRanges"))
 if (!is.null(chain)) {
  gr = unlist(rtracklayer::liftOver(gr, chain))
  if (!is.null(genome)) genome(gr) = genome
  }
 gr = gr[width(gr)==1]
 v = paste0("v_", seq_len(length(gr)))
 gg = unique(genome(gr))[1]
 data.frame(chr=as.character(seqnames(gr)), pos=GenomicRanges::start(gr), ref=gr$Reference_Allele,
    alt=gr$Tumor_Seq_Allele2, samp=names(gr), var=v)
}

#' Illustration of variant vizualization with TnT
#' @import curatedTCGAData
#' @import TnT
#' @import GenomeInfoDb
#' @import GenomicRanges
#' @import MultiAssayExperiment
#' @import ensembldb
#' @import RaggedExperiment
#' @import Biobase
#' @import EnsDb.Hsapiens.v75
#' @importFrom grDevices rainbow
#' @param mae instance of MultiAssayExperiment, defaults to NULL
#' @param sitecode character(1) TCGA code, defaults to "ACC"
#' @param chr character(1) chromosome for variant request, defaults to "11"
#' @param viewstart numeric(1) leftmost basepair to view, defaults to 6e7
#' @param viewend numeric(1) rightmost basepair to view, defaults to 7e7 
#' @param verbose logical(1) passed to curatedTCGAData::curatedTCGAData
#' @note Mutation component is extracted from MultiAssayExperiment, and
#' rowRanges method used to obtain addresses of recorded somatic mutations.
#' Non-SNV mutations are dropped.
#' @return TnT::TnTGenome
#' @export
TnTdemo = function(mae=NULL, sitecode="ACC", chr="11", viewstart=6e7,
        viewend=7e7, verbose=TRUE) {
if (is.null(mae)) {
mae = curatedTCGAData::curatedTCGAData(sitecode, "Mutation", dry.run=FALSE,
    verbose=verbose) 
}
ra = MultiAssayExperiment::experiments(mae)[[1]]
rr11 = rowRanges(ra)[which(as.character(seqnames(rowRanges(ra)))==chr)]
genome(rr11) = "GRCh37"
#
# obtain exons
#
ex = ensembldb::exons(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75)
#
# build lollipop track with variants
#
nonco = setdiff(rr11,ex)
nonco$value=1
nonco = nonco[width(nonco)==1]
co = intersect(rr11, ex)
co$value = 1
co = co[width(co)==1]
pt = PinTrack(nonco, color="darkgreen")
copt = PinTrack(co, color="red")
#
# build gene track
#
gene = ensembldb::genes(EnsDb.Hsapiens.v75)
ensGeneTrack <- TnT::FeatureTrack(gene, 
    tooltip = as.data.frame(gene),
    names = paste(gene$symbol, " (", 
      gene$gene_biotype, ")", sep = ""),
    color = TnT::mapcol(gene$gene_biotype, palette.fun = grDevices::rainbow))
#
# viz
#
mytt = GRanges(chr, IRanges(viewstart,viewend))
TnT::TnTGenome(list(pt,copt,ensGeneTrack), view.range=mytt)
}

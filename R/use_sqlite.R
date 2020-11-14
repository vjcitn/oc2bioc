#' transform internal _header tables to char vectors
#' @importFrom jsonlite fromJSON
#' @importFrom magrittr %>%
#' @rawNamespace import("dplyr", except=c("select", "intersect"))
#' @import ggplot2
#' @param con DBI connection
#' @param tblname character(1) table name
#' @return vector of character
header_tx = function (con, tblname) 
{
    dd <- ((con %>% dplyr::tbl(tblname) %>% dplyr::select(col_def)) %>% 
        as.data.frame)[[1]]
    za <- sapply(dd, jsonlite::fromJSON)
    as.character(za[3, ])
}

#' get SQLite content for a table
#' @param con DBI connection
#' @param tablename character(1) table name
#' @param as.data.frame logical(1)
#' @param use_header logical(1) if TRUE look in database for table that elaborates table colnames
get_oc_sqlite_content = function(con, tablename, as.data.frame=TRUE, use_header=TRUE) {
 ans = con %>% dplyr::tbl(tablename)
 if (!as.data.frame) return(ans)
 ans = as.data.frame(ans)
 if (use_header) {
   hd = header_tx(con, paste(tablename, "_header", sep=""))
   names(ans) = hd
   }
 ans
}

#' get an optionally thinned table
#' @param con an RSQLite connection
#' @param tablename character(1)
#' @param thincols character() open-cravat columns
#' @param use_header logical(1) 
get_oc_tab = function(con, tablename="variant", thincols = c("All Mappings", "Samples", "Gene Targets"),
    use_header=TRUE) {
 vt = get_oc_sqlite_content(con, tablename, as.data.frame=TRUE, use_header=use_header)
 nv = names(vt) # column names
 inds = na.omit(match(thincols, nv))
 if (length(inds)>0) vt = vt[,-inds]
 vt
}

#' get the SO frequency
#' @param con an RSQLite connection
get_SO_freq = function(con) {
 tb = get_oc_tab(con, tablename="variant")
 ans = tb %>% dplyr::select(`Sequence Ontology`) %>% group_by(`Sequence Ontology`) %>% summarise(n=n())
 frq = ans$n
 code = ans$`Sequence Ontology`
 somap = SO_map$term
 names(somap) = SO_map$abbr
 reco = somap[code]
 data.frame(SO=as.character(reco), freq=frq)
}

#' make a sequence ontology barplot over all variants
#' @param con SQLite connection to openCRAVAT .sqlite
#' @return a ggplot instance
#' @export
barplot_sequence_ontology = function(con) {
 vt = get_oc_sqlite_content(con, "variant")
 somap = SO_map$term
 names(somap) = SO_map$abbr
 vt$seq_onto = somap[vt$`Sequence Ontology`]
 vt$chrom = ordered(vt$Chrom, levels=paste("chr", c(1:22, "X", "Y"), sep=""))
# ggplot2::ggplot(vt, aes(y=chr2, fill=seq_onto)) + ggplot2::geom_bar() 
 kk = vt %>% dplyr::select(chrom, seq_onto) %>% dplyr::group_by(chrom, seq_onto) %>%
   dplyr::summarise(n=dplyr::n()) 
 ggplot2::ggplot(kk, aes(x=chrom, fill=seq_onto)) + geom_col(aes(y=n)) + coord_flip()
}

#' make a gene ontology barplot over all variants
#' @param con SQLite connection to openCRAVAT .sqlite
#' @param kpev character() evidence codes to retain, defaulting to EXP, TAS, IPI, IDA, IMP, IGI
#' @param ont character() ontology codes to retain, defaulting to MF and BP
#' @param ncat numeric(1) maximum number of codes to report
#' @return a ggplot instance
#' @export
barplot_gene_ontology = function(con, kpev=c("EXP", "TAS",
     "IPI", "IDA", "IMP", "IGI"), ont=c("MF", "BP"), ncat=15) {
 gg = get_GO_table(con, kpev=kpev, ont=ont)
 dd = tail(sort(table(gg$term)), ncat)
 gg = gg[gg$term %in% names(dd),]
 ggg = gg %>% dplyr::group_by(term, EVIDENCE) %>% summarise(n=dplyr::n())
 ggg$n = pmin(ggg$n, 1000)
 ggplot2::ggplot(ggg, aes(x=term, y=n, fill=EVIDENCE)) + ggplot2::geom_col() + coord_flip()
}

#' @importFrom AnnotationDbi select mapIds
#' @import GO.db
#' @import org.Hs.eg.db
get_GO_table = function(con, kpev=c("EXP", "TAS", "IPI", "IDA", "IMP", "IGI"), ont="MF") {
 vt = get_oc_sqlite_content(con, "variant")
 gg = unique(as.character(na.omit(vt$Gene)))
 tab = AnnotationDbi::select(org.Hs.eg.db, keys=gg, keytype="SYMBOL", columns=c("GO", "ONTOLOGY"))
 tab = tab[which(tab$EVIDENCE %in% kpev),]
 tab = tab[which(tab$ONTOLOGY %in% ont),]
 tms = AnnotationDbi::mapIds(GO.db, keys=unique(tab$GO), keytype="GOID", column="TERM")
 tab$term = tms[tab$GO]
 tab
}


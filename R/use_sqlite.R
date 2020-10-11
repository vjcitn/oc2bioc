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
get_oc_sqlite_content = function(con, tablename, as.data.frame=TRUE) {
 ans = con %>% dplyr::tbl(tablename)
 if (!as.data.frame) return(ans)
 hd = header_tx(con, paste(tablename, "_header", sep=""))
 ans = as.data.frame(ans)
 names(ans) = hd
 ans
}

#' get an optionally thinned table
get_oc_tab = function(con, tablename="variant", thincols = c("All Mappings", "Samples", "Gene Targets")) {
 vt = get_oc_sqlite_content(con, tablename, as.data.frame=TRUE)
 nv = names(vt) # column names
 inds = na.omit(match(thincols, nv))
 if (length(inds)>0) vt = vt[,-inds]
 vt
}

#' get the SO frequency
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
 vt$chr2 = ordered(vt$Chrom, levels=paste("chr", c(1:22, "X", "Y"), sep=""))
 ggplot2::ggplot(vt, aes(y=chr2, fill=seq_onto)) + ggplot2::geom_bar() 
}

#' make a gene ontology barplot over all variants
#' @param con SQLite connection to openCRAVAT .sqlite
#' @param kpev character() evidence codes to retain, defaulting to EXP and TAS
#' @param ont character() ontology codes to retain, defaulting to MF and BP
#' @return a ggplot instance
#' @export
barplot_gene_ontology = function(con, kpev=c("EXP", "TAS"), ont=c("MF", "BP"), ncat=15) {
 gg = get_GO_table(con, kpev=kpev, ont=ont)
 dd = tail(sort(table(gg$term)), ncat)
 gg = gg[gg$term %in% names(dd),]
 ggplot2::ggplot(gg, aes(y=term, fill=EVIDENCE)) + ggplot2::geom_bar() 
}

#' @importFrom AnnotationDbi select mapIds
#' @import GO.db
#' @import org.Hs.eg.db
get_GO_table = function(con, kpev=c("EXP", "TAS"), ont="MF") {
 vt = get_oc_sqlite_content(con, "variant")
 gg = unique(as.character(na.omit(vt$Gene)))
 tab = AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys=gg, keytype="SYMBOL", columns=c("GO", "ONTOLOGY"))
 tab = tab[which(tab$EVIDENCE %in% kpev),]
 tab = tab[which(tab$ONTOLOGY %in% ont),]
 tms = AnnotationDbi::mapIds(GO.db::GO.db, keys=unique(tab$GO), keytype="GOID", column="TERM")
 tab$term = tms[tab$GO]
 tab
}


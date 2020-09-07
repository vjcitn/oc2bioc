
loginOC = function(uname, passwd)
 httr::GET("https://run.opencravat.org/server/login", httr::authenticate(uname, passwd))

logoutOC = function()
 httr::GET("https://run.opencravat.org/server/logout")


#' single-variant query
#' @import httr
#' @param uname user name registered at run.opencravat.org; will use value of environment variable OCRAVAT_USERNAME if missing
#' @param passwd password for opencravat account; will use value of environment variable OCRAVAT_PASSWORD if missing
#' @param chr character(1) chromosome in chrn format
#' @param pos character(1) position, in hg38 coordinates
#' @param ref character(1) reference base
#' @param alt character(1) alternate base
#' @param annotators character() vector of annotator tags
#' @examples
#' queryOC # don't want to use auth
#' @export
queryOC = function(uname, passwd, chr="chr1", pos="12777320", ref="G", alt="T",
  annotators = c("clinvar", "dbsnp", "exac_gene", "go", "rvis")) {
  if (missing(uname)) uname = Sys.getenv("OCRAVAT_USERNAME")
  if (missing(passwd)) passwd = Sys.getenv("OCRAVAT_PASSWORD")
  if (nchar(uname)==0) stop("no username available")
  loginOC(uname=uname, passwd=passwd)
  on.exit(logoutOC())
  q = paste0("https://run.opencravat.org/submit/annotate?",
       "chrom=", chr, "&pos=", pos, "&ref_base=", ref, "&alt_base=", alt,
       "&annotators=", paste(annotators, collapse=","))
  httr::GET(q)
}



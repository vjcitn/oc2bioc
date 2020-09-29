#' build a runnable python program that posts a file of variants to OpenCRAVAT
#' @importFrom reticulate py_run_file
#' @param url 'submit/submit' OpenCRAVAT endpoint url
#' @param postfile path to file with variants in TSV format
#' @param annotators character() vector of annotator names, must be installed in server identified by `url`
#' @param reports character() vector of reporter names, must be available in server identified by `url`
#' @param assembly character(1) hg38 or hg37
#' @param note character(1) arbitrary string
#' @return a reticulate python.builtin.dict with component `r` for which `r$json()` will return
#' metadata about run; an attribute `serverurl` is added with the serverurl string,
#' omitting the `submit/submit...`
#' @note This is an extremely fragile and cumbersome approach.  It would be better
#' to compose using a JSON constructor, and to use httr::POST.  But the quoting requirements
#' and the use of open() in the request seem hard for such an approach.
#' @examples
#' if (interactive()) run_oc_req()  # then use reticulate run_py_*
#' @export
run_oc_req = function(url="http://0.0.0.0:8080/submit/submit",
   postfile="/tmp/abc.txt", annotators=c("clinvar", "segway_lung"),
   reports=c("text", "vcf"), assembly="hg38", note="test run") {
 md = function(x) dQuote(x, q=FALSE)
 ms = function(x) sQuote(x, q=FALSE)
 mmd = function(x) paste(dQuote(x, q=FALSE), collapse=", ")
 tf = tempfile()
 on.exit(rm(tf))

 tmpl = "import requests; r = requests.post(%s, files={'file_0':open(%s)}, data={'options': '{%s: [%s], %s: [%s], %s: %s, %s: %s}'})"
 cat(sprintf(tmpl, ms(url), ms(postfile), md("annotators"), mmd(annotators), md("reports"),
   mmd(reports), md("assembly"), md(assembly), md("note"), md(note)), "\n", file=tf)
 ans = reticulate::py_run_file(tf)
 attr(ans, "serverurl") = gsub("submit..*", "", url)
 ans
}

#' browse report from OpenCRAVAT request
#' @param x output of `run_oc_req`
#' @note need to check that run is complete
#' @export
browse_oc_result = function(x) {
  uhead = attr(x, "serverurl")
  id = x$r$json()$id
  browseURL(paste0(uhead, "result/index.html?job_id=", id))
#  http://0.0.0.0:8080/result/index.html?job_id=200924-075831
  }

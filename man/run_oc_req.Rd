% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_req.R
\name{run_oc_req}
\alias{run_oc_req}
\title{build a runnable python program that posts a file of variants to OpenCRAVAT}
\usage{
run_oc_req(
  url = "http://0.0.0.0:8080/submit/submit",
  postfile = "/tmp/abc.txt",
  annotators = c("clinvar", "segway_lung"),
  reports = c("text", "vcf"),
  assembly = "hg38",
  note = "test run"
)
}
\arguments{
\item{url}{'submit/submit' OpenCRAVAT endpoint url}

\item{postfile}{path to file with variants in TSV format}

\item{annotators}{character() vector of annotator names, must be installed in server identified by `url`}

\item{reports}{character() vector of reporter names, must be available in server identified by `url`}

\item{assembly}{character(1) hg38 or hg37}

\item{note}{character(1) arbitrary string}
}
\value{
a reticulate python.builtin.dict with component `r` for which `r$json()` will return
metadata about run; an attribute `serverurl` is added with the serverurl string,
omitting the `submit/submit...`
}
\description{
build a runnable python program that posts a file of variants to OpenCRAVAT
}
\note{
This is an extremely fragile and cumbersome approach.  It would be better
to compose using a JSON constructor, and to use httr::POST.  But the quoting requirements
and the use of open() in the request seem hard for such an approach.
}
\examples{
if (interactive()) run_oc_req()  # then use reticulate run_py_*
}

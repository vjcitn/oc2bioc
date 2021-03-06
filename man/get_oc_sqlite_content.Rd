% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/use_sqlite.R
\name{get_oc_sqlite_content}
\alias{get_oc_sqlite_content}
\title{get SQLite content for a table}
\usage{
get_oc_sqlite_content(con, tablename, as.data.frame = TRUE, use_header = TRUE)
}
\arguments{
\item{con}{DBI connection}

\item{tablename}{character(1) table name}

\item{as.data.frame}{logical(1)}

\item{use_header}{logical(1) if TRUE look in database for table that elaborates table colnames}
}
\description{
get SQLite content for a table
}

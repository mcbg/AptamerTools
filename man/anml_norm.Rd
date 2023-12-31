% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalization.R
\name{anml_norm}
\alias{anml_norm}
\title{ANML}
\usage{
anml_norm(ds, v, ref, ext, mad = "mad", iterations = 100, verbose = FALSE)
}
\arguments{
\item{ds}{data.table}

\item{v}{String of name of variable for values}

\item{ref}{String of name of variable for references}

\item{ext}{String giving extension given to scale factor variable name}

\item{mad}{String giving variable name of median absolute deviation of references}

\item{iteration}{Integer giving number of times to iterate}
}
\description{
ANML
}

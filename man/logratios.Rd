% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logratios.R
\name{logratios}
\alias{logratios}
\title{logratios}
\usage{
logratios(W, base, p = 0.05)
}
\arguments{
\item{W}{raw count data, with OTUs as columns}

\item{base}{base OTU value}

\item{p}{how much to purturb zero counts, defaults to 0.05}
}
\description{
This function transforms from count data to logratios
}

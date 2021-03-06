% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MCmat.R
\name{MCmat}
\alias{MCmat}
\title{MCmat}
\usage{
MCmat(Y, W, eY, N, Q, base, sigma, MCiter, stepsize = 1, poorman = FALSE)
}
\arguments{
\item{Y}{logratio matrix}

\item{W}{corresponding count matrix}

\item{eY}{current expected value of logratio matrix}

\item{N}{number of samples, or nrow of Y}

\item{Q}{number of OTUs minus the base, or ncol of Y}

\item{base}{OTU index used for base}

\item{sigma}{current estimate of sigma}

\item{MCiter}{number of MC samples to generate}

\item{stepsize}{variance used for MH samples, defaults to 1. Tweak to adjust acceptance ratio}

\item{poorman}{boolean of whether to just use simple diagonal inverse to calculate sigma inverse. Defaults to FALSE}
}
\description{
This function simulates MC step for an entire matrix. Should not need to be used by user directly. NOTE: currently written in parallel
}

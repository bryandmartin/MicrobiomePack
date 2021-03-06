% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/savePDF.R
\name{savePDF}
\alias{savePDF}
\title{savePDF}
\usage{
savePDF(filepath, img, latex = TRUE)
}
\arguments{
\item{filepath}{computer specific filepath for saving, defaults to current working directory}

\item{img}{command used to generate a plot}

\item{latex}{whether to use Computer Modern font, defaults to true. Note if there are issues may need font_import(), font_install('fontcm'), loadfonts()}
}
\description{
This function saves a pdf image to a given filepath, with LaTeX fonts
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/StoreFTRANSFACFile.R
\name{StoreFTRANSFACFile}
\alias{StoreFTRANSFACFile}
\title{Generation and storing a file of the standard TRANSFAC format}
\usage{
StoreFTRANSFACFile(TheObject)
}
\arguments{
\item{TheObject}{This argument is an object of the class which holds the information ready to be plotted.IDs, Scores and Froked_PWM are mednatory.}
}
\value{
This function stores a .txt file at working directory, and returns name of the file for more convenience.
}
\description{
This function generates a .txt file of the format TRANSFAC with slight modifications in positions column.
}

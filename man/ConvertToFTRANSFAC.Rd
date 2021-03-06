% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ConvertToFTRANSFAC.R
\name{ConvertToFTRANSFAC}
\alias{ConvertToFTRANSFAC}
\title{Generating proper matrix similar to TRANSFAC format of all matrices.}
\usage{
ConvertToFTRANSFAC(TheObject)
}
\arguments{
\item{TheObject}{This argument is an object of the class which holds the information ready to be plotted.}
}
\value{
This class receives a class Object which holds the plotting data, and updates it by adding the proper matrix of new format: FTRANSFAC.
}
\description{
This function generates a matrix of 5 column (Position,A,T, C, G) with redundant position numbers at Position column reflecting number of leafs and their PWMs.
}

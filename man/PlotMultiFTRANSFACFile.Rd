% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PlotMultiFTRANSFACFile.R
\name{PlotMultiFTRANSFACFile}
\alias{PlotMultiFTRANSFACFile}
\title{A function for storing .PDF of plots, by providing a .txt file of FPWMs concatination, in proper format.}
\usage{
PlotMultiFTRANSFACFile(File = "All.txt")
}
\arguments{
\item{File}{the directory of .txt file of multiple FPWMs merged in proper format.}
}
\value{
Stores number of PDF files regarding the number of FPWMs provided within the file.
}
\description{
This function reads an stored .txt file of multiple FTRANSFAC matrices and generates the associated plot for each set then stores the figure as a PDF file. Name of each files indicates from each line the infomration is being imported to result to given plot.
}

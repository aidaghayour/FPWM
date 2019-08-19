#' A function for storing .PDF of plots, by providing a .txt file of FPWMs concatination, in proper format.
#'
#' This function reads an stored .txt file of multiple FTRANSFAC matrices and generates the associated plot for each set then stores the figure as a PDF file. Name of each files indicates from each line the infomration is being imported to result to given plot.
#' @param File the directory of .txt file of multiple FPWMs merged in proper format.
#' @return Stores number of PDF files regarding the number of FPWMs provided within the file.
#' @export
PlotMultiFTRANSFACFile <- function(File="All.txt") {
Lines <- readLines(File)
EndLines <- grep(pattern =  "//",x = Lines)

skip=0
for (i in EndLines)
{
  file_name <- paste0(i, "_chunck.txt")
  content_all <- as.vector(reader::n.readLines(fn = File, skip = skip , n = (i-skip),header = FALSE))
  fileConn <- file(file_name)
  writeLines(content_all, fileConn)
  close(fileConn)
  
  pdf(paste0("starting from line ", skip , ".pdf"))
  X <- ReadFTRANSFACFile(File = file_name, Methylation = FALSE)
  FPWMPlotter(X, Methylation = FALSE)
  dev.off()
  file.remove(file_name)
  skip = i
}
}
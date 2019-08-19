#' Generation and storing a file of the standard TRANSFAC format
#'
#' This function generates a .txt file of the format TRANSFAC with slight modifications in positions column.
#' @param TheObject This argument is an object of the class which holds the information ready to be plotted.IDs, Scores and Froked_PWM are mednatory.
#' @return This function stores a .txt file at working directory, and returns name of the file for more convenience.
#' @export
StoreFTRANSFACFile <- function(TheObject)
{
  line1 <- paste0("Parent logo : ", TheObject@xid)
  line2 <- "XX"
  line3 <- paste0("Leaf logos : ", paste(TheObject@id, collapse = ','))
  line4 <- paste0("Overlapping scores : ", paste(TheObject@score, collapse = ','))
  line5 <- "XX"
  line6 <- paste0("Forking point : ", TheObject@sp)
  line7 <- paste0("PO\tA\tC\tG\tT")
  motif_matrix <- TheObject@forked
  line_matrix <- ""
  for (i in seq(1, nrow(motif_matrix), 1))
  {
    if(i < nrow(motif_matrix))
    {
      line_matrix <- paste0(line_matrix, motif_matrix[i,1], "\t", motif_matrix[i,2], "\t", motif_matrix[i,3], "\t",  motif_matrix[i,4], "\t", motif_matrix[i,5], "\n")
    }
    else
    {
      line_matrix <- paste0(line_matrix, motif_matrix[i,1], "\t", motif_matrix[i,2], "\t", motif_matrix[i,3], "\t", motif_matrix[i,4], "\t", motif_matrix[i,5])
    }
  }
  Total <- round((Reduce("+", TheObject@score)), digits=2)
  line_aftermatrix_1 <- "XX"
  line_aftermatrix_2 <- "CC program: Forked-PWM"
  line_aftermatrix_3 <- paste0("Source: intersectPeakMatrix of ",TheObject@xid, " and ", length(TheObject@matrix), " other matrices.")
  line_aftermatrix_4 <- paste0("Matirx built by ",Total, " % Of total peaks.")
  line_aftermatrix_5 <- "XX"
  line_aftermatrix_6 <- "//"
  file_name <- paste0(TheObject@xid, "___",length(TheObject@matrix) ,"-FTRANSFAC.txt")
  content_all <- c(line1, line2, line3, line4, line5, line6, line7,line_matrix, line_aftermatrix_1, line_aftermatrix_2, line_aftermatrix_3, line_aftermatrix_4,line_aftermatrix_5,line_aftermatrix_6)
  fileConn <- file(file_name)
  writeLines(content_all, fileConn)
  close(fileConn)
  return(file_name)
}

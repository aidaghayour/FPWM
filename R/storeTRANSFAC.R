#' A function for storing TRANSFAC files of the forked matrices.
#'
#' This function generates files of regular TRANSFAC format in order to further analysis and evaluation. Each file name holdes the name of Transfactor of interest, and the co-factor that is under analysis in the current matrix.
#' @param TheObject the input is the object of FPWM class that holds the raw matrices directly exported from TFregulomeR().
#' @export
storeTRANSFAC <- function(TheObject)
{
for (q in c(1:length(TheObject@id))){
line1 <- paste0("Parent logo : ", TheObject@xid)
line2 <- "XX"
line3 <- paste0("Overlapped with : ", paste(TheObject@id[[q]], collapse = ','))
line4 <- paste0("Overlapping score : ", paste(TheObject@score[[q]], collapse = ','))
line5 <- "XX"
line6 <- paste0("Forking point : ", TheObject@sp)
line7 <- paste0("PO\tA\tC\tG\tT")

a=sum(TheObject@matrix[[q]][1,])
b=sum(TheObject@parentmatrix[1,])
motif_matrix <- rbind(round(Parent/(b/a)),TheObject@matrix[[q]][(sp+1):nrow(TheObject@matrix[[q]]),])


line_matrix <- ""
for (i in seq(1, nrow(motif_matrix), 1))
{
  if(i < nrow(motif_matrix))
  {
    line_matrix <- paste0(line_matrix, i, "\t", motif_matrix[i,1], "\t", motif_matrix[i,2], "\t", motif_matrix[i,3], "\t",  motif_matrix[i,4], "\n")
  }
  else
  {
    line_matrix <- paste0(line_matrix, i, "\t", motif_matrix[i,1], "\t", motif_matrix[i,2], "\t", motif_matrix[i,3], "\t", motif_matrix[i,4], "\t")
  }
}
Total <- round(TheObject@score[[q]], digits=2)
line_aftermatrix_1 <- "XX"
line_aftermatrix_2 <- "CC program: Forked-PWM"
line_aftermatrix_3 <- paste0("Source: intersectPeakMatrix of ",TheObject@xid, " and ",TheObject@id[[q]] )
line_aftermatrix_4 <- paste0("Matirx built by ",Total, " % Of total peaks.")
line_aftermatrix_5 <- "XX"
line_aftermatrix_6 <- "//"
file_name <- paste0(TheObject@xid, " forked to ",TheObject@id[[q]], "-TRANSFAC.txt")
content_all <- c(line1, line2, line3, line4, line5, line6, line7,line_matrix, line_aftermatrix_1, line_aftermatrix_2, line_aftermatrix_3, line_aftermatrix_4,line_aftermatrix_5,line_aftermatrix_6)
fileConn <- file(file_name)
writeLines(content_all, fileConn)
close(fileConn)
retrun(file_name)
}}
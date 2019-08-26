#' A function to merge Beta Score matrices to generate a single matrix.
#'
#' This function takes the class object and creates a merge of exclusive Beta Score Matrices
#' by calculating the elementwise weighted average of them; up to the forking position. Weight of each matrix, is the overlapping percentage of intersectPeakmatrix.
#' @param TheObject is an object of S4 class that holds original matrices exported from the package TFregulomeR().
#' @param sp is the forking position. User can define up to which position it is required to merge two matrices using this argument.
#' @examples This function is called within ClassAssignment() function.
#' @return This function receives a class object, and returns an updated class object.
#' @export
BetaAdder <- function( TheObject, sp)
{
  message("\n\nMethylation Score Matrices are adding up....\n\n")
  X <-  Map('*',TheObject@betalevel,TheObject@score)
  S <- X[[1]][,1:sp]
  W <- TheObject@score[[1]]
for ( i in c(2:length(TheObject@id))){
      S <- X[[i]][,1:sp]+S
      W <- TheObject@score[[i]] + W
      }
  TheObject@parentbeta <- as.matrix(round(S/W))
  
  return(TheObject)
  message("\n\n....Methylation Score Matrices are added and assigned as parent node's matrix....\n\n")
}

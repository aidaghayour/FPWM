#' A function to merge motif matrices to generate one matrix as parent node.
#' 
#' This function takes the object and creates a merge of all matrices 
#' by claculating the elementwise addition of them, up to a user specified position (Forking Position).
#' @param TheObject is a object of S4 class that holds original matrices exported from the package TFregulomeR().
#' @param sp is the forking position. User can define up to which position it is required to merge matrices using this argument.
#' @examples This function is called whithin ClassAssignment() function. 
#' @return This function recieves a class object, and returns an updated class object by adding merged matrix to parentmatrix slot
#' @export
MatrixAdder <- function( TheObject, sp)
{ S <- TheObject@matrix[[1]][1:sp,]
  for ( i in c(2:length(TheObject@matrix))){
  S <- TheObject@matrix[[i]][1:sp,]+S}
  TheObject@parentmatrix <- S
  
  return(TheObject)
}
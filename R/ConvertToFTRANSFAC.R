#' Generating proper matrix similar to TRANSFAC format of all matrices.
#'
#' This function generates a matrix of 5 column (Position,A,T, C, G) with redundant position numbers at Position column reflecting number of leafs and their PWMs.
#' @param TheObject This argument is an object of the class which holds the information ready to be plotted.
#' @return This class receives a class Object which holds the plotting data, and updates it by adding the proper matrix of new format: FTRANSFAC.
#' @export
ConvertToFTRANSFAC <- function(TheObject)
{ 
  message("\n\n... Matrices are converting to Forked-TRANSFAC format....\n\n")
  Cnumber = length(TheObject@matrix)
  RowNum = nrow(TheObject@matrix[[1]])
  sp = TheObject@sp
  R = rep(c((sp + 1):RowNum), times = Cnumber)
  Step = (RowNum - sp)
  DF <- data.frame(colnames(c("PO","A","C","G","T")))
  DF[1:sp,"PO"]=c(1:sp)
  DF[(sp+1):(length(R)+sp),"PO"] <- R
  DF[1:sp,2:5]=TheObject@parentmatrix
  c = 1
  for(i in seq(sp+1, dim.data.frame(x = DF)[1], Step)){
    DF[i:((Step+i)-1),2:5] <- TheObject@matrix[[c]][(sp+1):RowNum,]
    c <- c+1
  }
  TheObject@forked <- DF
  message("\n\n...FPWM is stored at @forked slot of the object ....\n\n")
  return(TheObject)

}

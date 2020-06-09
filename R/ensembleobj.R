#' A function for generating the object for TF in different cell lines
#'
#' This function browses all the cell lines of a given TF and augoments the provided motif with added cell lines to represent the impact of cell line on motif structure
#' @param sp Forking point for final plot
#' @param CelllinesNumb Maximum number of cell lines to be considered
#' @param tfname The name of target tf in strings
#' @param tfID the targeted motif using MethMotif IDs as the target cell line under study.
#' @export
Ensembles <- function( sp, tfname, tfID, CelllinesNumb)
{
  tmp <- ensemblesfunc(tfname = tfname, tfID = tfID)
  
  
  Matrices <- tmp$Matrices[1:CelllinesNumb]
  Betas <- tmp$Betas[1:CelllinesNumb]
  ID <- tmp$ID[1:CelllinesNumb]
  
  MultiClass <- new("FPWMClassObj")
  message("\n....Assigning data to class slots!...." )
  MultiClass <- updateFPWMClassObj(MultiClass, id = ID,
                                   matrix=Matrices,
                                   betalevel = Betas,
                                   sp = sp)
  
  message("\n....Merging matrices and converting them to Forked-TRANSFAC format!...." )
  MultiClass@parentmatrix <- Matrices[[1]][1:sp,]
  MultiClass@parentbeta <- Betas[[1]][,1:sp]
  MultiClass <- ConvertToFTRANSFAC(MultiClass)
  N=length(MultiClass@betalevel)
  for (i in c(1:N)) {
    X<-MultiClass@betalevel[[i]]
    MultiClass@betalevel[[i]] <- X[,(sp+1):ncol(MultiClass@betalevel[[i]])]
  }
  MultiClass@xid <- "peak_id_x"
  message("\n....Methylation Score matrices are modified!...." )
  
  MultiClass <- ModifyBetaFormat(MultiClass)
  message("\n....Returning class object!...." )
  MultiClass@score <- as.list(c(1:CelllinesNumb))
  return(MultiClass)
}
#' A function for generating object of TFBStools holding PFM of each fokred matrix.
#'
#' This function generates a TFBStools object of each matrix present in FPWM class object, and returns a list containing all the objects.
#' @param TheObject is the object of the FPWM class. It needs to contain the matrices, IDs and parent matrix.
#' @export
ToTFBSTools <- function(TheObject)
{
  ls <- list()

  for (q in c(1:length(TheObject@id)))
  {
    a = sum(TheObject@matrix[[q]][1, ])
    b = sum(TheObject@parentmatrix[1, ])
    motif_matrix <-
      rbind(round(TheObject@parentmatrix / (b / a)), TheObject@matrix[[q]][(TheObject@sp + 1):nrow(TheObject@matrix[[q]]), ])
    
    assign(
      paste0("PWM_", TheObject@id[[q]]),
      TFBSTools::PFMatrix(
        profileMatrix = t(motif_matrix),
        ID = toString(TheObject@id[[q]]),
        strand = "*",
        name = paste0("Forked from ", TheObject@xid)
      )
    )
    ls[[q]] <- get(eval(paste0("PWM_", TheObject@id[[q]])))
  }
  return(ls)
}
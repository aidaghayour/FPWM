#' A function togenerate a class object then assign proper data exported from TFregulomeR to its slots.
#'
#' This function assigns proper data to their associated slots of a S4 classe. This information is either provided by user, or exported from TFregulomeR's dataware using user specified data.
#' @param sp This argument, defines from which point on, the matrix needs to be forked, or in the other words, up to which point two exclusive matrices need to be aggregated.
#' @param peak_id_x This argument holds an id of TFBS compatible with TFregulomeR(). This is the target peak ID wich will be employed by IntersecPeakMatrix of TFregulomeR to extract desired data.
#' @param peak_id_y_list This argument is a list of TF ID's which will be intersected with Peak_id_x. 
#' @param height An argument which allows user to customize the height of final graph relative to screen.
#' @param width An argument which allows user to customize the width of final graph relative to screen.
#' @examples
#'
#' @return This component, returns a class object which holds all the neccessary information for other functuins.
#' @export
ObjectGenerator <- function( sp,peak_id_y_list,peak_id_x, height=2, width=3)
{
  # sp = 5
  # peak_id_y_list = list("MM1_HSA_K562_CEBPD","MM1_HSA_K562_CEBPD","MM1_HSA_K562_ATF4","MM1_HSA_K562_ATF4","MM1_HSA_K562_JUND", "MM1_HSA_K562_ATF7","MM1_HSA_K562_ATF1","MM1_HSA_K562_JUN")
  # peak_id_x = "MM1_HSA_K562_CEBPB"
  
  for (i in c(1:length(peak_id_y_list))){
    message("\n....Accessing TFregulomeR!...." )
    assign(paste0("Motif",i),TFregulomeR::intersectPeakMatrix(motif_type ="TRANSFAC",peak_id_x = peak_id_x, motif_only_for_id_x = T, peak_id_y = peak_id_y_list[i]))}
  
    
    TFregulomeDataovlaplist <- list()
    TFregulomeDataBetalist <- list()
    TFregulomeDatamatrixlist <- list()
    for (i in c(1:length(peak_id_y_list))){
      x = as.name(paste0("Motif",i))
      y <- eval(x)
      TFregulomeDataovlaplist[[i]] = y[peak_id_x,peak_id_y_list[i][[1]]][[1]]@overlap_percentage_x
      TFregulomeDataBetalist[[i]] = y[peak_id_x,peak_id_y_list[i][[1]]][[1]]@MethMotif_x@MMBetaScore
      TFregulomeDatamatrixlist[[i]] = y[peak_id_x,peak_id_y_list[i][[1]]][[1]]@MethMotif_x@MMmotif@motif_matrix
      }
    
    
  
  MultiClass <- new("FPWMClassObj")
  message("\n....Assigning data to class slots!...." )
  MultiClass <- updateFPWMClassObj(MultiClass, id = peak_id_y_list,
                                 matrix=TFregulomeDatamatrixlist,
                                 betalevel = TFregulomeDataBetalist,
                                 score= TFregulomeDataovlaplist,
                                 sp = sp)
  message("\n....Adding up PWMs!...." )
  MultiClass <- MatrixAdder( MultiClass, sp)
  message("\n....Adding up Methylation Score matrices!...." )
  MultiClass <- BetaAdder(MultiClass, sp)
  message("\n....Merging matrices and converting them to Forked-TRANSFAC format!...." )
  MultiClass <- ConvertToFTRANSFAC(MultiClass)
  
  for (i in c(1:length(MultiClass@betalevel))) {
    X<-MultiClass@betalevel[[i]]
    MultiClass@betalevel[[i]] <- X[,(sp+1):ncol(MultiClass@betalevel[[i]])]
  }
  MultiClass@xid <- peak_id_x
  message("\n....Methylation Score matrices are modified!...." )
  MultiClass <- ModifyBetaFormat(MultiClass)
  message("\n....Returning class object!...." )
return(MultiClass)
}

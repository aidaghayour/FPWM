#' A function for converting Beta Score matrices into proper data frames.
#'
#' This function receives the S4 class object and converts Betalevel matrices into data frames for better plotting and browsing purposes.
#' @param TheObject is an object of S4 class that holds original matrices exported from the package TFregulomeR().
#' @examples This function is called within ClassAssignment() function.
#' @return This function receives a class object, and returns an updated class object by modifying Beta Score Matrices.
#' @export
ModifyBetaFormat <- function(TheObject)
{
  BS1 <- TheObject@parentbeta
  BS1 <-
    cbind(c("beta score < 10%", "beta score in between", "beta score > 90%"),
          BS1)
  BS1 <- rbind(c("position", c(1:(ncol(BS1) - 1))), BS1)


  M1 <- matrix(nrow = ((ncol(BS1) - 1) * 3), ncol = 3)


  pos1 <- rep(t(BS1[1, 2:ncol(BS1)]), times = 3)
  M1[, 1] <- pos1

  pos1 <- t(BS1[2, 2:ncol(BS1)])
  M1[1:(ncol(BS1) - 1), 2] <- pos1
  pos1 <- t(BS1[3, 2:ncol(BS1)])
  M1[ncol(BS1):((ncol(BS1) - 1) * 2), 2] <- pos1
  pos1 <- t(BS1[4, 2:ncol(BS1)])
  M1[(((ncol(BS1) - 1) * 2) + 1):((ncol(BS1) - 1) * 3), 2] <- pos1

  pos1 <- rep(t(BS1[2, 1]), times = (ncol(BS1) - 1))
  M1[1:(ncol(BS1) - 1), 3] <- pos1
  pos1 <- rep(t(BS1[3, 1]), times = (ncol(BS1) - 1))
  M1[ncol(BS1):((ncol(BS1) - 1) * 2), 3] <- pos1
  pos1 <- rep(t(BS1[4, 1]), times = (ncol(BS1) - 1))
  M1[(((ncol(BS1) - 1) * 2) + 1):((ncol(BS1) - 1) * 3), 3] <- pos1

  TheObject@parentbeta <- M1
  message("\n....Parent node's Methylation Score Matrix is modified to right format!...." )
  #####################################
  
  for (i in c(1:length(TheObject@betalevel))) {
    
    
  
    BS2 <- as.matrix(TheObject@betalevel[[i]])
    BS2 <-
      cbind(c("beta score < 10%", "beta score in between", "beta score > 90%"),
            BS2)
    BS2 <- rbind(c("position", c(TheObject@sp + 1:(ncol(
      BS2
    ) - 1))), BS2)
    M2 <- matrix(nrow = ((ncol(BS2) - 1) * 3), ncol = 3)
  
  
    pos2 <- rep(t(BS2[1, 2:ncol(BS2)]), times = 3)
    M2[, 1] <- pos2
    
    pos2 <- t(BS2[2, 2:ncol(BS2)])
    M2[1:(ncol(BS2) - 1), 2] <- pos2
    pos2 <- t(BS2[3, 2:ncol(BS2)])
    M2[ncol(BS2):((ncol(BS2) - 1) * 2), 2] <- pos2
    pos2 <- t(BS2[4, 2:ncol(BS2)])
    M2[(((ncol(BS2) - 1) * 2) + 1):((ncol(BS2) - 1) * 3), 2] <- pos2
    
    pos2 <- rep(t(BS2[2, 1]), times = (ncol(BS2) - 1))
    M2[1:(ncol(BS2) - 1), 3] <- pos2
    pos2 <- rep(t(BS2[3, 1]), times = (ncol(BS2) - 1))
    M2[ncol(BS2):((ncol(BS2) - 1) * 2), 3] <- pos2
    pos2 <- rep(t(BS2[4, 1]), times = (ncol(BS2) - 1))
    M2[(((ncol(BS2) - 1) * 2) + 1):((ncol(BS2) - 1) * 3), 3] <- pos2
    
    TheObject@betalevel[[i]] <- M2}
  message("\n....Leaf nodes' Methylation Score Matrix is modified to right format!...." )
  
  return(TheObject)
}

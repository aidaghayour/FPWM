#' A function for generating barchart and seq-logo of co-factors of selected TF
#'
#' This function generates a barchart of co-binding Percentage for each co-factor of selected TF, along with seq-logo for each of co-factors.
#' @param NumberofTop Number of top co-factors with higher co-binding Percentage to be illustrated
#' @param highestscore Co-binding Percentage wich will be the minimum percentage of the shown co-factors.
#' @param cell A character string, which is the name of cell under study.
#' @param TF A character string which will be the Transcription Factor of interest.
#' @param Local A logical value, which will read a local .CSV file in case of TRUE. The file should contain two columns: scores, columnnames which are the co-binding percentages and IDs respectively.
#' @param path The path to .CSV file in case Local=TRUE.
#' @param Methylation Is a logic argument which indicates if user wants Methylation Score to bu plotted on top of sequence logos or not.
#' @export
Barandseqlogo <- function(NumberofTop, highestscore,cell,TF,Local = FALSE,path="",Methylation = FALSE)
{
  if( Local==TRUE & path=="" ) stop('You chose to work locally, but never provided path of local file!')
  x_peak_id <- paste0("MM1_HSA_",cell,"_",TF)
  if (Local == TRUE){
    message("\n\nYou chose to work with locally available .CSV file\n\n")
    df<-read.csv(path)
    scores <- as.numeric(df[,"scores"])
    columnnames <- as.vector(df[,"columnnames"])
    df <- data.frame(scores = scores, columnnames = columnnames)
  } else {
    
    message("\n\nYou chose to export data from TFregulomeR\n\n")
    
    k562_record <- TFregulomeR::TFBSBrowser(cell_tissue_name = cell)
    intersecmatrix_CEBPB_and_k562 <-
      TFregulomeR::intersectPeakMatrix(peak_id_x = x_peak_id, peak_id_y = k562_record$ID, motif_only_for_id_x = TRUE)
    results <-
      TFregulomeR::intersectPeakMatrixResult(
        intersectPeakMatrix = intersecmatrix_CEBPB_and_k562,
        return_intersection_matrix = TRUE,
        save_MethMotif_logo = FALSE
      )
    resultdata <- results$intersection_matrix
    resultdata <- sort( resultdata, decreasing = TRUE)
    scores <- as.numeric(resultdata[1,])
    columnnames <- colnames(resultdata)
    df <- data.frame(scores = scores, columnnames = columnnames)
    
    write.csv(df,'co-factors.csv')
    message("\n\nData file named 'co-factors.csv' is stored in current working directory! \n\n")
  }##########
  
  i <- 1
  while (1){
    
    if (scores[i] <= highestscore) {
      index <- (i-1)
      break()
    }
    i <- i+1
  }
  
  if(index >= NumberofTop){
    message("Number of co-factors with co-binding Percentage higher than provided, is more than or equal to selected number of Co-factors to show")
    message(paste0("...... Showing ", NumberofTop, " top ones..."))
  }else
  {
    message(paste0("...... Showing ", index , " of co-factors, among ", NumberofTop, " top ones."))
    NumberofTop <- index
  }
  df <- df[c(1:NumberofTop),]
  colss<-df[,"columnnames"]
  toberemoved <- paste0("MM1_HSA_",cell,"_")
  Pattern <- paste0("^.*?", toberemoved)
  xlabels <- gsub(pattern = toberemoved,x=colss, replacement = "")
  xlabels <- paste0("CEBPB + ", xlabels)
  p <-
    ggplot2::ggplot(data = df, ggplot2::aes(
      x = reorder(columnnames,scores),
      y = scores
    )) + ggplot2::geom_bar(stat = "identity")  + ggplot2::scale_x_discrete(labels= rev(as.vector(xlabels))) + ggplot2::theme(axis.text.y =
                                                                                                                               ggplot2::element_text(color = "black", size=8, angle=0, vjust=0.5, hjust=.95) ) + ggplot2::theme(legend.position = "none") +
    ggplot2::ylab("co-binding (Percentage)") + ggplot2::theme(
      axis.title.y  = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank()
    ) +  ggplot2::coord_flip()
  
  if (Methylation == TRUE) {
    T = (NumberofTop*2)
  } else
  {
    T = NumberofTop
  }
  
  
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(T, 2)))
  vplayout <- function(x, y) grid::viewport(layout.pos.row = x, layout.pos.col = y)
  
  
  
  print(p, vp = vplayout(1:T, 2))
  
  
  for (i in c(1:NumberofTop)) {
    X <-
      TFregulomeR::intersectPeakMatrix(
        peak_id_x = x_peak_id,
        motif_only_for_id_x = TRUE,
        peak_id_y = columnnames[i],
        motif_only_for_id_y = TRUE
      )
    
    M <-
      X[x_peak_id, columnnames[i]][[1]]@MethMotif_x@MMmotif@motif_matrix
    assign(paste0("p",i),
           ggplot2::ggplot() + ggseqlogo::geom_logo(t(M)) + ggseqlogo::theme_logo() +
             ggplot2::theme(axis.text.x = ggplot2::element_blank(),axis.text.y = ggplot2::element_blank()) + ggplot2::ylab('') + ggplot2::scale_y_continuous(limits = c(0, 2)))
    
    
    if (Methylation == TRUE){
      barplot_color <- c("dodgerblue1", "darkorange1", "darkgreen")
      Betas = list()
      ####
      BM <- X[x_peak_id, columnnames[i]][[1]]@MethMotif_x@MMBetaScore
      #Format
      BM <-
        cbind(c("beta score < 10%", "beta score in between", "beta score > 90%"),
              BM)
      BM <- rbind(c("position", c(1:(ncol(BM) - 1))), BM)
      
      
      M1 <- matrix(nrow = ((ncol(BM) - 1) * 3), ncol = 3)
      
      
      pos1 <- rep(t(BM[1, 2:ncol(BM)]), times = 3)
      M1[, 1] <- pos1
      
      pos1 <- t(BM[2, 2:ncol(BM)])
      M1[1:(ncol(BM) - 1), 2] <- pos1
      pos1 <- t(BM[3, 2:ncol(BM)])
      M1[ncol(BM):((ncol(BM) - 1) * 2), 2] <- pos1
      pos1 <- t(BM[4, 2:ncol(BM)])
      M1[(((ncol(BM) - 1) * 2) + 1):((ncol(BM) - 1) * 3), 2] <- pos1
      
      pos1 <- rep(t(BM[2, 1]), times = (ncol(BM) - 1))
      M1[1:(ncol(BM) - 1), 3] <- pos1
      pos1 <- rep(t(BM[3, 1]), times = (ncol(BM) - 1))
      M1[ncol(BM):((ncol(BM) - 1) * 2), 3] <- pos1
      pos1 <- rep(t(BM[4, 1]), times = (ncol(BM) - 1))
      M1[(((ncol(BM) - 1) * 2) + 1):((ncol(BM) - 1) * 3), 3] <- pos1
      
      #Formet-End
      Betas[[i]] <- data.frame(M1)
      #Modif
      M1 <- data.frame(M1)
      M1 <- M1[order(as.numeric(as.character(M1$X1))), ]
      totals2 <-
        as.vector(by(as.numeric(as.character(M1$X2)),
                     as.numeric(as.character(M1$X1)),
                     sum))
      assign(paste0("labels",i),
             unlist(lapply(as.character(totals2), function(x)
               c(rep("", 2), x))))
      
      assign(paste0("posi",i),rep(totals2 + 1, each = 3))
      
      M1 <- data.frame(M1, pos = get(paste0("posi",i)), labels = get(paste0("labels",i)))
      
      M1$X1 <- as.character(M1$X1)
      M1$X1 <- factor(M1$X1, levels = unique(M1$X1))
      assign(paste0("N",i), M1)
      #Modif-End
      
    }}
  
  
  #Lim
  if (Methylation == TRUE){
    lim = c()
    for (i in c(1:NumberofTop)){
      lim[i] <- max(as.numeric(as.character(get(paste0("N", i))$labels)), na.rm = TRUE)
    }
    MaxLim <- round(max(lim) * 1.6)
    #Lim-End
    
    #Bps[i]
    for (q in c(1:NumberofTop)) {
      assign(paste0("BSp",q),
             
             ggplot2::ggplot(get(paste0("N", q)), ggplot2::aes(
               x = X1, y = X2, fill = X3
             ))
             + ggplot2::geom_bar(stat = "identity", position =
                                   "stack")
             + ggplot2::geom_text(ggplot2::aes(y = get(paste0("posi",q)), label = get(paste0("labels",q))), vjust = 0, size=3)
             + ggplot2::theme(
               legend.position = "none",
               axis.title.y = ggplot2::element_blank(),
               axis.title.x = ggplot2::element_blank(),
               axis.text.y = ggplot2::element_blank(),
               axis.ticks.y = ggplot2::element_blank(),
               axis.ticks.x = ggplot2::element_blank(),
               axis.text.x = ggplot2::element_blank(),
               plot.margin = ggplot2::margin(
                 t = 10,
                 r = 20,
                 b = 0,
                 l = 19,
                 unit = "pt"
               ),
               panel.grid.major = ggplot2::element_blank(),
               panel.grid.minor = ggplot2::element_blank(),
               panel.background = ggplot2::element_blank()
             ) + ggplot2::scale_x_discrete(breaks = NULL) + ggplot2::scale_y_discrete(breaks =
                                                                                        NULL, limits = as.character(0:MaxLim)) + ggplot2::scale_fill_manual(values = barplot_color)
      )
      
    }
    #Bps[i]-End
    
    
  }
  
  if(Methylation == FALSE){
    
    for (q in c(1:NumberofTop))
    {
      print(get(paste0("p",q)), vp = vplayout(q, 1))
      
    }
    
  }else{
    
    q=1
    for(j in seq(from =2, to=T,by=2)){
      print(get(paste0("p",q)), vp = vplayout((j), 1))
      print(get(paste0("BSp",q)), vp = vplayout(j-1, 1))
      q=q+1}
  }
  
}

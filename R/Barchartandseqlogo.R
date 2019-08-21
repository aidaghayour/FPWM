#' A function for generating barchart and seq-logo of co-factors of selected TF
#'
#' This function generates a barchart of co-binding Percentage for each co-factor of selected TF, along with seq-logo for each of co-factors.
#' @param NumberofTop Number of top co-factors with higher co-binding Percentage to be illustrated
#' @param highestscore Co-binding Percentage wich will be the minimum percentage of the shown co-factors.
#' @param cell A character string, which is the name of cell under study.
#' @param TF A character string which will be the Transcription Factor of interest.
#' @param Local A logical value, which will read a local .CSV file in case of TRUE. The file should contain two columns: scores, columnnames which are the co-binding percentages and IDs respectively.
#' @param path The path to .CSV file in case Local=TRUE.
#' @export
Barandseqlogo <- function(NumberofTop, highestscore,cell,TF,Local = FALSE,path="")
{
  x_peak_id <- paste0("MM1_HSA_",cell,"_",TF)
if (Local == TRUE){
  df<-read.csv(path)
  scores <- as.numeric(df[,"scores"])
  columnnames <- as.vector(df[,"columnnames"])
  df <- data.frame(scores = scores, columnnames = columnnames)
} else {



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
p <-
  ggplot2::ggplot(data = df, ggplot2::aes(
    x = reorder(columnnames,scores),
    y = scores
  )) + ggplot2::geom_bar(stat = "identity")  + ggplot2::scale_x_discrete(labels= rev(as.vector(xlabels))) + ggplot2::theme(axis.text.y =
                                                               ggplot2::element_text(color = "black", size=8, angle=0, vjust=0.5, hjust=.95) ) + ggplot2::theme(legend.position = "none") + ggplot2::xlab("TFregulomeR IDs of Co-factors") +
  ggplot2::ylab("co-binding (Percentage)") + ggplot2::theme(
    axis.ticks.y = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank()
  ) +  ggplot2::coord_flip()


grid::pushViewport(grid::viewport(layout = grid::grid.layout(NumberofTop, 2)))
vplayout <- function(x, y) grid::viewport(layout.pos.row = x, layout.pos.col = y)

print(p, vp = vplayout(1:NumberofTop, 2))
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
  p1 <-
    ggplot2::ggplot() + ggseqlogo::geom_logo(t(M)) + ggseqlogo::theme_logo() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),axis.text.y = ggplot2::element_blank()) + ggplot2::ylab('') + ggplot2::scale_y_continuous(limits = c(0, 2))
  
  print(p1, vp = vplayout(i, 1))
}
}

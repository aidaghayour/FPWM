#' A function for exporting the motif matrix and augomenting it with additional cell lines
#'
#' This function 
#' @param tfname The name of target tf in strings
#' @param the targeted motif using MethMotif IDs as the target cell line under study.
#' @export
ensemblesfunc <- function(tfname = "JUN",tfID = "MM1_HSA_K562_JUN"){
  
JUN_record <- dataBrowser(tf = tfname, species = "human")
K562_JUN_all_peaks <- loadPeaks(id = tfID,
                                includeMotifOnly = TRUE)



N <- nrow(JUN_record)
Matrices <- list()
Betas <- list()
Score <- list()
ID <- as.list(JUN_record$ID)

for (i in seq(1,N,1)){
  cell_i <- JUN_record$cell_tissue_name[i]
  id_i <- JUN_record$ID[i]
  common_peak_i <- commonPeaks(target_peak_id = tfID,
                               motif_only_for_target_peak = TRUE,
                               compared_peak_id = id_i,
                               motif_only_for_compared_peak = TRUE)
  common_peak_res_i <- commonPeakResult(commonPeaks = common_peak_i,
                                        return_common_peak_sites = TRUE)
  tmp <- paste0(eval(tfID),"_common_peaks")
  common_peak_regions_i <- eval(parse(text = paste0("common_peak_res_i$common_peak_list$",tmp)))
  K562_JUN_all_peaks[,cell_i] <- 0
  K562_JUN_all_peaks[(K562_JUN_all_peaks$id %in% common_peak_regions_i$id),cell_i] <- 1
 }

K562_JUN_all_peaks$sum <- apply(K562_JUN_all_peaks[,JUN_record$cell_tissue_name],1,sum)

# Peak numbers, MethMotif logs and read enrichments for N K562 JUN sub-ensembles for Supplementary Figure 2 and Figure 1A
sub_ensemble_peak_num <- c()
sub_ensemble_read_score <- as.data.frame(matrix(nrow = nrow(K562_JUN_all_peaks),
                                                ncol = 2))
colnames(sub_ensemble_read_score) <- c("cell_type_num", "read_fold_change")
sub_ensemble_read_score$cell_type_num <- K562_JUN_all_peaks$sum
for (i in seq(1,N,1)){
  peak_subset_i <- K562_JUN_all_peaks[which(K562_JUN_all_peaks$sum==i),]
  sub_ensemble_peak_num <- c(sub_ensemble_peak_num, nrow(peak_subset_i))
  sub_ensemble_read_score[which(sub_ensemble_read_score$cell_type_num==i),'read_fold_change'] <-
    peak_subset_i$tag_fold_change
  
  # MethMotif logo
  common_peak_i <- commonPeaks(user_target_peak_list = list(peak_subset_i[,1:5]),
                               user_target_peak_id = tfID,
                               compared_peak_id = tfID,
                               motif_only_for_compared_peak = TRUE)
  common_peak_res_i <- commonPeakResult(commonPeaks = common_peak_i,
                                        save_MethMotif_logo = TRUE)
  Matrices[[i]] <- common_peak_i[tfID,"common"][[1]]@MethMotif@MMmotif@motif_matrix
  Betas[[i]] <- common_peak_i[tfID,"common"][[1]]@MethMotif@MMBetaScore
  #Score[[i]] <- common_peak_i[tfID,"common"][[1]]@common_percentage
  file.rename(paste0(tfID,"_common_peaks-logo-entropy.pdf"),
              paste0(tfID,"_in_peaks_shared_by_",i,"_cell_types.pdf"))
}
rt <- list()
rt[["Matrices"]] <- Matrices
rt[["Betas"]] <- Betas
rt[["ID"]] <- ID
#rt[["Score"]] <- Score
return(rt)
}

#' Generating and storing a .txt file named "All.txt" which contains multiple FPWMs, contatinated together, respecting TRANSFAC format.
#'
#' This function generates a .txt file which holds number of the data structures needed for one set of plotting in FTRANSFAC format.
#' @param List_sp List of forking position numbers for each one of FPWM.
#' @param Listof_peak_id_y_list A list of lists. Each list within this list, is a set of IDs which are going to form one FPWM plot.
#' @param List_peak_id_x A list of IDs.  The ID in List_peak_id_x[i] will be employed to form multiple IntersecPeakMatrices with all the IDs exisiting in Listof_peak_id_y_list[i].
#' @return This function stores a .txt file at working directory, and returns name of the file for more convenience.
#' @export
StoreMultiTRANSFACFile <- function( List_sp,Listof_peak_id_y_list, List_peak_id_x){

# List_peak_id_x <- list("MM1_HSA_K562_CEBPB", "MM1_HSA_K562_CEBPD")
# List_peak_id_y2 <- list("MM1_HSA_K562_CEBPD","MM1_HSA_K562_CEBPD","MM1_HSA_K562_CEBPD","MM1_HSA_K562_CEBPD")
# List_peak_id_y1 <- list("MM1_HSA_K562_ATF4","MM1_HSA_K562_ATF4","MM1_HSA_K562_ATF4","MM1_HSA_K562_ATF4")
# Listof_peak_id_y_list <- list(List_peak_id_y1,List_peak_id_y2)
# List_sp <- list(5,6,4,5)
dir.create(path = "./All")
setwd("./All")

for (i in c(1:length(List_peak_id_x)))
{
  Object <- ObjectGenerator(peak_id_x= List_peak_id_x[[i]][1], peak_id_y_list = Listof_peak_id_y_list[[i]], sp= List_sp[[i]][1])
  StoreFTRANSFACFile(Object)
}

file_list <- list.files()
file.create("All.txt")
All <- readLines("All.txt")
for (file in file_list){
  x <- readLines(file)
  All <- append(x,All)
}
write(All, "../All.txt")
setwd("..")
unlink("All", recursive=TRUE)
}
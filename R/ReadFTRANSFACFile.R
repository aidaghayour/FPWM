#' A function for generating a class object from a local file in proper format
#'
#' This function reads an stored .txt file of FTRANSFAC format and constructs a class object from it. As default, the returned class Object does not contain Methylation Score matrices. If needed, files exported from TFregulomeR() with the same name and format should be provided before setting MEthylation==TRUE.
#' @param File the directory of .txt file
#' @param Methylation a logical argument which indicates if Methylation Score files are provided and needed to be included in Object or not.
#' @return A class object for plotting. The Methylation Score matrices can be optionally ommitted or not.
#' @export
ReadFTRANSFACFile <- function( File = "MM1_HSA_K562_CEBPB___4-FTRANSFAC.txt", Methylation = FALSE)
{
  Lines <- readLines(File)
  Forked_Matrix_Start_point <- grep(pattern =  "PO",x = Lines)
  Forked_Matrix_end_point <- (Forked_Matrix_Start_point - 2) + min(grep(pattern =  "XX",x = Lines[Forked_Matrix_Start_point:length(Lines)]))
  
  NumberOfRows <- (Forked_Matrix_end_point - Forked_Matrix_Start_point)
  SkipRows <- (Forked_Matrix_Start_point - 1 )
  
  ForkedMatrix <- read.table(sep = "\t", file = File, nrows = NumberOfRows, skip =SkipRows, header = TRUE )
  #####
  
  ParentID <- word(Lines[grep(pattern =  "Parent",x = Lines)],-1)
  ParentID <-gsub(" ", "",ParentID, fixed = TRUE)
  LeafIDs <- as.list(strsplit(sub(".*:","",Lines[grep(pattern =  "Leaf",x = Lines, ignore.case=TRUE)]), ",")[[1]])
  LeafIDs <-as.list(gsub(" ", "",LeafIDs, fixed = TRUE))
  Overlappscore <- as.list(as.numeric(strsplit(sub(".*:","",Lines[grep(pattern =  "overlapp",x = Lines,ignore.case=TRUE)]), ",")[[1]]))
  
  
  ForkingPoint<- as.numeric(word(Lines[grep(pattern =  "Forking",x = Lines,ignore.case = TRUE)],-1)[1])
  
  TheObject <- new("FPWMClassObj")
  TheObject <- updateFPWMClassObj(TheObject, xid =ParentID, id = LeafIDs,
                                  forked=ForkedMatrix,
                                  score= Overlappscore,
                                  sp = ForkingPoint)
  
  betalevel = list()
  j = 1
  
  if (Methylation == TRUE)
  {
    for (i in LeafIDs) {
      
    
      X <- read.table(list.files( path = '.' , pattern = paste0("*", i,"-", ".*Score*")), skip = 1)
      X[,1:4] = NULL
      colnames(X) <- c(1:dim(X)[2])
      betalevel[[j]] <- X
      j = j+1
    
  }
    TheObject@betalevel <- betalevel
    TheObject <- BetaAdder(TheObject, sp= TheObject@sp)
    for (i in c(1:length(TheObject@betalevel))) {
      X<-TheObject@betalevel[[i]]
      TheObject@betalevel[[i]] <- X[,(ForkingPoint+1):ncol(TheObject@betalevel[[i]])]
    }
    TheObject <- ModifyBetaFormat(TheObject)
  }
  
 
  
  
  
  
  return(TheObject)
  
}
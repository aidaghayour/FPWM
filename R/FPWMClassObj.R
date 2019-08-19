FPWMClassObj <- setClass(
  # Set the name for the class
  "FPWMClassObj",
  # define the slots
  slots = c(
    xid = "character",
    id = "list",
    matrix="list",
    betalevel = "list",
    score= "list",
    sp = "numeric",
    parentmatrix = "matrix",
    parentbeta = "matrix",
    forked = "data.frame"
    
  ),
  # Set the default values for the slots.
  prototype = list(
    xid = "",
    id = list(),
    matrix=list(),
    betalevel = list(),
    score= list(),
    sp = 5,
    forked = data.frame(),
    parentbeta = matrix(),
    parentmatrix = matrix() )
)


setGeneric(name = "updateFPWMClassObj",
           def = function(theObject,xid,
                          id, matrix, betalevel,score, sp,  forked, parentbeta,parentmatrix)
           {
             standardGeneric("updateFPWMClassObj")
           }
)


setMethod(f = "updateFPWMClassObj",
          signature(theObject = "FPWMClassObj"),
          definition = function(theObject,xid,
                                id, matrix, betalevel,score, sp, forked, parentbeta,parentmatrix)
          {
            if (missing(theObject))
            {
              stop("No FPWMClassObj object.Please use 'theObject = '")
            }
            
            if (!missing(xid))
            {
              if (!is.character(xid))
              {
                stop("'xid' should be character!")
              }
              theObject@xid <- xid
            }
            
            if (!missing(id))
            {
              if (!is.list(id))
              {
                stop("'id' should be list!")
              }
              theObject@id <- id
            }
            
            
              if (!missing(matrix))
              {
                if (!is.list(matrix))
                {
                  stop("'matrix' should be list!")
                }
                theObject@matrix <- matrix
              }
            
            
              if (!missing(betalevel))
              {
                if (!is.list(betalevel))
                {
                  stop("'betalevel' should be list!")
                }
                theObject@betalevel <- betalevel
              }
              
            
              if (!missing(score))
              {
                if (!is.list(score))
                {
                  stop("'score' should be list!")
                }
                theObject@score <- score
              }
            
            
            if (!missing(sp))
            {
              if (!is.numeric(sp))
              {
                stop("'sp' should be numeric!")
              }
              theObject@sp <- sp
            }
            
            
            
            if (!missing(forked))
            {
              if (!is.data.frame(forked))
              {
                stop("'forked' should be data frme!")
              }
              theObject@forked <- forked
            }
            
            
            if (!missing(parentbeta))
            {
              if (!is.matrix(parentbeta))
              {
                stop("'parentbeta' should be matrix!")
              }
              theObject@parentbeta <- parentbeta
            }
            
            
            if (!missing(parentmatrix))
            {
              if (!is.matrix(parentmatrix))
              {
                stop("'parentmatrix' should be matrix!")
              }
              theObject@parentmatrix <- parentmatrix
            }

            return(theObject)
          })
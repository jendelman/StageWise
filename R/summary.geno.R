#' Summarize information about genomic relationships
#'
#' Summarize information about genomic relationships
#' 
#' First column contains diagonal elements of the G matrix. For objects with dominance, column 2 is diagonal elements of the D matrix, and column 3 is the genomic inbreeding coefficient (scaled row-sum of the D matrix).
#' 
#' @param object object of \code{\link{class_geno}}
#'
#' @return matrix
#'
#' @include class_geno.R
#' @name summary.geno
#' @exportMethod summary
NULL

setGeneric("summary")
setMethod("summary",c(object="class_geno"),
          definition=function(object){
            if (class(object)=="class_genoD") {          
              x <- cbind(diagG=diag(object@G),diagD=diag(object@D),
                     Fg=object@Fg)
            } else {
              x <- matrix(diag(object@G),ncol=1)
              colnames(x) <- "diagG"
            }
            rownames(x) <- rownames(object@G)
            return(x)
          })
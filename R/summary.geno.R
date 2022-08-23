#' Summarize information about genomic relationships
#'
#' Summarize information about genomic relationships
#' 
#' First column contains diagonal elements of the G matrix. For objects with dominance, second column contains the genomic inbreeding coefficients (scaled row-sum of the D matrix).
#' 
#' @param object object of \code{\link{class_geno}}
#'
#' @return data frame
#'
#' @include class_geno.R
#' @name summary.geno
#' @exportMethod summary
NULL

setGeneric("summary")
setMethod("summary",c(object="class_geno"),
          definition=function(object){
            if (class(object)=="class_genoD") {          
              x <- data.frame(diagG=diag(object@G),Fg=object@Fg)
            } else {
              x <- data.frame(diagG=diag(object@G))
            }
            rownames(x) <- rownames(object@G)
            return(x)
          })
#' G matrix quantile
#'
#' G matrix quantile
#' 
#' Unlike the S3 method, \code{prob} must have length = 1
#' 
#' @param x object of \code{\link{class_geno}}
#' @param prob probability
#' @return data frame with the quantile of the G matrix coefficients for each id
#'
#' @include class_geno.R
#' @name quantile
#' @exportMethod quantile
NULL

setGeneric("quantile")
setMethod("quantile",c(x="class_geno"),
          definition=function(x,prob){
            stopifnot(length(prob)==1)
            Gtmp <- as.matrix(x@G)
            diag(Gtmp) <- NA
            out <- data.frame(id=rownames(x@G),
                       G=apply(Gtmp,1,quantile,probs=prob,na.rm=T))
            rownames(out) <- NULL
            return(out)
          })
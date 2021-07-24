#' G matrix quantile
#'
#' G matrix quantile
#' 
#' Unlike the generic method, this method only accepts a single probability, not a vector
#' 
#' @param x object of \code{\link{class_geno}}
#' @param probs probability
#' @return data frame with the quantile of the G matrix coefficients for each id
#'
#' @include class_geno.R
#' @name quantile
NULL

setMethod("quantile",c(x="class_geno"),
          definition=function(x,probs){
            stopifnot(length(probs)==1)
            Gtmp <- as.matrix(x@G)
            diag(Gtmp) <- NA
            out <- data.frame(id=rownames(x@G),
                       G=apply(Gtmp,1,quantile,probs=probs,na.rm=T))
            rownames(out) <- NULL
            return(out)
          })
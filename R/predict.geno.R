#' Predict additive (breeding) values from marker effects
#'
#' Predict additive (breeding) values from marker effects
#' 
#' Use the \code{\link{blup}} function with \code{what="marker"} to generate the data frame for \code{marker.effects}.  
#' 
#' @param object object of \code{\link{class_geno}}
#' @param marker.effects data frame with columns "marker","add.effect"
#' @return data frame with columns "id", "BV"
#'
#' @include class_geno.R
#' @name predict
#' @exportMethod predict
NULL

setGeneric("predict")
setMethod("predict",c(object="class_geno"),
          definition=function(object,marker.effects){
            stopifnot(marker.effects$marker==colnames(object@coeff))
            GEBV <- object@coeff %*% marker.effects$add.effect
            out <- data.frame(id=rownames(GEBV),BV=as.numeric(GEBV))
            rownames(out) <- NULL
            return(out)
          })
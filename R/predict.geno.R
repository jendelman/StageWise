#' Predict individual values from marker effects
#'
#' Predict individual values from marker effects
#' 
#' Use the \code{\link{blup}} function with what="AM" or "DM" to generate the data frame for \code{marker.effects}.  
#' 
#' @param object object of \code{\link{class_geno}}
#' @param marker.effects data frame with columns "marker" and "effect"
#' @return data frame with columns "id" and "value"
#'
#' @include class_geno.R
#' @name predict
#' @exportMethod predict
NULL

setGeneric("predict")
setMethod("predict",c(object="class_geno"),
          definition=function(object, marker.effects){
            stopifnot(marker.effects$marker %in% colnames(object@coeff))
            if (attr(marker.effects,"what")=="AM") {
              pred <- as.numeric(object@coeff[,marker.effects$marker] 
                                 %*% matrix(marker.effects$effect,ncol=1))
            } else {
              stopifnot(class(object)=="class_genoD")
              pred <- as.numeric(object@coeff.D[,marker.effects$marker] 
                                 %*% matrix(marker.effects$effect,ncol=1))
            }
            data.frame(id=rownames(object@coeff), value=pred)
          })
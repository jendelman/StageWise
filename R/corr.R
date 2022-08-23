#' Trait correlations
#' 
#' Trait correlations
#' 
#' Use either the argument \code{traits} or \code{effect}, not both. Using \code{traits} leads to a partitioning of the total correlation between those two traits, based on path analysis. Using \code{effect} displays the correlation between all traits for that effect. Use the \code{summary} command to see the names of the possible effects.
#' 
#' @param vars object of \code{\link{class_var}} from \code{\link{Stage2}}
#' @param traits pair of traits
#' @param effect name of effect
#' 
#' @return matrix
#' @export

corr <- function(vars,traits=NULL,effect=NULL) {
  stopifnot(is.null(traits)|is.null(effect))
  stopifnot(!is.null(traits)|!is.null(effect))
  stopifnot(ncol(vars@resid) > 1)
  
  ix <- which(!is.na(vars@vars[1,1,]))
  if (!is.null(traits)) {
    stopifnot(length(traits)==2)
    x <- summary(vars)
    stopifnot(traits %in% colnames(x$PVE))
    y <- lapply(as.list(rownames(x$PVE)),function(x){corr(vars,effect=x)})
    corrs <- sapply(y,function(z){z[traits[1],traits[2]]})
    z <- apply(cbind(x$PVE[,traits],corrs),1,function(v){v[3]*sqrt(v[1])*sqrt(v[2])})
    tmp <- matrix(c(z,sum(z)),ncol=1)
    colnames(tmp) <- "Corr"
    rownames(tmp) <- c(names(z),"Total")
    round(tmp,3)
  } else {
    stopifnot(effect %in% dimnames(vars@vars)[[3]][ix])
    round(cov2cor(vars@vars[,,effect]),3)
  }
}

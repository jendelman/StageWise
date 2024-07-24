#' Report dominance parameters
#'
#' Report dominance parameters
#' 
#' The dominance variance (Vd) and baseline heterosis (b) are quantified relative to additive variance (Va) and std. dev. (SDa), respectively. As of v1.11, the variances are scaled to the population (previously, it was just the variance components). For a multi-trait/loc model, \code{index.coeff} specifies the coefficients of the standardized true values (see also \code{\link{blup}}), with \code{gamma} indicating the relative weight of non-additive to additive genetic merit for the standardization (see also \code{\link{gain}}). 
#' 
#' @param params list returned by \code{\link{Stage2}}
#' @param index.coeff merit index coefficients
#' @param gamma contribution of non-additive values for genetic merit
#'
#' @return data frame with estimates 
#'
#' @export

dominance <- function(params, index.coeff=NULL, gamma=0) {
  if (!is.element("heterosis",names(params)))
    stop("Dominance model was not used")

  stopifnot("scale" %in% names(params))
  
  id <- grep("dominance",params$vc$name,fixed=T)
  ia <- grep("additive",params$vc$name,fixed=T)
  vc <- params$vc[,-1]
  rownames(vc) <- params$vc$name

  nhet <- nrow(params$heterosis) 
  if (nhet > 1) {
    stopifnot(!is.null(index.coeff))
    if (colnames(params$heterosis)[1]=="trait") {
      traits <- params$heterosis$trait
      Va <- f.cov.trait(vc[ia,],traits,us=(nhet>2))
      Vd <- f.cov.trait(vc[id,],traits,us=(nhet>2))
    } else {
      #multi-location
      traits <- params$heterosis$loc
      add <- f.cov.loc(vc[ia,], traits)
      Va <- add$cov.mat
      tmp <- matrix(vc[id,1],ncol=1)
      rownames(tmp) <- traits
      Vd <- tcrossprod(sqrt(tmp))
    }
    
    stopifnot(sort(names(index.coeff))==sort(traits))
    index.coeff <- index.coeff[match(traits,names(index.coeff))]
    trait.scale <- sqrt(diag(Va)+gamma^2*diag(Vd))
    index.coeff <- index.coeff/trait.scale
    
    Va.mean <- as.numeric(crossprod(index.coeff,Va%*%index.coeff))
    Vd.mean <- as.numeric(crossprod(index.coeff,Vd%*%index.coeff))
    b.mean <- sum(params$heterosis$estimate*index.coeff)
  } else {
    Vd.mean <- params$vc$estimate[id]
    Va.mean <- params$vc$estimate[ia]
    b.mean <- params$heterosis$estimate
  }
    
  Vd <- Vd.mean*params$scale$dom + b.mean^2*params$scale$heterosis
  Va <- Va.mean*params$scale$add

  ans <- data.frame(ratio=c("Vd/Va","b/SDa"),
                      estimate=c(Vd/Va, b.mean/sqrt(Va)))
  return(ans)
} 

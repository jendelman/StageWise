#' Report dominance parameters
#'
#' Report dominance parameters
#' 
#' The dominance variance (Vd) and baseline heterosis (b) are quantified relative to additive variance (Va) and std. dev. (SDa), respectively. For single traits, the estimate and SE of the ratios are calculated based on the delta method (Rice 2007, p. 162-166). For a multi-trait/loc model, the result is just the ratio of the estimates, with \code{index.coeff} specifying the coefficients of the standardized true values (see also \code{\link{blup}}) and \code{gamma} specifying the relative weight for non-additive to additive genetic merit (see also \code{\link{gain}}). 
#' 
#' @param params list returned by \code{\link{Stage2}}
#' @param digits number of digits for rounding
#' @param index.coeff merit index coefficients
#' @param gamma contribution of non-additive values for genetic merit
#'
#' @return data frame with estimates and SE
#'
#' @references Rice JA (2007) Mathematical Statistics and Data Analysis, 3rd ed. Duxbury, Pacific Grove.
#' @export

dominance <- function(params, digits=2, index.coeff=NULL, gamma=0) {
  if (!is.element("heterosis",names(params)))
    stop("Dominance model was not used")

  id <- grep("dominance",params$vc$name,fixed=T)
  ia <- grep("additive",params$vc$name,fixed=T)
  vc <- params$vc[,-1]
  rownames(vc) <- params$vc$name

  nhet <- nrow(params$heterosis) 
  if (nhet > 1) {
    stopifnot(!is.null(index.coeff))
    if (colnames(params$heterosis)[1]=="trait") {
      traits <- params$heterosis$trait
      #Vd <- Va <- matrix(0,nrow=nhet,ncol=nhet,dimnames=list(traits,traits))
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
    
    ans <- data.frame(ratio=c("Vd/Va","b/SDa"),
                      estimate=round(c(Vd.mean/Va.mean,b.mean/sqrt(Va.mean)),digits))
    return(ans)
    
  } else {
  
    Vd.mean <- params$vc$estimate[id]
    Vd.var <- params$vc$SE[id]^2
    Va.mean <- params$vc$estimate[ia]
    Va.var <- params$vc$SE[ia]^2
    SDa.mean <- sqrt(Va.mean)*(1-2*Va.var/Va.mean^2)
    SDa.var <- Va.var/(4*Va.mean)
    b.mean <- params$heterosis$estimate
    b.var <- params$heterosis$SE^2
    
    ratio1.mean <- Vd.mean/Va.mean*(1+Va.var/Va.mean^2)
    ratio1.SE <- sqrt(Va.var*Vd.mean^2/Va.mean^2+Vd.var)/Va.mean
    ratio2.mean <- b.mean/SDa.mean*(1+SDa.var/SDa.mean^2)
    ratio2.SE <- sqrt(SDa.var*b.mean^2/SDa.mean^2+b.var)/SDa.mean
    
    ans <- data.frame(ratio=c("Vd/Va","b/SDa"),
                      estimate=round(c(ratio1.mean,ratio2.mean),digits),
                      SE=round(c(ratio1.SE,ratio2.SE),digits))
  
    return(ans)
  }
}

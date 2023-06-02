#' Report dominance parameters
#'
#' Report dominance parameters
#' 
#' The dominance variance (Vd) and baseline heterosis (b) are quantified relative to additive variance (Va) and std. dev. (SDa), respectively. The estimate and SE of the ratios are calculated based on the delta method (Rice 2006, p. 162-166). Only implemented for single trait/location models.
#' 
#' @param params list returned by \code{\link{Stage2}}
#' @param digits number of digits for rounding
#'
#' @return data frame with estimates and SE
#'
#' @references Rice JA (2007) Mathematical Statistics and Data Analysis, 3rd ed. Duxbury, Pacific Grove.
#' @export

dominance <- function(params, digits=2) {
  if (!is.element("heterosis",names(params)))
    stop("Dominance model was not used")
  if (nrow(params$heterosis) > 1)
    stop("Unavailable for multi-location/trait models")
  
  id <- grep("dominance",params$vc$name,fixed=T)
  ia <- grep("additive",params$vc$name,fixed=T)
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
  
  # if (length(ia) > 1) {
  #  ans <- list('Vd/Va'=data.frame(trait=params$heterosis$trait,
  #                           estimate=round(ratio1.mean,digits),
  #                           SE=round(ratio1.SE,digits)),
  #              'b/SDa'=data.frame(trait=params$heterosis$trait,
  #                            estimate=round(ratio2.mean,digits),
  #                            SE=round(ratio2.SE,digits)))
  # } else {
  ans <- data.frame(ratio=c("Vd/Va","b/SDa"),
                    estimate=round(c(ratio1.mean,ratio2.mean),digits),
                    SE=round(c(ratio1.SE,ratio2.SE),digits))
  
  return(ans)
}

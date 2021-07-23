#' S4 class for variances
#' 
#' @slot add additive
#' @slot g.resid genetic residual
#' @slot resid residual
#' @slot meanG mean of diagonal of G
#' @slot meanOmega mean of diagonal of Omega
#' @slot fixed.marker.var variance of marker fixed effects
#' @slot fixed.marker.cov contribution of marker fixed effects to additive covariance between locations
#' 
#' @importFrom methods new
#' @export
class_var <- setClass("class_var",slots=c(add="Matrix",g.resid="Matrix",resid="Matrix",
                                            meanG="numeric",meanOmega="numeric",
                                            fixed.marker.var="matrix",fixed.marker.cov="matrix"))

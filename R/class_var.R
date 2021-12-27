#' S4 class for variances
#' 
#' @slot add additive
#' @slot g.resid genetic residual
#' @slot resid residual
#' @slot meanG mean of diagonal of G
#' @slot meanOmega mean of diagonal of Omega
#' @slot fixed.marker.cov contribution of marker fixed effects to additive covariance between locations
#' 
#' @import methods
#' @import Matrix
#' @export
class_var <- setClass("class_var",slots=c(add="Matrix",g.resid="Matrix",resid="Matrix",
                                            meanG="numeric",meanOmega="matrix",
                                            fixed.marker.cov="array"))

#' S4 class for variances
#' 
#' @slot add additive
#' @slot g.iid iid genotype effect
#' @slot dom dominance 
#' @slot resid residual
#' @slot diagG average diagonal element of the G matrix
#' @slot diagD average diagonal element of the D matrix
#' @slot vars variances for reporting
#' @slot B var-cov matrix of fixed effects for gain
#' @slot fix.eff.marker names of fixed effect markers
#' 
#' @import methods
#' @import Matrix
#' @export
class_var <- setClass("class_var",slots=c(add="Matrix",
                                          g.iid="Matrix",
                                          dom="Matrix",
                                          resid="Matrix",
                                          diagG="numeric",
                                          diagD="numeric",
                                          vars="array",
                                          B="matrix",
                                          fix.eff.marker="character"))

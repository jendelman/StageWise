#' S4 class for variances
#' 
#' @slot geno1 first genetic effect
#' @slot geno2 second genetic effect
#' @slot model 0=no markers, 1=add, 2=add+g.resid, 3=add+dom
#' @slot resid residual
#' @slot vars variances for reporting
#' @slot B var-cov matrix of fixed effects for gain
#' @slot fix.eff.marker names of fixed effect markers
#' 
#' @import methods
#' @import Matrix
#' @export
class_var <- setClass("class_var",slots=c(geno1="Matrix",
                                          geno2="Matrix",
                                          model="integer",
                                          resid="Matrix",
                                          vars="array",
                                          B="matrix",
                                          fix.eff.marker="character"))

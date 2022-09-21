#' S4 class to prepare for blup
#' 
#' @slot id genotype identifiers
#' @slot ploidy ploidy
#' @slot var.u variance of random effects
#' @slot var.uhat variance of BLUPs
#' @slot avg.env average fixed effect of the environments
#' @slot heterosis regression coefficients for inbreeding
#' @slot fixed.marker fixed marker effects
#' @slot B var-cov matrix for fixed effects
#' @slot random random effect estimates
#' @slot geno1.var first var-cov matrix from \code{\link{class_var}}
#' @slot geno2.var second var-cov matrix from \code{\link{class_var}}
#' @slot model model from \code{\link{class_var}}
#' 
#' @import methods
#' @import Matrix
#' @export
class_prep <- setClass("class_prep",slots=c(id="character",ploidy="integer",var.u="Matrix",
                                            var.uhat="Matrix",avg.env="numeric",heterosis="numeric",
                                            fixed.marker="numeric",B="matrix",random="numeric",
                                            geno1.var="Matrix",geno2.var="Matrix",model="integer"))

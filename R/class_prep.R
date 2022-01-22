#' S4 class to prepare for blup
#' 
#' @slot id genotype identifiers
#' @slot var.u variance of random effects
#' @slot var.u.inv inverse of var.u
#' @slot var.uhat variance of BLUPs
#' @slot fixed fixed effect estimates
#' @slot random random effect estimates
#' @slot add matrix of additive variances from \code{\link{class_var}}
#' @slot loc.env data frame with loc, env
#' @slot trait.env data frame with trait, env
#' @slot index.scale sqrt of genetic variances for the locations/traits
#' @slot fixed.marker names of fixed effect markers
#' 
#' @import methods
#' @import Matrix
#' @export
class_prep <- setClass("class_prep",slots=c(id="character",var.u="Matrix",var.u.inv="Matrix",var.uhat="Matrix",
                                            fixed="numeric",random="numeric",add="Matrix",
                                            loc.env="data.frame",trait.env="data.frame",
                                            index.scale="numeric",fixed.marker="character"))

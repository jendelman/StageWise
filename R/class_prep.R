#' S4 class to prepare for blup
#' 
#' @slot y Stage 1 BLUEs
#' @slot Z Incidence matrix for random effects
#' @slot id genotype identifiers
#' @slot var.u variance of random effects
#' @slot Vinv inverted covariance matrix of the Stage 1 BLUEs
#' @slot Pmat P matrix from Searle
#' @slot fixed fixed effect estimates
#' @slot random random effect estimates
#' @slot add logical whether additive effects predicted
#' @slot loc.env data frame with loc, env
#' @slot index.scale sqrt of genetic variances for the locations/traits
#' @slot fixed.marker names of fixed effect markers
#' 
#' @import methods
#' @import Matrix
#' @export
class_prep <- setClass("class_prep",slots=c(y="numeric",Z="Matrix",id="character",var.u="Matrix",Vinv="Matrix",Pmat="Matrix",
                                            fixed="numeric",random="numeric",add="logical",loc.env="data.frame",
                                            index.scale="numeric",fixed.marker="character"))

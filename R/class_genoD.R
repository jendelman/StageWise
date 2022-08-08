#' S4 class for marker genotype data with dominance
#' 
#' @slot ploidy ploidy
#' @slot map Marker map positions
#' @slot coeff Coefficients of the additive marker effects (dim: indiv x marker)
#' @slot scale Scaling factor between markers and indiv for additive effects
#' @slot G Additive relationship matrix (from markers and potentially also pedigree)
#' @slot eigen.G list of eigenvalues and eigenvectors for G
#' @slot coeff.D coefficients of the dominance marker effects (dim: indiv x marker)
#' @slot scale.D Scaling factor between markers and indiv for dominance effects
#' @slot D Dominance relationship matrix
#' @slot eigen.D list of eigenvalues and eigenvectors for D
#' @slot Fg genomic inbreeding coefficient
#' 
#' @import methods
#' @import Matrix
#' @include class_geno.R
#' @export
class_genoD <- setClass("class_genoD",
                        slots=c(coeff.D="Matrix",scale.D="numeric",D="Matrix",eigen.D="list",Fg="numeric"),
                        contains = "class_geno")

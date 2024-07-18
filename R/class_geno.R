#' S4 class for marker genotype data
#' 
#' @slot ploidy If mixed ploidy, then a vector equal to pop size; otherwise a single integer
#' @slot map Marker map positions
#' @slot coeff Coefficients of the marker effects (dim: indiv x marker)
#' @slot scale Scaling factor between markers and indiv, vector of length equal to pop size
#' @slot G Additive relationship matrix (from markers and potentially also pedigree)
#' @slot eigen.G list of eigenvalues and eigenvectors
#' 
#' @import methods
#' @import Matrix
#' @export
class_geno <- setClass("class_geno",slots=c(ploidy="integer",map="data.frame",
                                            coeff="Matrix",scale="numeric",
                                            G="Matrix",eigen.G="list"))

#' S4 class for marker genotype data
#' 
#' @slot coeff Coefficients of the marker effects (dim: indiv x marker)
#' @slot scale Scaling factor between markers and indiv
#' @slot G Additive (genomic) relationship matrix
#' @slot Ginv Inverse G matrix
#' @slot map Marker map positions
#' 
#' @importFrom methods new
#' @import Matrix
#' @export
class_geno <- setClass("class_geno",slots=c(coeff="Matrix",scale="numeric",G="Matrix",Ginv="Matrix",map="data.frame"))

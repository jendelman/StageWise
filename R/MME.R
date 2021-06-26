#' S4 class for solving the mixed model equations
#' 
#' @slot data data frame with id, env, blue, trait (optional)
#' @slot kernels list of variance-covariance matrices for the genetic effects
#' @slot Rmat residual variance-covariance matrix

#' 
#' @importFrom methods new
#' @export
MME <- setClass("MME",slots=c(data="data.frame",kernels="list",Rmat="Matrix"))

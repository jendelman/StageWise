#' Genomic inbreeding coefficient
#'
#' Genomic inbreeding coefficient
#' 
#' Under the additive model, the inbreeding coefficient comes from the diagonal elements of the G matrix according to F = (G-1)/(ploidy-1). For dominance, the inbreeding coefficient is the scaled row-sum of the dominance coefficient matrix. 
#' 
#' @param geno object of \code{\link{class_geno}}
#'
#' @return data frame with F[G] and (when dominance is present) F[D]
#'
#' @export

inbreeding <- function(geno) {
  stopifnot(inherits(geno,"class_geno"))
  
  x <- data.frame(F.G=(diag(geno@G)-1)/(geno@ploidy-1))
  if (inherits(geno,"class_genoD"))
    x$F.D <- geno@Fg

  pops <- attr(geno@scale,"pop")
  np <- length(unique(pops))
  if (np > 1) {
    x$pop <- pops
  }
  rownames(x) <- rownames(geno@G)
  return(x)
}
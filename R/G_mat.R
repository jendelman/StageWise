#' Additive genomic relationships
#' 
#' Coefficients and relationship matrix for additive effects with bi-allelic markers
#' 
#' Additive effects are based on the traditional orthogonal decomposition of genetic variance in panmictic populations (Fisher 1918; Kempthorne 1957; Endelman et al. 2018). The G matrix is computed from the coefficients and scaling factor according to G = tcrossprod(coeff/scale). Missing genotype data is replaced with the population mean. 
#' 
#' @references Fisher (1918) Trans. Roy. Soc. Edin. 52:399-433.
#' @references Kempthorne (1957) An Introduction to Genetic Statistics.
#' @references Endelman et al. (2018) Genetics 209:77-87.
#'  
#' @param geno Matrix of allele dosages (markers x indiv)
#' @param ploidy 2 or 4
#' 
#' @return List containing
#' \describe{
#' \item{coeff}{Coefficients of the marker effects (dim: indiv x marker)}
#' \item{scale}{Scaling factor between markers and indiv}
#' \item{mat}{G matrix}
#' }
#' 
#' @export

G_mat <- function(geno,ploidy) {
  m <- nrow(geno)
  n <- ncol(geno)
  coeff <- scale(t(geno),scale=F)[1:n,1:m]
  coeff[is.na(coeff)] <- 0
  p <- apply(geno,1,mean,na.rm=T)/ploidy
  scale <- sqrt(sum(ploidy*p*(1-p)))
  return(list(coeff=coeff,scale=scale,mat=tcrossprod(coeff/scale)))
}

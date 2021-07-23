#' Read marker genotype data
#' 
#' Read marker genotype data
#' 
#' First three columns of the file: marker, chrom, position. Subsequent columns contain the allele dosage for individuals/clones, coded 0,1,2,...ploidy (fractional values are allowed). Additive coefficients are computed by subtracting the population mean from each marker, and the additive (genomic) relationship matrix is computed as G = tcrossprod(coeff)/scale. The scale parameter ensures the mean of the diagonal elements of G equals 1 under panmictic equilibrium. Missing genotype data is replaced with the population mean. For numerical conditioning, the inverse of G is computed after adding 1e-6 to the diagonal elements.
#'  
#' @param filename Name of CSV file
#' @param ploidy 2,4,6,etc. (even numbers)
#' @param map TRUE/FALSE
#' 
#' @return Variable of class \code{\link{class_geno}}.
#' 
#' @export
#' @importFrom utils read.csv

read_geno <- function(filename,ploidy,map) {
  data <- read.csv(file = filename,check.names=F)
  if (map) {
    map <- data[,1:3]
    colnames(map) <- c("marker","chrom","position")
    geno <- as.matrix(data[,-(1:3)])
  } else {
    map <- data.frame(marker=character(0),
                      chrom=character(0),
                      position=numeric(0))
    geno <- as.matrix(data[,-1])
  }
  rownames(geno) <- data[,1]
  id <- colnames(geno)
  coeff <- Matrix(scale(t(geno),scale=F),
                  dimnames=list(id,map$marker))
  coeff[which(is.na(coeff))] <- 0
  p <- apply(geno,1,mean,na.rm=T)/ploidy
  scale <- sum(ploidy*p*(1-p))
  G <- tcrossprod(coeff)/scale
  n <- nrow(G)
  Ginv <- solve(as(G+Diagonal(n=n)*1e-6,"dpoMatrix"))
  
  return(new(Class="class_geno",
             coeff=coeff,scale=scale,G=G,Ginv=Ginv,map=map))
}

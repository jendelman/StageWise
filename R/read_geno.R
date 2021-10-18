#' Read marker genotype data
#' 
#' Read marker genotype data
#' 
#' When \code{map=TRUE}, first three columns of the file are marker, chrom, position. When \code{map=FALSE}, the first column is marker. Subsequent columns contain the allele dosage for individuals/clones, coded 0,1,2,...ploidy (fractional values are allowed). The input file for diploids can also be coded using {-1,0,1} (fractional values allowed). Additive coefficients are computed by subtracting the population mean from each marker, and the additive (genomic) relationship matrix is computed as G = tcrossprod(coeff)/scale. The scale parameter ensures the mean of the diagonal elements of G equals 1 under panmictic equilibrium. Missing genotype data is replaced with the population mean. For numerical conditioning, eigenvalues of G smaller than \code{eigen.tol} are replaced by \code{eigen.tol}. Monomorphic markers are removed.
#'  
#' @param filename Name of CSV file
#' @param ploidy 2,4,6,etc. (even numbers)
#' @param map TRUE/FALSE
#' @param eigen.tol 1e-9
#' 
#' @return Variable of class \code{\link{class_geno}}.
#' 
#' @export
#' @importFrom utils read.csv

read_geno <- function(filename,ploidy,map,eigen.tol=1e-9) {
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
  m <- nrow(geno)
  if ((min(geno[1:min(m,1000),],na.rm=T)==-1) & (ploidy==2)) {
    geno <- geno + 1
  }
  rownames(geno) <- data[,1]
  id <- colnames(geno)
  p <- apply(geno,1,mean,na.rm=T)/ploidy
  ix <- which(p > 0 & p < 1)
  stopifnot(length(ix) > 0)
  geno <- geno[ix,]
  p <- p[ix]
  if (nrow(map)>0)
    map <- map[ix,]
  
  coeff <- Matrix(scale(t(geno),scale=F),
                  dimnames=list(id,rownames(geno)))
  coeff[which(is.na(coeff))] <- 0
  
  scale <- sum(ploidy*p*(1-p))
  G <- tcrossprod(coeff)/scale
  eigen.G <- eigen(G,symmetric=TRUE)
  class(eigen.G) <- "list"
  eigen.G$values <- ifelse(eigen.G$values < eigen.tol,
                         eigen.tol,eigen.G$values)
  eigen.G$vectors <- Matrix(eigen.G$vectors,dimnames=list(id,id))
  G <- tcrossprod(eigen.G$vectors%*%Diagonal(x=sqrt(eigen.G$values)))

  return(new(Class="class_geno",map=map,coeff=coeff,scale=scale,
             G=G,eigen.G=eigen.G))
}

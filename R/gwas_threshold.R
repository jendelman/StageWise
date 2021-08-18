#' Compute GWAS discovery threshold
#' 
#' Compute GWAS discovery threshold
#' 
#' Uses a Bonferroni-type correction based on an effective number of markers that accounts for LD  (Moskvina and Schmidt, 2008). 
#' 
#' @param geno object of \code{\link{class_geno}}
#' @param alpha genome-wide significance level
#' @param exclude.chrom chromosomes to exclude
#' @param n.core number of cores to use
#' 
#' @references Moskvina V, Schmidt KM (2008) On multiple-testing correction in genome-wide association studies. Genetic Epidemiology 32:567-573. doi:10.1002/gepi.20331
#' 
#' @return -log10(p) threshold
#' 
#' @export
#' @importFrom stats cor
#' @import Matrix
#' @importFrom methods as
#' @importFrom parallel makeCluster clusterExport parLapply stopCluster
#' 
gwas_threshold <- function(geno,alpha=0.05,exclude.chrom=NULL,n.core=1) {

	stopifnot(inherits(geno,"class_geno"))
	stopifnot(nrow(geno@map)>0)
	chrom <- as.character(unique(geno@map$chrom))
	if (!is.null(exclude.chrom)) {
	  stopifnot(exclude.chrom %in% chrom)
	  exclude.chrom <- as.character(exclude.chrom)
	  chrom <- setdiff(chrom,exclude.chrom)
	}
	n.chrom <- length(chrom)
	r2 <- vector("list",n.chrom)
	names(r2) <- chrom
	x <- split(geno@map$marker,geno@map$chrom)
	x <- x[names(x) %in% chrom]
	
	f1 <- function(markers) {
	  r2 <- as(cor(as.matrix(geno@coeff[,markers]))^2,"dspMatrix")
	  Keff(r2,alpha)
	}
	
	if (n.core==1) {
	  y <- sapply(x,f1) 
	} else {
	  cl <- makeCluster(n.core)
	  clusterExport(cl=cl,varlist = NULL)
	  y <- unlist(parLapply(cl=cl,X=x,fun=f1))
	  stopCluster(cl)
	}
	return(-log10(alpha/sum(y)))
}

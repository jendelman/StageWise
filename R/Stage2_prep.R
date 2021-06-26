#' Prepare data for single-trait, Stage 2 analysis of multi-environment trials
#' 
#' Prepare data for single-trait, Stage 2 analysis of multi-environment trials
#' 
#' Designed to prepare data for \code{\link{Stage2}} based on output from \code{\link{Stage1}}. The argument \code{blue.vcov} is a list containing a matrix for each environment; the matrix contains BLUEs in column 1 and their variance-covariance in the remaining columns. The rownames are the id and env concatenated with ":". 
#' 
#' @param blue.vcov list containing the blue.vcov matrix from Stage 1 
#' @param trait trait name 
#' @param exclude.id vector of individuals to exclude
#' 
#' @return a list containing
#' \describe{
#' \item{blue}{data frame of BLUEs}
#' \item{Omega}{variance-covariance matrix of BLUEs}
#' }
#' 
#' @export

Stage2_prep <- function(blue.vcov,trait,exclude.id=character(0)) {
  
  stopifnot(is.list(blue.vcov))
  stopifnot(!sapply(blue.vcov,is.null))
  n.env <- length(blue.vcov)
  vcov <- vector("list",n.env)
  blue.all <- NULL
  for (i in 1:n.env) {
    tmp <- strsplit(rownames(blue.vcov[[i]]),split=":",fixed=T)
    id.i <- sapply(tmp,"[",1)
    env.i <- sapply(tmp,"[",2)
    ix <- which(!(id.i %in% exclude.id))
    blue <- data.frame(id=id.i[ix],env=env.i[ix],blue=as.numeric(blue.vcov[[i]][ix,1]))
    colnames(blue) <- c("id","env",trait)
    blue.all <- rbind(blue.all,blue)
    vcov[[i]] <- blue.vcov[[i]][ix,1+ix]
  }
  Omega <- direct_sum(vcov)
  attr(Omega,"INVERSE") <- FALSE
  rownames(Omega) <- colnames(Omega) <- 1:nrow(Omega)
  return(list(blue=blue.all,Omega=Omega))
}


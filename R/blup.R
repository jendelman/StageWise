#' BLUP
#' 
#' BLUP
#' 
#' Argument \code{what="id"} leads to prediction of breeding values (BV) and genotypic values (GV), including the average fixed effect of the environments and any fixed effect markers.  For \code{what="marker"}, environment fixed effects are not included in the BLUP. Argument \code{index.coeff} should be a named vector, matching the names of the locations or traits. Index coefficients are assumed to imply relative weights for the different locations or traits; as such, they are divided by the square root of the genetic variance estimate for that location/trait and then rescaled to have unit sum.
#' 
#' @param data object of \code{\link{class_prep}} from \code{\link{blup_prep}}
#' @param geno object of \code{\link{class_geno}} from \code{\link{read_geno}}
#' @param what "id" or "marker"
#' @param index.coeff index coefficients for the locations or traits
#' @param gwas.ncore Integer indicating number of cores to use for GWAS (default is 0 for no GWAS). Requires \code{what="markers"}.
#' 
#' @return Data frames of BLUPs
#' 
#' @import Matrix
#' @importFrom stats formula pnorm
#' @importFrom parallel makeCluster clusterExport parCapply stopCluster
#' @export

blup <- function(data,geno=NULL,what,index.coeff=NULL,gwas.ncore=0L) {
  
  stopifnot(what %in% c("id","marker"))
  n.id <- length(data@id)
  n.mark <- length(data@fixed.marker)
  
  if (nrow(data@loc.env) > 0) {
    locations <- sort(unique(data@loc.env$loc))
    data@index.scale <- data@index.scale[locations]
    data@loc.env$loc <- factor(data@loc.env$loc,levels=locations)
    n.loc <- length(locations)
    if (is.null(index.coeff)) {
      index.coeff <- 1/data@index.scale
    } else {
      ix <- match(locations,names(index.coeff))
      if (any(is.na(ix))) {stop("Check names of locations in index.coeff")}
      index.coeff <- index.coeff[ix]/data@index.scale
    }
    index.coeff <- index.coeff/sum(index.coeff)
    n.env <- length(data@fixed) - n.loc*n.mark
    fix.value <- sum(index.coeff*tapply(data@fixed[data@loc.env$env],data@loc.env$loc,mean))
  } else {
    n.loc <- 1
    index.coeff <- 1
    n.env <- length(data@fixed) - n.loc*n.mark
    fix.value <- mean(data@fixed[1:n.env])
  }
  
  if (n.mark > 0) {
    a <- kronecker(Diagonal(n=n.mark),Matrix(index.coeff,nrow=1)) %*% 
      Matrix(data@fixed[n.env+1:(n.loc*n.mark)],ncol=1)
    fix.value <- fix.value + as.numeric(as.matrix(geno@coeff[,data@fixed.marker]) %*% a)
  }
  
  if (nrow(data@add) > 0 & is.null(geno)) {
    stop("Missing geno argument")
  }
  if (nrow(data@add)==0 & !is.null(geno)) {
    stop("Marker data was not used in blup_prep")
  }  
  
  M <- kronecker(Diagonal(n=n.id),Matrix(index.coeff,nrow=1))
  if (what=="id") {
    if (is.null(geno)) {
      out <- data.frame(id=data@id,GV=as.numeric(M%*%Matrix(data@random,ncol=1)) + fix.value)
      
      GZt <- tcrossprod(data@var.u,data@Z)
      var.uhat <- crossprod(tcrossprod(chol(data@Pmat),GZt))
      numer <- diag(M%*% tcrossprod(var.uhat,M))
      denom <- diag(M%*% tcrossprod(data@var.u,M))
      out$GV.r2 <- numer/denom
    } else {
      out <- data.frame(id=data@id,BV=as.numeric(cbind(M,0*M) %*% Matrix(data@random,ncol=1)) + fix.value,
                           GV=as.numeric(cbind(M,M) %*% Matrix(data@random,ncol=1)) + fix.value)
      GZt <- tcrossprod(data@var.u,cbind(data@Z,data@Z))
      var.uhat <- crossprod(tcrossprod(chol(data@Pmat),GZt))
      numer <- diag(cbind(M,0*M)%*% tcrossprod(var.uhat,cbind(M,0*M)))
      denom <- diag(cbind(M,0*M)%*% tcrossprod(data@var.u,cbind(M,0*M)))
      out$BV.r2 <- numer/denom
      numer <- diag(cbind(M,M)%*% tcrossprod(var.uhat,cbind(M,M)))
      denom <- diag(cbind(M,M)%*% tcrossprod(data@var.u,cbind(M,M)))
      out$GV.r2 <- numer/denom
    }
    return(out)
  } 
  
  #markers
  if (is.null(geno)) {
    stop("Cannot predict marker effects without genotype data")
  }
  
  m <- ncol(geno@coeff)
  if (n.loc > 1) {
     Ct <- data@Z %*% kronecker(geno@coeff,Matrix(index.coeff,ncol=1))
  } else {
     Ct <- data@Z %*% geno@coeff
  }
  V.alpha <- sum(index.coeff*diag(data@add))/geno@scale
  Ct <- Ct * V.alpha
  add.effect <- as.numeric(crossprod(Ct,data@Pmat %*% Matrix(data@y,ncol=1)))
  
  if (nrow(geno@map)==0) {
      out <- data.frame(marker=colnames(geno@coeff),add.effect=add.effect)
  } else {
      out <- data.frame(geno@map,add.effect=add.effect)
  }
  
  f.se <- function(x,P) {sqrt(crossprod(x,P%*%x))}

  if (gwas.ncore > 0) {
    if (gwas.ncore == 1) {
      se <- apply(Ct,2,f.se,P=data@Pmat)
      std.effect <- out$add.effect/se
    } else {
      cl <- makeCluster(gwas.ncore)
      clusterExport(cl=cl,varlist=NULL)
      se <- parCapply(cl=cl,x=Ct,f.se,P=data@Pmat)
      stopCluster(cl)
      std.effect <- out$add.effect/sapply(se,as.numeric)
    } 
    out$gwas.score <- -log10(pnorm(q=abs(std.effect),lower.tail=FALSE)*2)
  }
  
  if (n.mark > 0) {
    k <- match(data@fixed.marker,out$marker)
    out$add.effect[k] <- out$add.effect[k] + as.numeric(a)
  }
  
  return(out)
}

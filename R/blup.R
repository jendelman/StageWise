#' BLUP
#' 
#' BLUP
#' 
#' Argument \code{index.coeff} should be a named vector, matching the names of the locations or traits. Index coefficients are normalized to have unit sum.
#' 
#' @param data object of \code{\link{class_prep}} from \code{\link{blup_prep}}
#' @param geno object of \code{\link{class_geno}} from \code{\link{read_geno}}
#' @param what "id" or "marker"
#' @param index.coeff index coefficients for the locations or traits
#' @param gwas.ncore Number of cores to use for GWAS (default is 0 = no GWAS) if \code{what="markers"}
#' 
#' @return Data frames of BLUPs
#' 
#' @import Matrix
#' @importFrom stats formula pnorm
#' @importFrom parallel makeCluster clusterExport parCapply stopCluster
#' @export

blup <- function(data,geno=NULL,what,index.coeff=NULL,gwas.ncore=0) {
  
  stopifnot(what %in% c("id","marker"))
  fix.value <- mean(data@fixed)
  n.id <- length(data@id)
  
  #if (n.mark > 0) {
  #  a <- kronecker(Diagonal(n=n.mark),Matrix(index.coeff,nrow=1)) %*% beta[n.env+1:(n.loc*n.mark),]
  #  fix.value <- fix.value + as.numeric(as.matrix(geno@coeff[,fix.eff.markers]) %*% a)
  #}
  if (nrow(data@loc.env) > 0) {
    locations <- unique(data@loc.env$loc)
    n.loc <- length(locations)
    if (is.null(index.coeff)) {
      index.coeff <- rep(1/n.loc,n.loc)
      names(index.coeff) <- locations
    } else {
      index.coeff <- index.coeff/sum(index.coeff)
      ix <- match(locations,names(index.coeff))
      if (any(is.na(ix))) {
        stop("Check names of locations in index.coeff")
      }
      index.coeff <- index.coeff[ix]
    }
  } else {
    n.loc <- 1
    index.coeff <- 1
  }
  
  if (data@add & is.null(geno)) {
    stop("Missing geno argument")
  }
  if (!data@add & !is.null(geno)) {
    stop("Marker data was not used in blup_prep")
  }  
  M <- kronecker(Diagonal(n=n.id),Matrix(index.coeff,nrow=1))
  if (what=="id") {
    if (is.null(geno)) {
      out <- data.frame(id=data@id,GV=as.numeric(M%*%Matrix(data@random,ncol=1)) + fix.value)
      
      GZt <- tcrossprod(data@var.u,data@Z)
      var.uhat <- forceSymmetric(GZt %*% tcrossprod(data@Pmat,GZt))
      numer <- diag(M%*% tcrossprod(var.uhat,M))
      denom <- diag(M%*% tcrossprod(data@var.u,M))
      out$GV.r2 <- numer/denom
    } else {
      out <- data.frame(id=data@id,BV=as.numeric(cbind(M,0*M) %*% Matrix(data@random,ncol=1)) + fix.value,
                           GV=as.numeric(cbind(M,M) %*% Matrix(data@random,ncol=1)) + fix.value)
      GZt <- tcrossprod(data@var.u,cbind(data@Z,data@Z))
      var.uhat <- forceSymmetric(GZt %*% tcrossprod(data@Pmat,GZt))
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
  
# b <- Vinv %*% (data$BLUE - X %*% beta)
# m <- ncol(geno@coeff)
# if (n.loc > 1) {
#   Ct <- Z %*% kronecker(geno@coeff,Matrix(index.coeff,ncol=1)) 
# } else {
#   Ct <- Z %*% geno@coeff
# }
# if (nrow(geno@map)==0) {
#   out[[2]] <- data.frame(marker=colnames(geno@coeff),
#                          add.effect=as.numeric(crossprod(Ct,b)))
# } else {
#   out[[2]] <- data.frame(geno@map,
#                          add.effect=as.numeric(crossprod(Ct,b)))
# }
# 
# if (gwas.ncore > 0) {
#   VXW <- Vinv%*%X%*%W
#   if (gwas.ncore > 1) {
#     cl <- makeCluster(gwas.ncore)
#     clusterExport(cl=cl,varlist=NULL)
#     se <- parCapply(cl=cl,x=Ct,function(z){
#       sqrt(crossprod(z,Vinv%*%z) - crossprod(z,VXW%*%z))
#     })
#     std.effect <- out2[[2]]$add.effect/sapply(se,as.numeric)
#     stopCluster(cl)
#   } else {
#     se <- apply(Ct,2,function(z){
#       sqrt(crossprod(z,Vinv%*%z) - crossprod(z,VXW%*%z))
#     })
#     std.effect <- out2[[2]]$add.effect/se
#   }
#   out[[2]]$gwas.score <- -log10(pnorm(q=abs(std.effect),lower.tail=FALSE)*2)
# }
# return(out)
}

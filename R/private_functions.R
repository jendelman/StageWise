sigdig <- function(x,digits=3) {
  m <- digits - ceiling(log10(max(abs(x),na.rm=T)))
  if (m > 0) {
    round(x,m)
  } else {
    signif(round(x),digits)
  }
}

pvar <- function(mu=0,V=NULL,weights=NULL) {
  if (!is.null(V)) {
    n <- nrow(V)
  } else {
    n <- length(mu)
  }
  if (is.null(weights)) {    
    weights <- rep(1,n)
  } 
  weights <- weights/sum(weights)
  if (!is.null(V)) {
    x <- sum(diag(V)*weights) - matrix(weights,nrow=1)%*%V%*%matrix(weights,ncol=1) + sum(mu^2*weights) - sum(mu*weights)^2
  } else {
    x <- sum(mu^2*weights) - sum(mu*weights)^2
  }
  as.numeric(x)
}

kron <- function(eigen.A, B) {
  #returns Q such that QtQ is solution
  #eigen.B <- eigen(B)
  tmp <- svd(B)
  eigen.B <- list(values=tmp$d,vectors=tmp$u)
  V1 <- kronecker(Diagonal(x=sqrt(eigen.A$values)),Diagonal(x=sqrt(eigen.B$values)))
  V1.inv <- kronecker(Diagonal(x=1/sqrt(eigen.A$values)),Diagonal(x=1/sqrt(eigen.B$values)))
  V2 <- as(as(Matrix(eigen.B$vectors,dimnames=list(rownames(B),rownames(B))),"generalMatrix"),"unpackedMatrix")
  V3 <- kronecker(eigen.A$vectors,V2,make.dimnames=T)
  return(list(mat=tcrossprod(V1,V3), inv=tcrossprod(V1.inv,V3)))
}

Keff <- function(r2,alpha) {
  m <- nrow(r2)
  if (m > 1) {
    Q <- sqrt(r2)
    Q[upper.tri(Q,diag=T)] <- NA
    rmax <- apply(Q[-1,],1,max,na.rm=T)
    kappa <- sqrt(1-rmax^(-1.31*log10(alpha)))
    return(1+sum(kappa))
  } else {
    return(1)
  }
}

get_x <- function(map) {
  #takes a map with chrom and position and returns x axis values for plotting multiple chromosomes
  a <- tapply(map[,2],map[,1],max)
  n <- length(a)
  m <- tapply(map[,2],map[,1],length)
  b <- c(0,apply(array(1:(n-1)),1,function(k){sum(a[1:k])}))
  x <- map[,2] + rep(b,times=m)
  return(x)
}

coerce_dpo <- function(x) {
  x2 <- try(as(x,"dpoMatrix"),silent=TRUE)
  if (inherits(x2,"dpoMatrix")) {
    tmp <- try(chol(x2),silent=TRUE)
    if (inherits(tmp,"Cholesky")) {
      return(x2)
    }
  }
  
  d <- Diagonal(x=1/sqrt(diag(x)))
  x2 <- crossprod(d,x%*%d)
  eg <- eigen(x2,symmetric=TRUE)
  thresh <- .Machine$double.eps*10
  repeat {
    lambda <- ifelse(eg$values < thresh,thresh,eg$values)
    K <- Matrix(Diagonal(x=sqrt(diag(x))) %*% eg$vectors %*% Diagonal(x=sqrt(lambda)))
    x3 <- tcrossprod(K)
    dimnames(x3) <- dimnames(x)
    tmp <- chol(x3)
    if (inherits(tmp,"Cholesky")) {
      return(x3)
    }
    thresh <- thresh*10
  }
}

fu <- function(x) {
  n.trait <- nrow(x)
  V <- Variable(n.trait,n.trait,PSD=TRUE)
  solved <- solve(Problem(Minimize(norm(V-x,type="F"))),solver="SCS")
  V1 <- solved$getValue(V)
  dimnames(V1) <- dimnames(x)
  coerce_dpo(V1)
}

f.cov.trait <- function(vc,traits,us) {
  n.trait <- length(traits)
  if (us)  {
    cov.mat <- matrix(0,nrow=n.trait,n.trait)
    dimnames(cov.mat) <- list(traits,traits)
    tmp <- expand.grid(traits,traits)
    tmp <- apply(tmp,2,as.character)[,c(2,1)]
    tmp2 <- apply(tmp,1,paste,collapse=":")
    ix <- sapply(as.list(tmp2),grep,x=rownames(vc),fixed=T)
    iv <- which(sapply(ix,length)>0)
    ix <- unlist(ix[iv])
    tmp <- tmp[iv,]
    cov.mat[cbind(tmp[,1],tmp[,2])] <- vc[ix,1]
    cov.mat[upper.tri(cov.mat,diag=F)] <- cov.mat[lower.tri(cov.mat,diag=F)]
  } else {
    iu <- apply(array(traits),1,grep,x=rownames(vc),fixed=T)
    cov.mat <- diag(vc[iu,1])
    dimnames(cov.mat) <- list(traits,traits)
    if (vc[1,4]!="B") {
      cov.mat[1,2] <- cov.mat[2,1] <- vc[1,1]*sqrt(vc[iu[1],1]*vc[iu[2],1])
    }
  }
  return(coerce_dpo(cov.mat))
}

f.cov.loc <- function(vc,locs) {
  n.loc <- length(locs)
  vcnames <- rownames(vc)
  if (n.loc==2) {
    iu <- apply(array(locs),1,grep,x=vcnames,fixed=T)
    cov.mat <- diag(vc[iu,1])
    dimnames(cov.mat) <- list(locs,locs)
    cov.mat[1,2] <- cov.mat[2,1] <- vc[1,1]*sqrt(vc[iu[1],1]*vc[iu[2],1])
    fa.mat <- t(chol(cov.mat))
  } else {
    psi <- diag(vc[apply(array(paste0(locs,"!var")),1,grep,x=vcnames,fixed=T),1])
    fa.mat <- matrix(0,nrow=n.loc,ncol=2)
    rownames(fa.mat) <- locs
    ix <- apply(array(paste0(locs,"!fa1")),1,grep,x=vcnames,fixed=T)
    fa.mat[,1] <- vc[ix,1]
    ix <- apply(array(paste0(locs,"!fa2")),1,grep,x=vcnames,fixed=T)
    fa.mat[,2] <- vc[ix,1]
    cov.mat <- tcrossprod(fa.mat) + psi
  }
  
  #rotate
  tmp <- svd(fa.mat)
  fa.mat <- tmp$u %*% diag(tmp$d)
  #scale
  D <- diag(1/sqrt(diag(cov.mat)))
  dimnames(D) <- list(locs,locs)
  
  return(list(cov.mat=coerce_dpo(cov.mat),
              loadings=D%*%fa.mat))
}

cov_to_cor <- function(x) {
  d <- Diagonal(x=1/sqrt(diag(x)))
  x2 <- crossprod(d,x%*%d)
  dimnames(x2) <- dimnames(x)
  return(as.matrix(x2))
}

direct_sum <- function(x) {
  n <- length(x) 
  m <- sapply(x,nrow)
  m.cumulative <- apply(array(1:n),1,function(k){sum(m[1:k])})
  
  z <- expand.grid(col=1:m[1],row=1:m[1])
  out <- data.frame(row=z$row,col=z$col,value=as.vector(x[[1]]))
  if (n > 1) {
    for (i in 2:n) {
      z <- expand.grid(col=(1:m[i])+m.cumulative[i-1],row=(1:m[i])+m.cumulative[i-1])
      tmp <- data.frame(row=z$row,col=z$col,value=as.vector(x[[i]]))
      out <- rbind(out,tmp)
    }
  }
  out <- out[out$col <= out$row,]
  dimO <- max(m.cumulative)
  dname <- as.character(1:dimO)
  return(new("dsTMatrix",uplo="L",i=out[,1]-1L,j=out[,2]-1L,x=out[,3],
             Dim=c(dimO,dimO),Dimnames=list(dname,dname),factors=list()))
}

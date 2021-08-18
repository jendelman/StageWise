kron <- function(eigen.A,B) {
  eigen.B <- eigen(B)
  V <- as(Matrix(eigen.B$vectors,dimnames=list(rownames(B),rownames(B))),"dgeMatrix")
  half <- kronecker(eigen.A$vectors,V,make.dimnames=T) %*%
    kronecker(Diagonal(x=sqrt(eigen.A$values)),Diagonal(x=sqrt(eigen.B$values)))
  return(list(half=half,full=tcrossprod(half)))
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

uniplot <- function(z) {
  x <- seq(-1,1,by=0.01)
  y1 <- sqrt(1-x^2)
  y2 <- -y1
  plot.data <- data.frame(name=rownames(z),x=z[,1],y=z[,2])
  p <- ggplot(data=data.frame(x,y1,y2)) + geom_line(mapping=aes(x=.data$x,y=.data$y1),colour="red") + geom_line(mapping=aes(x=.data$x,y=.data$y2),colour="red") + theme_bw() + coord_fixed(ratio=1) +
    geom_point(mapping=aes(x=.data$x,y=.data$y),data=plot.data) + xlab("") + ylab("") + 
    theme(axis.text = element_blank(),axis.ticks = element_blank()) 
  p + geom_text_repel(aes(x=.data$x,y=.data$y,label=.data$name),nudge_x=-0.1,plot.data,segment.colour="blue",force=2)
}

coerce_dpo <- function(x) {
  d <- Diagonal(x=1/sqrt(diag(x)))
  x2 <- crossprod(d,x%*%d)
  eg <- eigen(x2,symmetric=TRUE)
  lambda <- ifelse(eg$values < .Machine$double.eps*2,.Machine$double.eps*2,eg$values)
  d <- Diagonal(x=sqrt(diag(x)))
  x3 <- tcrossprod(Matrix(d %*% eg$vectors %*% Diagonal(x=sqrt(lambda))))
  dimnames(x3) <- dimnames(x)
  return(x3)
}

partition <- function(x) {
  Va <- mean(x[upper.tri(x,diag=F)])
  Va <- max(0,Va)
  VaL <- mean(diag(x)) - Va
  return(c(Va,VaL))
}

f.id <- function(vc,lvls,keyword) {
  vcnames <- rownames(vc)
  ix1 <- grep(keyword,vcnames,fixed=T)
  ix2 <- lapply(as.list(lvls),grep,x=vcnames,fixed=T)
  ix <- sapply(ix2,intersect,y=ix1)
  return(vc[ix,1])
}

f.us.trait <- function(vc,traits) {
  tmp <- expand.grid(traits,traits)
  tmp <- apply(tmp,2,as.character)[,c(2,1)]
  tmp2 <- apply(tmp,1,paste,collapse=":")
  ix <- sapply(as.list(tmp2),grep,x=rownames(vc),fixed=T)
  iv <- which(sapply(ix,length)>0)
  ix <- unlist(ix[iv])
  tmp <- tmp[iv,]
  n.lvl <- length(traits)
  cov.mat <- matrix(0,nrow=n.lvl,n.lvl)
  dimnames(cov.mat) <- list(traits,traits)
  cov.mat[cbind(tmp[,1],tmp[,2])] <- vc[ix,1]
  cov.mat[upper.tri(cov.mat,diag=F)] <- cov.mat[lower.tri(cov.mat,diag=F)]
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
  D <- diag(1/sqrt(diag(cov.mat)))
  dimnames(D) <- list(locs,locs)
  return(list(cov.mat=coerce_dpo(cov.mat),plot=uniplot(D%*%fa.mat)))
}

cov_to_cor <- function(x) {
  d <- Diagonal(x=1/sqrt(diag(x)))
  x2 <- crossprod(d,x%*%d)
  dimnames(x2) <- dimnames(x)
  return(as.matrix(x2))
}

cor2cov <- function(x) {
  n <- nrow(x)-1
  rho.mat <- matrix(x[1],nrow=n,ncol=n)
  diag(rho.mat) <- 1
  scale.mat <- diag(sqrt(x[-1]))
  dimnames(scale.mat) <- list(rownames(x)[-1],rownames(x)[-1])
  y <- coerce_dpo(scale.mat%*% rho.mat%*%scale.mat)
  dimnames(y) <- list(rownames(x)[-1],rownames(x)[-1])
  return(y)
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

uniplot <- function(z) {
  x <- seq(-1,1,by=0.01)
  y1 <- sqrt(1-x^2)
  y2 <- -y1
  plot.data <- data.frame(name=rownames(z),x=z[,1],y=z[,2])
  p <- ggplot(data=data.frame(x,y1,y2)) + geom_line(mapping=aes(x=x,y=y1),colour="red") + geom_line(mapping=aes(x=x,y=y2),colour="red") + theme_bw() + coord_fixed(ratio=1) +
    geom_point(mapping=aes(x=x,y=y),data=plot.data) + xlab("") + ylab("") + 
    theme(axis.text = element_blank(),axis.ticks = element_blank()) 
  p + geom_text_repel(aes(x=x,y=y,label=name),plot.data)
}

f.id <- function(vc,lvls,keyword) {
  vcnames <- rownames(vc)
  ix1 <- grep(keyword,vcnames,fixed=T)
  ix2 <- lapply(as.list(lvls),grep,x=vcnames,fixed=T)
  ix <- sapply(ix2,intersect,y=ix1)
  return(vc[ix,1])
}

f.cov.trait <- function(vc,traits) {
  tmp <- expand.grid(traits,traits)
  tmp <- apply(tmp,2,as.character)[,c(2,1)]
  tmp <- tmp[tmp[,1] >= tmp[,2],]
  tmp2 <- apply(tmp,1,paste,collapse=":")
  ix <- sapply(as.list(tmp2),grep,x=rownames(vc),fixed=T)
  n.lvl <- length(traits)
  cov.mat <- matrix(0,nrow=n.lvl,n.lvl)
  dimnames(cov.mat) <- list(traits,traits)
  cov.mat[cbind(tmp[,1],tmp[,2])] <- vc[ix,1]
  cov.mat[upper.tri(cov.mat,diag=F)] <- cov.mat[lower.tri(cov.mat,diag=F)]
  return(cov.mat)
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
  return(list(cov.mat=cov.mat,plot=uniplot(D%*%fa.mat)))
}


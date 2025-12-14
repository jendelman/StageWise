#' Read marker genotype data
#' 
#' Read marker genotype data
#' 
#' When \code{map=TRUE}, first three columns of the file are marker, chrom, position. When \code{map=FALSE}, the first column is marker. Subsequent columns contain the allele dosage for individuals/clones, coded 0,1,2,...ploidy (fractional values are allowed). The input file for diploids can also be coded using {-1,0,1} (fractional values allowed). Additive coefficients are computed by subtracting the population mean from each marker, and the additive (genomic) relationship matrix is computed as G = tcrossprod(coeff)/scale. The scale parameter ensures the mean of the diagonal elements of G equals 1 under panmictic equilibrium. Missing genotype data is replaced with the population mean. 
#' 
#' G can be blended with the pedigree relationship matrix (A) by providing a pedigree data frame in \code{ped} and blending parameter \code{w}. The blended relationship matrix is H = (1-w)G + wA. The first three columns of \code{ped} are id, parent1, parent2. Missing parents must be coded NA. An optional fourth column in binary (0/1) format can be used to indicate which individuals should be included in the H matrix, but this option cannot be combined with dominance. If there is no fourth column, only genotyped individuals are included. If a vector of w values is provided, the function returns a list of \code{\link{class_geno}} objects. 
#' 
#' If the A matrix is not used, then G is blended with the identity matrix (times the mean diagonal of G) to improve numerical conditioning for matrix inversion. The default for w is 1e-5, which is somewhat arbitrary and based on tests with the vignette dataset. The D matrix is also blended with the identity matrix using 1e-5 for numerical conditioning.
#' 
#' When \code{dominance=FALSE}, non-additive effects are captured using a residual genetic effect, with zero covariance. If \code{dominance=TRUE}, a (digenic) dominance covariance matrix is used instead. 
#' 
#' The argument \code{min.minor.allele} specifies the minimum number of individuals that must contain the minor allele. Markers that do not meet this threshold are discarded.
#' 
#' Optional argument \code{pop.file} gives the name of a CSV file with two columns: id,pop. If the populations have different ploidy, this is indicated using a named vector for \code{ploidy}.
#'  
#' @param filename Name of CSV file with marker allele dosage
#' @param ploidy 2,4,6,etc. (even numbers)
#' @param map TRUE/FALSE
#' @param min.minor.allele threshold for marker filtering (see Details)
#' @param w blending parameter (see Details)
#' @param ped optional, pedigree data frame with 3 or 4 columns (see Details)
#' @param dominance TRUE/FALSE whether to include dominance covariance (see Details)
#' @param pop.file CSV file defining populations
#' 
#' @return Variable of class \code{\link{class_geno}}.
#' 
#' @export
#' @importFrom utils capture.output
#' @import Matrix
#' @importFrom AGHmatrix Amatrix
#' @importFrom data.table fread

read_geno <- function(filename, ploidy, map, min.minor.allele=5, 
                      w=1e-5, ped=NULL, dominance=FALSE, pop.file=NULL) {
  if (any(w < 1e-5))
    stop("Blending parameter should not be smaller than 1e-5")
  
  data <- fread(file = filename, check.names=F, sep=",", header=T)
  if (map) {
    map <- data.frame(data[,1:3])
    colnames(map) <- c("marker","chrom","position")
    geno <- as.matrix(data[,-(1:3)])
    rownames(geno) <- as.character(map$marker)
  } else {
    map <- data.frame(marker=character(0),
                      chrom=character(0),
                      position=numeric(0))
    geno <- as.matrix(data[,-1])
    tmp <- data.frame(data[,1])
    colnames(tmp) <- "marker"
    rownames(geno) <- as.character(tmp$marker)
  }
  
  m <- nrow(geno)
  id <- colnames(geno)
  
  n.ploidy <- length(unique(ploidy))
  
  if (!is.null(pop.file)) {
    pdat <- read.csv(pop.file,colClasses = c("character","character")) 
    colnames(pdat) <- c("id","pop")
    pops <- unique(pdat$pop)
    np <- length(pops)
    ix <- match(id,pdat$id,nomatch = 0)
    stopifnot(ix > 0)
    pdat <- pdat[ix,]
    if (n.ploidy > 1) {
      iu <- match(pops,names(ploidy),nomatch=0)
      stopifnot(iu>0)
      ploidy <- ploidy[pops]
    } else {
      ploidy <- rep(ploidy,np)
      names(ploidy) <- pops
    }
    pdat$ploidy <- ploidy[pdat$pop]
    
  } else {
    stopifnot(length(ploidy)==1)
    #check for recoding -1,0,1 to 0,1,2
    if ((min(geno[1:min(m,1000),],na.rm=T)==-1) & (ploidy==2)) {
      geno <- geno + 1
    }
    
    np <- 1
    pops <- "1"
    pdat <- data.frame(id=id,pop="1",ploidy=ploidy)
  }
  
  pmat <- t(geno)/pdat$ploidy
  p2 <- matrix(as.numeric(NA),nrow=m,ncol=np)
  colnames(p2) <- pops
  
  for (i in 1:np) {
    p2[,i] <- apply(pmat[pdat$pop==pops[i],],2,mean,na.rm=T)
  }
  popn <- as.numeric(table(pdat$pop)[pops])
  p <- apply(t(p2)*popn,2,sum)/sum(popn)
  flip <- which(p > 0.5)
  pmat[,flip] <- 1-pmat[,flip]
  n.minor <- apply(pmat,2,function(z){sum(z > 0.1,na.rm=T)})
  
  # n.minor <- apply(geno,1,function(z){
  #   p <- mean(z,na.rm=T)/ploidy
  #   tab <- table(factor(round(z),levels=0:ploidy))
  #   if (p > 0.5) {
  #     sum(tab) - tab[ploidy+1]
  #   } else {
  #     sum(tab) - tab[1]
  #   }
  # })
  ix <- which(n.minor >= min.minor.allele)
  m <- length(ix)
  stopifnot(m > 0)
  cat(sub("X",min.minor.allele,"Minor allele threshold = X genotypes\n"))
  cat(sub("X",m,"Number of markers = X\n"))
  geno <- geno[ix,]
  p2 <- p2[ix,,drop=F]
  n <- ncol(geno)
  cat(sub("X",n,"Number of genotypes = X\n"))
  
  #p <- apply(geno,1,mean,na.rm=T)/ploidy
  #p <- p[ix]
  if (nrow(map)>0)
    map <- map[ix,]
  
  coeff <- matrix(0,nrow=n,ncol=m,
                  dimnames=list(id,rownames(geno)))
  scale.param <- numeric(n)
  names(scale.param) <- id
  attr(scale.param,"pop") <- pdat$pop
  
  for (i in 1:np) {
    ix <- which(pdat$pop==pops[i])
    pmat <- matrix(1,ncol=1,nrow=popn[i]) %*% matrix(p2[,i],nrow=1)
    scale.param[ix] <- ploidy[i]*sum(p2[,i]*(1-p2[,i]))
    coeff[ix,] <- t(geno[,ix])-ploidy[i]*pmat
  }
  coeff[which(is.na(coeff))] <- 0
  coeff <- Matrix(coeff)
  
  #scale <- ploidy*sum(p*(1-p))
  #G <- tcrossprod(coeff)/scale
  G <- tcrossprod(coeff/sqrt(scale.param))
  
  nw <- length(w)
  H <- vector("list",nw)
  Hinv <- vector("list",nw)
  
  if (is.null(ped)) {
    for (i in 1:nw) {
      H[[i]] <- (1-w[i])*G + w[i]*mean(diag(G))*Diagonal(n=nrow(G))
    }
  } else {
    if (n.ploidy > 1)
      stop("Pedigree option not available for mixed ploidy")
    
    if (ncol(ped)==4) {
      if (dominance)
        stop("Dominance cannot be used with ungenotyped individuals.")
      
      colnames(ped) <- c("id","parent1","parent2","H")
      
      id3 <- ped$id[which(ped$H==1)]
      n3 <- length(id3)
      cat(sub("X",n3,"Individuals in the H matrix: X\n"))
    } else {
      colnames(ped) <- c("id","parent1","parent2")
    }
    for (i in 1:3) {
      ped[,i] <- as.character(ped[,i])
    }
    if (!all(id %in% ped$id)) {
      stop("Some genotyped individuals are not in the pedigree file.")
    }
    
    ped2 <- data.frame(id=1:nrow(ped),
                       parent1=match(ped$parent1,ped$id,nomatch=0),
                       parent2=match(ped$parent2,ped$id,nomatch=0))
    invisible(capture.output(A <- as(Amatrix(ped2,ploidy=ploidy[1]),"symmetricMatrix")))
    rownames(A) <- ped$id
    colnames(A) <- ped$id
    
    if (ncol(ped)==4) {
      id2 <- setdiff(id3,id)
      iv <- match(id3,c(id,id2))
    } else {
      id2 <- character(0)
    }
    n2 <- length(id2)
    A <- A[c(id,id2),c(id,id2)]
    A11 <- A[id,id]
    
    if (n2 > 0) {
      A.inv <- solve(A)
      A11.inv <- solve(A11)
      coeff <- coeff[intersect(id,id3),]
      scale.param <- scale.param[intersect(id,id3)]
    }
    
    for (i in 1:nw) {
      Gw <- (1-w[i])*G + w[i]*A11
      
      if (n2==0) {
        H[[i]] <- Gw
      } else {
        tmp <- as(bdiag(solve(Gw)-A11.inv,Matrix(0,nrow=n2,ncol=n2)) + A.inv,"symmetricMatrix")
        Hinv[[i]] <- tmp[iv,iv]
      }
    }
    #id <- c(id,id2)
  }
  
  if (n.ploidy > 1) {
    ploidy <- as.integer(pdat$ploidy)
    names(ploidy) <- pdat$id
    attr(ploidy,"pop") <- pdat$pop
  } else {
    ploidy <- as.integer(pdat$ploidy[1])
    names(ploidy) <- NULL
  }
  
  G.list <- vector("list",nw)
  for (i in 1:nw) {
    if (!is.null(H[[i]])) {
      eigen.G <- eigen(H[[i]],symmetric=TRUE)
      eigen.G$vectors <- Matrix(eigen.G$vectors,dimnames=list(id,id))
    } else {
      eigen.G <- eigen(Hinv[[i]],symmetric=TRUE)
      eigen.G$values <- 1/eigen.G$values
      eigen.G$vectors <- Matrix(eigen.G$vectors,dimnames=list(id3,id3))
      H[[i]] <- tcrossprod(eigen.G$vectors%*%Diagonal(n=nrow(Hinv[[i]]),x=sqrt(eigen.G$values)))
    }
    class(eigen.G) <- "list"
    
    G.list[[i]] <- new(Class="class_geno",ploidy=ploidy,map=map,coeff=coeff,
                       scale=scale.param,G=H[[i]], eigen.G=eigen.G)
  }
  
  if (dominance) {
    coeff.D <- matrix(0,nrow=n,ncol=m,dimnames=list(id,rownames(geno)))
    scaleD.param <- numeric(n)
    names(scaleD.param) <- id
    for (i in 1:np) {
      ix <- which(pdat$pop==pops[i])
      pmat <- matrix(1,ncol=1,nrow=popn[i]) %*% matrix(p2[,i],nrow=1)
      scaleD.param[ix] <- 4*choose(pdat$ploidy[ix],2)*sum(p2[,i]^2*(1-p2[,i])^2)
      #Pmat <- kronecker(matrix(p2[,i],nrow=1,ncol=m),matrix(1,ncol=1,nrow=popn[i]))
      coeff.D[ix,] <- -2*choose(pdat$ploidy[ix],2)*pmat^2 + 2*(pdat$ploidy[ix]-1)*pmat*t(geno[,ix]) - 
        t(geno[,ix])*(t(geno[,ix])-1)
    }
    
    coeff.D[is.na(coeff.D)] <- 0
    coeff.D <- Matrix(coeff.D)
    #scale.D <- 4*choose(ploidy,2)*sum(p^2*(1-p)^2)
    D <- tcrossprod(coeff.D/sqrt(scaleD.param))
    D <- (1-1e-5)*D + (1e-5)*mean(diag(D))*Diagonal(n=nrow(D))
    eigen.D <- eigen(D,symmetric=TRUE)
    eigen.D$vectors <- Matrix(eigen.D$vectors,dimnames=list(id,id))
    class(eigen.D) <- "list"
    Fg <- apply(coeff.D,1,sum)/(-scale.param*(pdat$ploidy-1))
    names(Fg) <- id
    output <- lapply(G.list,function(x){new(Class="class_genoD",ploidy=x@ploidy,map=x@map,
                                            coeff=x@coeff,scale=x@scale,G=x@G,eigen.G=x@eigen.G,
                                            coeff.D=coeff.D,scale.D=scaleD.param,D=D,eigen.D=eigen.D,
                                            Fg=Fg)})
  } else {
    output <- G.list
  }
  
  if (nw==1) {
    return(output[[1]])
  } else {
    return(output)
  }
}

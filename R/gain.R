#' Genetic gain
#' 
#' Genetic gain for breeding values
#' 
#' Optional argument \code{restricted} is a data frame with columns "trait" and "sign", where the options for sign are "=",">","<", representing equal to zero, non-negative, and non-positive.  
#' 
#' @param input either object of \code{\link{class_prep}} or quad.mat returned by this function
#' @param traits optional, plots ellipse tradeoff
#' @param coeff optional, index coefficients expressed in genetic standard deviation units
#' @param restricted data frame of restricted traits see Details
#' @param solver name of convex solver (default is "ECOS")
#' 
#' @return List containing
#' \describe{
#' \item{quad.mat}{quadratic matrix for the ellipsoid}
#' \item{plot}{ellipse plot}
#' \item{result}{data frame with gain and coefficients for the traits}
#' }
#' 
#' @import ggplot2
#' @importFrom ggforce geom_ellipse
#' @import CVXR

gain <- function(input,traits=NULL,coeff=NULL,restricted=NULL,solver="ECOS") {
  
  if (class(input)[1]=="class_prep") {
    # need to compute var.bhat
    n1 <- ncol(input@add)
    n2 <- ncol(input@g.iid)
    n.trait <- max(n1,n2)
    if (n.trait==1) 
      stop("Input is not multi-trait")
    
    B <- matrix(0,nrow=n.trait,ncol=n.trait)
    
    if (input@ploidy==0) {
      #no marker data

      nt <- nrow(input@var.uhat)
      n <- nt/n.trait
      ix <- split(1:nt,f=rep(1:n.trait,times=n))
      
      for (i in 1:n.trait)
        for (j in i:n.trait) {
          ixi <- ix[[i]]
          ixj <- ix[[j]]
          L <- input@var.uhat[ixi,ixj]
          B[i,j] <- B[j,i] <- mean(diag(L)) - mean(L)
        }
      B <- B + input@B
      trait.names <- colnames(input@g.iid)
      trait.scale <- sqrt(diag(input@g.iid))
    } else {
      
      #additive and non-additive
      nt <- nrow(input@var.uhat)/2
      n <- nt/n.trait
      ix <- split(1:nt,f=rep(1:n.trait,times=n))
      if (ncol(input@dom) > 1) {
        gamma <- (input@ploidy/2 - 1)/(input@ploidy - 1)
        trait.scale <- sqrt(diag(input@add)+gamma^2*diag(input@dom))
      } else {
        gamma <- 0
        trait.scale <- sqrt(diag(input@add))
      }
      M <- cbind(diag(n),gamma*diag(n))
      
      for (i in 1:n.trait)
        for (j in i:n.trait) {
          ixi <- ix[[i]]
          ixj <- ix[[j]]
          L <- M %*% input@var.uhat[c(ixi,nt+ixi),c(ixj,nt+ixj)] %*% t(M)
          B[i,j] <- B[j,i] <- mean(diag(L)) - mean(L)
        }
      B <- B + input@B
      trait.names <- colnames(input@add)
    }
    
    D <- diag(trait.scale)
    quad.mat <- D%*%solve(B)%*%D
    dimnames(quad.mat) <- list(trait.names,trait.names)
  } else {
    quad.mat <- input
  } 
  
  if (!is.null(coeff)) {
    
    iv <- match(trait.names,names(coeff),nomatch=0)
    if (any(iv==0))
      stop("Trait names in coefficients do not match input")
    coeff <- coeff[iv]
    
    x <- Variable(n.trait)
    constraints <- list(quad_form(x,quad.mat) <= 1)
    
    if (!is.null(restricted)) {
      rest <- match(restricted$trait,trait.names,nomatch=0)
      if (any(rest==0))
        stop("Restricted trait names do not match input")
      nr <- length(rest)
      if (nr > 0) 
        coeff[rest] <- 0
      
      ix <- which(restricted$sign=="=")
      nix <- length(ix)
      if (nix > 0) {
        A <- matrix(0,nrow=length(ix),ncol=n.trait)
        A[cbind(1:nix,match(restricted$trait[ix],names(coeff)))] <- 1
        constraints <- c(constraints,list(A%*%x==0))
      }
      
      ix <- which(restricted$sign=="<")
      nix <- length(ix)
      if (nix > 0) {
        A <- matrix(0,nrow=length(ix),ncol=n.trait)
        A[cbind(1:nix,match(restricted$trait[ix],names(coeff)))] <- 1
        constraints <- c(constraints,list(A%*%x<=0))
      }
      
      ix <- which(restricted$sign==">")
      nix <- length(ix)
      if (nix > 0) {
        A <- matrix(0,nrow=length(ix),ncol=n.trait)
        A[cbind(1:nix,match(restricted$trait[ix],names(coeff)))] <- 1
        constraints <- c(constraints,list(A%*%x>=0))
      }
    } 
    
    v <- matrix(coeff,nrow=1)
    objective <- Maximize(v%*%x)
    problem <- Problem(objective,constraints)
    result <- solve(problem,solver=solver)
    if (result$status!="optimal")
      stop("Convex solver failed")
    x.opt <- as.numeric(result$getValue(x))
    c.opt <- as.numeric(quad.mat %*% x.opt)
    c.opt <- c.opt/sqrt(sum(c.opt^2))
    names(c.opt) <- names(x.opt) <- trait.names
    result <- data.frame(trait=trait.names,
                         response=round(x.opt,3),
                         coeff=round(c.opt,3))
    out <- list(quad.mat=quad.mat,result=result)
  } else {
    out <- list(quad.mat=quad.mat)
  }
  
  if (!is.null(traits)) {
    #plot ellipse
    ix <- match(traits,rownames(quad.mat),nomatch=0)
    if (any(ix==0))
      stop("Trait names not found")
    
    eg <- eigen(quad.mat[ix,ix])
    lens <- 1/sqrt(eg$values)
    angle <- atan(eg$vectors[2,2]/eg$vectors[1,2])
    p <- ggplot() + geom_ellipse(aes(x0=0,y0=0, a = lens[2], b=lens[1], angle = angle)) + coord_fixed() + theme_bw() + xlab(traits[1]) + ylab(traits[2]) 
    
    if (!is.null(coeff)) {
      if (all(c.opt[ix]!=0)) {
        angle2 <- atan(c.opt[ix[2]]/c.opt[ix[1]])
        if (c.opt[ix[2]] > 0 & angle2 < 0)
          angle2 <- pi + angle2
        x <- matrix(c(cos(angle2),sin(angle2)),ncol=1)
        r <- 1/as.numeric(sqrt(t(x) %*% quad.mat[ix,ix] %*% x))
        p <- p + geom_spoke(aes(x=0,y=0,angle=angle2,radius=r),col="red",lty=2)
      } else {
        
        #to do
        p <- p + geom_segment(aes(x=0,y=0,xend=0,yend=1),linetype=2,col="red")   
        
      }
      p <- p + geom_segment(aes(x=0,y=0,xend=x.opt[ix[1]],yend=x.opt[ix[2]]),col="blue") + geom_point(aes(x=x.opt[ix[1]],y=x.opt[ix[2]]),col="blue")
    }
    out <- c(out,list(plot=p))
  }
  return(out)
}
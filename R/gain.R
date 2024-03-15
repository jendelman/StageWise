#' Genetic gain
#' 
#' Genetic gain calculations
#' 
#' Either \code{merit} or \code{desired} can be used, not both. The former specifies the relative contribution of each trait to genetic merit, while the latter specifies the relative desired gain in genetic standard deviation units. All traits must be specified. Optional argument \code{restricted} is a data frame with columns "trait" and "sign", where the options for sign are "=",">","<", representing equal to zero, non-negative, and non-positive. When \code{desired} is used, the \code{restricted} argument is ignored.
#' 
#' The argument \code{gamma} controls the definition of genetic merit. (See notation in the journal publication.) The default is NULL, which implies breeding values. For purely additive values, use gamma = 0. For total genotypic value, use gamma = 1.
#' 
#' Note that this function assumes a selection index of BLUPs, not phenotypes. 
#' 
#' @param input either object of \code{\link{class_prep}} or quad.mat returned by this function
#' @param merit named vector of merit coefficients, in genetic standard deviation units
#' @param desired named vector of desired gains, in genetic standard deviation units
#' @param restricted data frame of restricted traits, see Details
#' @param traits optional vector with exactly 2 trait names, to plot elliptical response
#' @param gamma contribution of non-additive values for genetic merit
#' @param solver name of convex solver (default is "ECOS")
#' 
#' @return List containing
#' \describe{
#' \item{quad.mat}{quadratic matrix for the ellipsoid}
#' \item{plot}{ellipse plot}
#' \item{table}{data frame with the response and index coefficients for all traits}
#' }
#' 
#' @import ggplot2
#' @importFrom ggforce geom_ellipse
#' @import CVXR
#' @export

gain <- function(input, merit=NULL, desired=NULL, restricted=NULL,
                 traits=NULL, gamma=NULL, solver="ECOS", ...) {
  
  vararg <- list(...)
  if ("coeff" %in% names(vararg)) {
    merit <- vararg$coeff
  }
  stopifnot(is.null(merit) | is.null(desired))
  if (!is.null(traits))
    stopifnot(length(traits)==2)
  
  if (inherits(input,"class_prep")) {
    # need to compute var.bhat
    n.trait <- ncol(input@geno1.var)
    if (n.trait==1) 
      stop("Input is not multi-trait")
    
    B <- matrix(0,nrow=n.trait,ncol=n.trait)
    
    if (input@model < 2L) {
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
      trait.names <- colnames(input@geno1.var)
      trait.scale <- sqrt(diag(input@geno1.var))
    } else {
      
      #additive and non-additive
      nt <- nrow(input@var.uhat)/2
      n <- nt/n.trait
      ix <- split(1:nt,f=rep(1:n.trait,times=n))
      if (is.null(gamma)) {
        if (input@model==3L) {
          gamma <- (input@ploidy/2 - 1)/(input@ploidy - 1)  
        } else {
          gamma <- 0
        }
      }
      trait.scale <- sqrt(diag(input@geno1.var)+gamma^2*diag(input@geno2.var))
      M <- cbind(diag(n),gamma*diag(n))
      
      for (i in 1:n.trait)
        for (j in i:n.trait) {
          ixi <- ix[[i]]
          ixj <- ix[[j]]
          L <- M %*% input@var.uhat[c(ixi,nt+ixi),c(ixj,nt+ixj)] %*% t(M)
          B[i,j] <- B[j,i] <- mean(diag(L)) - mean(L)
        }
      B <- B + input@B
      trait.names <- colnames(input@geno1.var)
    }
    
    D <- diag(trait.scale)
    quad.mat <- D%*%solve(B)%*%D
    dimnames(quad.mat) <- list(trait.names,trait.names)
  } else {
    quad.mat <- input
    trait.names <- rownames(quad.mat)
    n.trait <- length(trait.names)
  } 
  
  if (!is.null(merit)) {
    
    iv <- match(trait.names,names(merit),nomatch=0)
    if (any(iv==0))
      stop("Trait names in merit coefficients do not match input")
    merit <- merit[iv]
    
    x <- Variable(n.trait)
    constraints <- list(quad_form(x,quad.mat) <= 1)
    
    if (!is.null(restricted)) {
      stopifnot(colnames(restricted)==c("trait","sign"))
      rest <- match(restricted$trait,trait.names,nomatch=0)
      if (any(rest==0))
        stop("Restricted trait names do not match input")
      nr <- length(rest)
      if (nr > 0) 
        merit[rest] <- 0
      
      ix <- which(restricted$sign=="=")
      nix <- length(ix)
      if (nix > 0) {
        A <- matrix(0,nrow=length(ix),ncol=n.trait)
        A[cbind(1:nix,match(restricted$trait[ix],names(merit)))] <- 1
        constraints <- c(constraints,list(A%*%x==0))
      }
      
      ix <- which(restricted$sign=="<")
      nix <- length(ix)
      if (nix > 0) {
        A <- matrix(0,nrow=length(ix),ncol=n.trait)
        A[cbind(1:nix,match(restricted$trait[ix],names(merit)))] <- 1
        constraints <- c(constraints,list(A%*%x<=0))
      }
      
      ix <- which(restricted$sign==">")
      nix <- length(ix)
      if (nix > 0) {
        A <- matrix(0,nrow=length(ix),ncol=n.trait)
        A[cbind(1:nix,match(restricted$trait[ix],names(merit)))] <- 1
        constraints <- c(constraints,list(A%*%x>=0))
      }
    } 
    
    v <- matrix(merit,nrow=1)
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
                         index=round(c.opt,3))
    rownames(result) <- NULL
    out <- list(quad.mat=quad.mat,table=result)
  } else {
    if (!is.null(desired)) {
      
      iv <- match(trait.names,names(desired),nomatch=0)
      if (any(iv==0))
        stop("Trait names in desired gains do not match input")
      desired <- desired[iv]
      c.opt <- as.numeric(quad.mat %*% matrix(desired,ncol=1))
      desired <- desired/sqrt(sum(desired*c.opt))
      c.opt <- c.opt/sqrt(sum(c.opt^2))
      result <- data.frame(trait=trait.names,
                           response=round(desired,3),
                           index=round(c.opt,3))
      rownames(result) <- NULL
      out <- list(quad.mat=quad.mat,table=result)
    } else {
      out <- list(quad.mat=quad.mat)
    }
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
    
    if (!is.null(merit)) {
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
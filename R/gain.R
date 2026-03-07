#' Genetic gain
#' 
#' Genetic gain calculations
#' 
#' The standardized genetic response vector, \eqn{x=H^{-1}g}{x=H^-1 g}, is constrained to an ellipsoid with equation \eqn{x'Qx = i^2}{x'Qx=i^2}, where i is selection intensity and H is a diagonal matrix of genetic standard deviations. The matrix of the quadratic form \eqn{Q=U'U}{Q=U'U}, where \eqn{U=RT^{-1}H}{U=R(T^-1)H}; R is the upper triangular cholesky factor of P = variance-covariance matrix of the selection values (e.g., BLUPs or phenotypes); and T = covariance between the genetic and selection values. For multi-trait BLUPs T=P, while under phenotypic selection T is the genetic variance-covariance matrix.
#' 
#' The first argument \code{input} can be a list of three matrices (P,T,H) or a variable of class \code{\link{class_prep}}, in which case the matrices are computed internally. 
#'
#' Either \code{merit} or \code{desired} can be used, not both. The former specifies the relative contribution of each trait to genetic merit, while the latter specifies the relative desired gain in genetic standard deviation units. All traits must be specified. Optional argument \code{restricted} is a data frame with columns "trait" and "sign", where the options for sign are "=",">","<", representing equal to zero, non-negative, and non-positive. When \code{desired} is used, the \code{restricted} argument is ignored.
#' 
#' The argument \code{gamma} controls the definition of genetic merit. (See notation in the journal publication.) The default is NULL, which implies breeding values. For purely additive values, use gamma = 0. For total genotypic value, use gamma = 1.
#' 
#' In the figure, the red dot is the intersection of the merit or desired gain vector with the ellipsoid, which is also the response under the desired gains index. The optimal merit response is shown with a blue dot.
#' 
#' In the function output, response values are in units of \eqn{i \sigma_g}{i sigma_g}, so multiplying by i gives the response for different intensities.
#' 
#' This function was originally designed to complement \code{\link{blup}}, which is based on an index of the form \eqn{I=b'BLUP[g]}{I=b'BLUP[g]}. For merit function \eqn{M=a'H^{-1}g}{M=a'(H^-1)g}, the index coefficients under BLUP are \eqn{b=H^{-1}a}{b=(H^-1)a}. If \eqn{c}{c} is the vector of coefficients passed to \code{\link{blup}}, internally the function uses \eqn{b=H^{-1}c}{b=(H^-1)c}. In other words, the argument should be scaled like \eqn{a}{a}. To return index coefficients on this scale, use \code{scale.coeff}=TRUE. To return the unscaled coefficients, use \code{scale.coeff}=FALSE. The default value, NULL, will use TRUE when the input is of class \code{\link{class_prep}} and FALSE when the input is a list of three matrices.
#' 
#' @param input variable of \code{\link{class_prep}} or list of matrices P,T,H (see Details)
#' @param merit named vector of merit coefficients, in genetic standard deviation units
#' @param desired named vector of desired gains, in genetic standard deviation units
#' @param restricted data frame of restricted traits, see Details
#' @param traits optional vector with exactly 2 trait names, to plot elliptical response
#' @param gamma contribution of non-additive values for genetic merit
#' @param scale.coeff should index coeff be scaled by genetic sd (see Details)
#' @param solver name of convex solver (default is "ECOS")
#' 
#' @return List containing
#' \describe{
#' \item{matrices}{for input \code{\link{class_prep}}, a list with P,T,H}
#' \item{plot}{ellipse plot}
#' \item{table}{data frame with the response and index coefficients for all traits}
#' }
#' 
#' @import ggplot2
#' @importFrom ggforce geom_ellipse
#' @import CVXR
#' @export

gain <- function(input, merit=NULL, desired=NULL, restricted=NULL,
                 traits=NULL, gamma=NULL, scale.coeff=NULL, solver="ECOS", ...) {
  
  intensity <- 1
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
    
    H <- diag(trait.scale)
    dimnames(H) <- dimnames(B) <- list(trait.names,trait.names)
    P <- Tmat <- B
  } else {
    stopifnot(inherits(input,"list") & length(input)==3)
    P <- input[[1]]
    trait.names <- rownames(P)
    Tmat <- input[[2]][trait.names,trait.names]
    H <- input[[3]][trait.names,trait.names]
    n.trait <- length(trait.names)
  } 
  Ginv <- solve(Tmat)
  U <- chol(P) %*% Ginv %*% H
  quad.mat <- crossprod(U)
  dimnames(quad.mat) <- list(trait.names,trait.names)
  
  if (!is.null(merit)) {
    iv <- match(trait.names,names(merit),nomatch=0)
    if (any(iv==0))
      stop("Trait names in merit coefficients do not match input")
    merit <- merit[iv]
    
    x <- Variable(n.trait)
    constraints <- list(quad_form(x,quad.mat) <= intensity^2)
    
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
  } 
  if (!is.null(desired)) {
    iv <- match(trait.names,names(desired),nomatch=0)
    if (any(iv==0))
      stop("Trait names in desired gains do not match input")
    desired <- desired[iv]
    x.opt <- desired*intensity/sqrt(as.numeric(crossprod(desired, quad.mat%*%desired)))
  } 
  if (inherits(input,"class_prep")) {
    out <- list(matrices=list(P=P,T=Tmat,H=H))
    if (is.null(scale.coeff))
      scale.coeff <- TRUE
  } else {
    out <- list()
    if (is.null(scale.coeff))
      scale.coeff <- FALSE
  }
  if (!is.null(desired)|!is.null(merit)) {
    c.opt <- crossprod(Ginv,H%*%x.opt)/intensity
    if (scale.coeff)
      c.opt <- c.opt*diag(H)
    c.opt <- c.opt/sqrt(sum(c.opt^2))
    result <- data.frame(trait=trait.names,
                         response=round(x.opt,3),
                         index.coeff=round(c.opt,3))
    rownames(result) <- NULL
    out <- c(out,list(table=result))
  } 

  if (!is.null(traits)) {
    #plot ellipse
    ix <- match(traits,trait.names,nomatch=0)
    if (any(ix==0))
      stop("Trait names not found")
    
    eg <- eigen(quad.mat[ix,ix])
    lens <- intensity/sqrt(eg$values)
    angle <- atan(eg$vectors[2,2]/eg$vectors[1,2])
    p <- ggplot() + geom_ellipse(aes(x0=0,y0=0, a = lens[2], b=lens[1], angle = angle)) + 
      coord_fixed() + theme_bw() + xlab(traits[1]) + ylab(traits[2]) 
    
    if (!is.null(desired)) {
      p <- p + geom_segment(aes(x=0,y=0,xend=x.opt[ix[1]],yend=x.opt[ix[2]]),col="red",lty=2) + 
          geom_point(aes(x=x.opt[ix[1]],y=x.opt[ix[2]]),col="red")
    }
    if (!is.null(merit)) {
      x2 <- merit*intensity/sqrt(as.numeric(crossprod(merit, quad.mat%*%merit)))
      p <- p + geom_segment(aes(x=0,y=0,xend=x2[ix[1]],yend=x2[ix[2]]),col="red",lty=2) + 
        geom_point(aes(x=x2[ix[1]],y=x2[ix[2]]),col="red")
      p <- p + geom_segment(aes(x=0,y=0,xend=x.opt[ix[1]],yend=x.opt[ix[2]]),col="blue") + 
        geom_point(aes(x=x.opt[ix[1]],y=x.opt[ix[2]]),col="blue")
    }
    out <- c(out,list(plot=p))
  }
  return(out)
}
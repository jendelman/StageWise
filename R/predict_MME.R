#' Compute BLUPs by solving the Mixed Model Equations
#' 
#' Compute BLUPs by solving the Mixed Model Equations
#' 
#' Use the function \code{\link{Stage2}} to create the object of class \code{\link{MME}}. BLUPs are computed at the average value of the fixed effects. If \code{weights} is used, the names must exactly match the names of the kernels in \code{data}. Using the argument \code{mask}, the phenotypes for a subset of the population can be masked, to enable cross-validation.
#' 
#' @param data variable of class \code{\link{MME}}
#' @param weights named vector of weights for the genetic effects in BLUP. Default is 1 for all effects.
#' @param mask (optional) data frame with column "id" and optional columns "env", "trait" 
#' 
#' @return For single trait analysis, function returns a data frame with columns: id,blup,r2. For multi-trait analysis, a list is returned containing
#' \describe{
#' \item{blup}{data frame of blups}
#' \item{r2}{data frame of reliabilities}
#' }
#' 
#' @import Matrix
#' @importFrom stats model.matrix
#' @export
#' 
predict_MME <- function(data,weights=NULL,mask=NULL) {
  
  n.random <- length(data@kernels)
  stopifnot(n.random > 0)
  if (is.null(weights)) {
    weights <- rep(1,n.random)
    names(weights) <- names(data@kernels)
  }
  stopifnot(any(weights > 0))
  stopifnot(n.random==length(weights))
  stopifnot(setequal(names(weights),names(data@kernels)))
  weights <- Matrix(weights[match(names(weights),names(data@kernels))],nrow=1,ncol=n.random)
  
  if (!is.null(mask)) {
    stopifnot(is.element("id",colnames(mask)))
    ix <- which(as.character(data@data$id) %in% as.character(mask$id))
    if ("env" %in% colnames(mask)) {
      ix2 <- which(as.character(data@data$env) %in% as.character(mask$env))
      ix <- intersect(ix,ix2)
    }
    if ("trait" %in% colnames(mask)) {
      ix2 <- which(as.character(data@data$trait) %in% as.character(mask$trait))
      ix <- intersect(ix,ix2)
    }
    data@data <- data@data[-ix,]
    data@Rmat <- data@Rmat[-ix,-ix]
  }

  if ("trait" %in% colnames(data@data)) {
    multi.trait <- TRUE
    data@data$trait <- factor(as.character(data@data$trait))
    traits <- levels(data@data$trait)
    n.trait <- length(traits)
    Z <- Matrix(model.matrix(~id:trait-1,data@data))
    colnames(Z) <- sub("trait","",colnames(Z),fixed=T)
    X <- Matrix(model.matrix(~env:trait-1,data@data))
    colnames(X) <- sub("trait","",colnames(X),fixed=T)
  } else {
    multi.trait <- FALSE
    Z <- Matrix(model.matrix(~id-1,data@data))
    X <- Matrix(model.matrix(~env-1,data@data))
  }
  colnames(Z) <- sub("id","",colnames(Z),fixed=T)
  colnames(X) <- sub("env","",colnames(X),fixed=T)
  
  Z2 <- kronecker(Matrix(1,nrow=1,ncol=n.random),Z)
  RX <- solve(data@Rmat,X)
  RZ <- solve(data@Rmat,Z2)
  Ry <- solve(data@Rmat,data@data$blue)
  Q11 <- crossprod(X,RX)
  Q12 <- crossprod(X,RZ)
  Q22 <- crossprod(Z2,RZ)+direct_sum(lapply(data@kernels,solve))
  Q <- rbind(cbind(Q11,Q12),cbind(t(Q12),Q22))
  b1 <- as.numeric(crossprod(X,Ry))
  b2 <- as.numeric(crossprod(Z2,Ry))
  b <- c(b1,b2)
  C <- solve(Q)
  ans <- as.numeric(C%*%b)
  
  n.fix <- ncol(X)
  C22 <- C[-(1:n.fix),-(1:n.fix)]
  
  tmp <- apply(X,2,mean)*ans[1:n.fix]
  id <- levels(data@data$id)
  n <- length(id)
  M <- kronecker(weights,Diagonal(n))
  if (multi.trait) {
    trait.X <- sapply(strsplit(colnames(X),split=":",fixed=T),"[",2)
    avg.fixed <- lapply(as.list(traits),function(trait){sum(tmp[trait.X %in% trait])*n.trait})
    
    trait.Z <- sapply(strsplit(colnames(Z),split=":",fixed=T),"[",2)
    trait.Z <- rep(trait.Z,n.random)
    iz <- lapply(as.list(traits),function(trait){which(trait.Z %in% trait)})
    
    blup <- lapply(iz,function(iz){as.numeric(M%*%ans[-(1:n.fix)][iz])})
    blup <- mapply(FUN=function(blup,avg.fixed){blup+avg.fixed},blup=blup,avg.fixed=avg.fixed)
    colnames(blup) <- traits
    
    numer <- lapply(iz,function(iz){diag(M%*%C22[iz,iz]%*%t(M))})
    denom <- lapply(iz,function(iz){diag(M%*%direct_sum(data@kernels)[iz,iz]%*%t(M))})
    r2 <- mapply(FUN=function(numer,denom){1-numer/denom},numer=numer,denom=denom)
    colnames(r2) <- traits
    
    return(list(blup=data.frame(id=as.character(id),blup),r2=data.frame(id=as.character(id),r2)))
  } else {
    avg.fixed <- sum(tmp)
    blup <- as.numeric(M%*%ans[-(1:n.fix)])+avg.fixed
    
    numer <- diag(M%*%C22%*%t(M))
    denom <- diag(M%*%direct_sum(data@kernels)%*%t(M))
    r2 <- 1-numer/denom
    return(data.frame(id=as.character(id),blup=blup,r2=r2))
  }

  
}

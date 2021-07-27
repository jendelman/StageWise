#' Prepare data for BLUP 
#' 
#' Prepare data for BLUP
#' 
#' The argument \code{mask} can be used to mask all observations for a particular individual, or by including a column named 'env', only the observations in particular environments can be masked. This is useful for cross-validation to test the accuracy of predicting into new environments. 
#' 
#' @param data data frame of BLUEs from Stage 1 
#' @param vcov list of variance-covariance matrices for the BLUEs
#' @param geno object of \code{\link{class_geno}} from \code{\link{read_geno}}
#' @param vars object of \code{\link{class_var}} from \code{\link{Stage2}}
#' @param mask (optional) data frame with column "id" and optional columns "env", "trait" 
#' 
#' @return Object of \code{\link{class_prep}}
#' 
#' @import methods
#' @import Matrix
#' @importFrom stats model.matrix
#' @export

blup_prep <- function(data,vcov=NULL,geno=NULL,vars,mask=NULL) {
  
  stopifnot(inherits(data,"data.frame"))
  stopifnot(inherits(vars,"class_var"))
  data$id <- as.character(data$id)
  data$env <- as.character(data$env)
  
  if (!is.null(mask)) {
    stopifnot(is.element("id",colnames(mask)))
    ix <- which(data$id %in% as.character(mask$id))
    if ("env" %in% colnames(mask)) {
      ix2 <- which(data$env %in% as.character(mask$env))
      ix <- intersect(ix,ix2)
    }
    if ("trait" %in% colnames(mask)) {
      ix2 <- which(data$trait %in% as.character(mask$trait))
      ix <- intersect(ix,ix2)
    }
    data <- data[-ix,]
  }
  
  if (!is.null(geno)) {
    stopifnot(inherits(geno,"class_geno"))
    stopifnot(!is.na(vars@meanG))
    id <- intersect(data$id,rownames(geno@G))
    data <- data[data$id %in% id,]
    id <- rownames(geno@G)
  } else {
    id <- unique(data$id)
  }
  
  n.env <- length(unique(data$env))
  
  if ("loc" %in% colnames(data)) {
    data$loc <- as.character(data$loc)
    loc.env <- unique(data[,c("loc","env")])
    Rvec <- tapply(loc.env$loc,loc.env$env,function(z){vars@resid[z,1]})
    data$loc <- factor(data$loc)
    locations <- rownames(vars@resid)
    n.loc <- length(locations)
  } else {
    n.loc <- 1
    Rvec <- rep(vars@resid[1,1],n.env)
    loc.env <- data.frame(loc=character(0),env=character(0))
  }
  
  if (!is.null(vcov)) {
    id.ix <- lapply(vcov,function(y){which(rownames(y) %in% data$id)})
    omega.list<- mapply(FUN=function(Q,ix){as(Q[ix,ix],"dpoMatrix")},Q=vcov,ix=id.ix)
  } else {
    omega.list <- lapply(table(data$env),function(k){Matrix(0,nrow=k,ncol=k)})
  }
  
  #does not work for multi-trait yet
  tmp <- mapply(function(a,b){solve(as(a+Diagonal(nrow(a))*b,"dpoMatrix"))},omega.list,as.list(Rvec))
  Rinv <- bdiag(tmp)
  tmp <- mapply(function(a,b){as(a+Diagonal(nrow(a))*b,"dpoMatrix")},omega.list,as.list(Rvec))
  Rmat <- bdiag(tmp)

  data$id <- factor(data$id,levels=id)
  n.id <- length(id)
  data$env <- factor(data$env)
  
  n.mark <- nrow(vars@fixed.marker.var)
  if (n.mark > 0) {
    fix.eff.markers <- rownames(vars@fixed.marker.var)
    dat2 <- data.frame(id=rownames(geno@coeff),as.matrix(geno@coeff[,fix.eff.markers]))
    colnames(dat2) <- c("id",fix.eff.markers)
    data <- merge(data,dat2,by="id")
  }
  
  if ("trait" %in% colnames(data)) {
    data$trait <- factor(as.character(data$trait))
    traits <- levels(data$trait)
    n.trait <- length(traits)
    stopifnot(n.trait > 1)
    stopifnot(traits==rownames(vars@cov.trait))
    
    Z <- sparse.model.matrix(~id:trait-1,data,sep="__")
    colnames(Z) <- sub("trait__","",colnames(Z),fixed=T)
    colnames(Z) <- sub("id__","",colnames(Z),fixed=T)
    
    X <- sparse.model.matrix(~env:trait-1,data,sep="__")
    colnames(X) <- sub("trait__","",colnames(X),fixed=T)
    colnames(X) <- sub("env__","",colnames(X),fixed=T)
    
  } else {
    n.trait <- 1
    X <- sparse.model.matrix(~env-1,data,sep = "__")
    colnames(X) <- sub("env__","",colnames(X),fixed=T)
    
    if (n.loc > 1) {
      Z <- sparse.model.matrix(~id:loc-1,data,sep="__")
      colnames(Z) <- sub("loc__","",colnames(Z),fixed=T)
      colnames(Z) <- sub("id__","",colnames(Z),fixed=T)
      loc.id <- data.frame(as.matrix(expand.grid(locations,id))[,c(2,1)])
      colnames(loc.id) <- c("id","loc")
      Znames <- apply(loc.id,1,paste,collapse=":")
      Z <- Z[,Znames]
      eigen.I <- list(values=rep(1,n.id),vectors=as(Diagonal(n.id,1),"dgeMatrix"))
      dimnames(eigen.I$vectors) <- list(id,id)
      
      if (is.null(geno)) {
        tmp <- kron(eigen.A=eigen.I, B=vars@g.resid)
        Gmat.half <- tmp$half
        Gmat <- tmp$full
      } else {
        Gmat1 <- kron(eigen.A=geno@eigen.G, B=vars@add)
        Gmat2 <- kron(eigen.A=eigen.I, B=cor2cov(vars@g.resid))
        Gmat <- as(bdiag(Gmat1$full,Gmat2$full),"symmetricMatrix")
      } 
    } else {
      Z <- sparse.model.matrix(~id-1,data,sep="__")
      colnames(Z) <- sub("id__","",colnames(Z),fixed=T)
      #Z <- Z[,id]
      
      tmp <- list(half=Diagonal(n=n.id,x=sqrt(as.numeric(vars@g.resid))))
      dimnames(tmp$half) <- list(id,id)
      tmp$full <- tcrossprod(tmp$half)
    
      if (is.null(geno)) {
        Gmat.half <- tmp$half
        Gmat <- tmp$full
      } else {
        Gmat2 <- tmp
        Gmat1 <- list(half=geno@eigen.G$vectors %*% 
                        Diagonal(x=sqrt( geno@eigen.G$values * as.numeric(vars@add) )),
                      full=geno@G * as.numeric(vars@add))
        Gmat <- as(bdiag(Gmat1$full,Gmat2$full),"symmetricMatrix")
      }
    }
  }

  if (n.mark > 0) {
    if (n.trait > 1) {
      q <- paste(fix.eff.markers,"trait",sep=":")
    } else {
      if (n.loc > 1) {
        q <- paste(fix.eff.markers,"loc",sep=":")
      } else {
        q <- fix.eff.markers
      }
    }
    if (n.mark > 1) {
      q <- paste(q,collapse="+")
    } 
    q <- paste0("~",q)
    q <- paste0(q,"-1")
    tmp <- sparse.model.matrix(formula(q),data,sep="__")
    if (n.trait > 1)
      colnames(tmp) <- sub("trait__","",colnames(tmp),fixed=T)
    if (n.loc > 1) {
      colnames(tmp) <- sub("loc__","",colnames(tmp),fixed=T)
      tmp2 <- expand.grid(factor(fix.eff.markers,levels=fix.eff.markers,ordered=T),locations)
      tmp2 <- tmp2[order(tmp2$Var1),]
      marker.loc <- apply(tmp2,1,paste,collapse=":")
      colnames(tmp) <- marker.loc
    }
    X <- cbind(X,tmp)
  }
  
  n <- nrow(Z)
  m <- ncol(Z)
  if (is.null(geno)) {
    GZ <- tcrossprod(Gmat,Z)
  } else {
    m <- 2*m
    GZ <- rbind(tcrossprod(Gmat1$full,Z),tcrossprod(Gmat2$full,Z))
  } 
  
  #which way to invert V
  success <- TRUE
  if ((m < n) & (n.loc==1)) {
    cholR <- chol(Rinv)
    GZR <- tcrossprod(GZ,cholR)
    Q <- tcrossprod(GZR) + Gmat
    Qinv <- try(solve(Q),silent=TRUE)
    if (class(Qinv)!="try-error") {
      Vinv <- Rinv - crossprod(chol(Qinv)%*%GZR%*%cholR)
    } else {
      success <- FALSE
      warning("switching to direct inversion of V")
    }
  } 
     
  if ((m >= n) | !success | (n.loc > 1)) {
    if (is.null(geno)) {
      V <- tcrossprod(Z%*%Gmat.half) + Rmat
    } else {
      V <- tcrossprod(Z%*%Gmat1$half) + tcrossprod(Z%*%Gmat2$half) + Rmat
    }
    Vinv <- try(solve(V),silent=TRUE)
    if (class(Vinv)=="try-error") {
      stop("V not invertible")
    }  
  }

  W <- forceSymmetric(solve(crossprod(X,Vinv%*%X)))
  XV <- crossprod(X,Vinv)
  Pmat <- Vinv - crossprod(chol(W)%*%XV)
  fixed <- as.numeric(W%*%XV%*%data$BLUE)
  names(fixed) <- colnames(X)
  random <- as.numeric(GZ%*%Pmat%*%data$BLUE)
  new(Class="class_prep",y=data$BLUE,id=id,Z=Z,var.u=Gmat,Pmat=Pmat,Vinv=Vinv,
      fixed=fixed,random=random,add=!is.null(geno),loc.env=loc.env,
      fixed.marker=as.character(rownames(vars@fixed.marker.var)))
}    
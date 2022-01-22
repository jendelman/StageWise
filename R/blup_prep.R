#' Prepare data for BLUP 
#' 
#' Prepare data for BLUP
#' 
#' The \code{method} argument can be used to control how the linear system is solved. "MME" leads to inversion of the MME coefficient matrix, while "Vinv" leads to inversion of the overall var-cov matrix for the response vector. If NULL, the software uses whichever method involves inverting the smaller matrix. If the number of random effects (m) is less than the number of BLUEs (n), "MME" is used.
#' 
#' @param data data frame of BLUEs from Stage 1 
#' @param vcov list of variance-covariance matrices for the BLUEs
#' @param geno object of \code{\link{class_geno}} from \code{\link{read_geno}}
#' @param vars object of \code{\link{class_var}} from \code{\link{Stage2}}
#' @param mask (optional) data frame with possible columns "id","env","trait" 
#' @param method (optional) "MME", "Vinv", NULL (defaut). see Details
#' 
#' @return Object of \code{\link{class_prep}}
#' 
#' @import methods
#' @import Matrix
#' @importFrom stats model.matrix
#' @export

blup_prep <- function(data,vcov=NULL,geno=NULL,vars,mask=NULL,method=NULL) {
  
  stopifnot(inherits(data,"data.frame"))
  stopifnot(inherits(vars,"class_var"))
  data$id <- as.character(data$id)
  data$env <- as.character(data$env)
  if (!is.null(method)) {
    method <- toupper(method)
    stopifnot(method %in% c("MME","VINV"))
  }
  
  n.loc <- 1
  n.trait <- 1
  if ("loc" %in% colnames(data)) {
    locations <- rownames(vars@resid)
    n.loc <- length(locations)
    data$loc <- as.character(data$loc)
    stopifnot(n.loc > 1)
  } 
  if ("trait" %in% colnames(data)) {
    traits <- rownames(vars@resid)
    n.trait <- length(traits)
    data$trait <- as.character(data$trait)
    stopifnot(n.trait > 1)
  } 
  
  if (length(vars@meanG)>0) {
    stopifnot(!is.null(geno))
  }
  
  data <- data[!is.na(data$BLUE),]
  
  if (!is.null(mask)) {
    tmp <- intersect(colnames(mask),c("id","env","trait"))
    stopifnot(length(tmp)>0)
    if (length(tmp) > 1) {
      tmp2 <- apply(mask[,tmp],1,paste,collapse=":")
      ix <- which(apply(data[,tmp],1,paste,collapse=":") %in% tmp2)
    } else {
      tmp2 <- mask[,tmp]
      ix <- which(data[,tmp] %in% tmp2)
    }
    if (length(ix) > 0) {
      data <- data[-ix,]
    }
  }
  
  if (!is.null(geno)) {
    stopifnot(inherits(geno,"class_geno"))
    stopifnot(length(vars@meanG)>0)
    id <- intersect(data$id,rownames(geno@G))
    data <- data[data$id %in% id,]
    id <- rownames(geno@G)
  } else {
    id <- unique(data$id)
  }
  
  if (n.trait > 1) {
    tmp <- paste(data$id,data$trait,sep=":")
  } else {
    tmp <- data$id
  }
  if (!is.null(vcov)) {
    stopifnot(nrow(vars@meanOmega) > 0)
    data$env <- factor(data$env,levels=names(vcov))
    omega.list <- mapply(FUN=function(Q,ix){
                          ix2 <- match(ix,rownames(Q))
                          as(Q[ix2,ix2,drop=FALSE],"dpoMatrix")
                        },Q=vcov,ix=split(tmp,data$env))
    #names(omega.list) <- names(vcov)
  } else {
    data$env <- factor(data$env)
    omega.list <- lapply(split(tmp,data$env),function(rnames){
                                n <- length(rnames)
                                Q <- Matrix(0,nrow=n,ncol=n,dimnames=list(rnames,rnames))
                                return(Q)})
  }
  
  #redo envs because some may have been dropped
  data$env <- as.character(data$env)
  loc.env <- data.frame(loc=character(0),env=character(0))
  trait.env <- data.frame(trait=character(0),env=character(0))
  if (n.loc > 1) 
    loc.env <- unique(data[,c("loc","env")])
  if (n.trait > 1)
    trait.env <- unique(data[,c("trait","env")])
  
  envs <- unique(data$env)
  n.env <- length(envs)
  omega.list <- omega.list[envs]
  data$env <- factor(data$env,levels=envs)
  n.obs <- sapply(omega.list,nrow)
  
  #Rlist
  if (n.trait==1) {
    if (n.loc > 1) {
      tmp <- vars@resid[as.character(data$loc)[match(envs,as.character(data$env))],1]
    } else {
      tmp <- rep(vars@resid[1,1],n.env)
    }
    Rlist <- mapply(FUN=function(n,v){Diagonal(n=n)*v},n=as.list(n.obs),v=as.list(tmp))
  } else {
    #multi-trait
    tmp <- split(data$id,data$env)
    tmp2 <- lapply(tmp,function(id) {
      n.id <- length(id)
      eigen.I <- list(values=rep(1,n.id),vectors=as(Diagonal(n.id,1),"dgeMatrix"))
      dimnames(eigen.I$vectors) <- list(id,id)
      crossprod(kron(eigen.I,vars@resid)$mat)
    })
    Rlist <- mapply(function(Q,rnames){
                    Q[rnames,rnames]},
                    Q=tmp2,rnames=lapply(omega.list,rownames))
  }

  n.id <- length(id)
  data$id <- factor(data$id,levels=id)
  
  if (n.trait > 1)
    data$trait <- factor(data$trait,levels=traits)
  if (n.loc > 1)
    data$loc <- factor(data$loc,levels=locations)
  
  n.mark <- dim(vars@fixed.marker.cov)[3]
  fix.eff.markers <- dimnames(vars@fixed.marker.cov)[[3]]
  if (n.mark > 0) {
    dat2 <- data.frame(id=rownames(geno@coeff),as.matrix(geno@coeff[,fix.eff.markers]))
    colnames(dat2) <- c("id",fix.eff.markers)
    data <- merge(data,dat2,by="id")
  }
  
  if (n.trait==1) {
    if (n.env > 1) {
      X <- sparse.model.matrix(~env-1,data,sep = "__")
      colnames(X) <- sub("env__","",colnames(X),fixed=T)
    } else {
      X <- sparse.model.matrix(~1,data)
      colnames(X) <- envs
    }
    
    if (n.loc > 1) {
      Z <- sparse.model.matrix(~id:loc-1,data,sep="__")
      colnames(Z) <- sub("loc__","",colnames(Z),fixed=T)
      colnames(Z) <- sub("id__","",colnames(Z),fixed=T)
      loc.id <- data.frame(as.matrix(expand.grid(locations,id))[,c(2,1)])
      colnames(loc.id) <- c("id","loc")
      Znames <- apply(loc.id,1,paste,collapse=":")
      Z <- Z[,Znames] #loc within id
      eigen.I <- list(values=rep(1,n.id),vectors=as(Diagonal(n.id,1),"dgeMatrix"))
      dimnames(eigen.I$vectors) <- list(id,id)
      
      if (is.null(geno)) {
        Gmat <- kron(eigen.A=eigen.I, B=vars@g.resid)
        index.scale <- sqrt(diag(as.matrix(vars@g.resid)))
      } else {
        Gmat1 <- kron(eigen.A=geno@eigen.G, B=vars@add)
        Gmat2 <- kron(eigen.A=eigen.I, B=coerce_dpo(tcrossprod(sqrt(vars@g.resid))))
        Gmat <- list(mat=bdiag(Gmat1$mat,Gmat2$mat),inv=bdiag(Gmat1$inv,Gmat2$inv))
        index.scale <- sqrt(diag(as.matrix(vars@add)))
      } 
    } else {
      index.scale <- numeric(0)
      Z <- sparse.model.matrix(~id-1,data,sep="__")
      colnames(Z) <- sub("id__","",colnames(Z),fixed=T)

      tmp <- list(mat=Diagonal(n=n.id,x=sqrt(as.numeric(vars@g.resid))),
                  inv=Diagonal(n=n.id,x=1/sqrt(as.numeric(vars@g.resid))))
      dimnames(tmp$inv) <- dimnames(tmp$mat) <- list(id,id)
    
      if (is.null(geno)) {
        Gmat <- tmp
      } else {
        Gmat1 <- list(mat=tcrossprod(Diagonal(x=sqrt(geno@eigen.G$values * as.numeric(vars@add))),
                                     geno@eigen.G$vectors),
                      inv=tcrossprod(Diagonal(x=1/sqrt(geno@eigen.G$values * as.numeric(vars@add))),
                                     geno@eigen.G$vectors))
        Gmat <- list(mat=bdiag(Gmat1$mat,tmp$mat),inv=bdiag(Gmat1$inv,tmp$inv))
      }
    }
  } else {
    #multi-trait
    Z <- sparse.model.matrix(~id:trait-1,data,sep="__")
    colnames(Z) <- sub("trait__","",colnames(Z),fixed=T)
    colnames(Z) <- sub("id__","",colnames(Z),fixed=T)
    
    if (n.env > 1) {
      X <- sparse.model.matrix(~env:trait-1,data,sep="__")
      colnames(X) <- sub("trait__","",colnames(X),fixed=T)
      colnames(X) <- sub("env__","",colnames(X),fixed=T)
    } else {
      X <- sparse.model.matrix(~trait-1,data,sep="__")
      colnames(X) <- sub("trait__","",colnames(X),fixed=T)
      colnames(X) <- paste(envs,colnames(X),sep=":")
    }
    
    trait.id <- data.frame(as.matrix(expand.grid(traits,id))[,c(2,1)])
    colnames(trait.id) <- c("id","trait")
    Znames <- apply(trait.id,1,paste,collapse=":")
    Z <- Z[,Znames] #trait within id
    eigen.I <- list(values=rep(1,n.id),vectors=as(Diagonal(n.id,1),"dgeMatrix"))
    dimnames(eigen.I$vectors) <- list(id,id)
    
    if (is.null(geno)) {
      Gmat <- kron(eigen.A=eigen.I, B=vars@g.resid)
      index.scale <- sqrt(diag(as.matrix(vars@g.resid)))
    } else {
      Gmat1 <- kron(eigen.A=geno@eigen.G, B=vars@add)
      Gmat2 <- kron(eigen.A=eigen.I, B=vars@g.resid)
      Gmat <- list(mat=bdiag(Gmat1$mat,Gmat2$mat),inv=bdiag(Gmat1$inv,Gmat2$inv))
      index.scale <- sqrt(diag(as.matrix(vars@add)))
    } 
  }
  
  if (!is.null(geno)) {
    Z <- cbind(Z,Z)
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
    if (n.trait > 1) {
      colnames(tmp) <- sub("trait__","",colnames(tmp),fixed=T)
      tmp2 <- expand.grid(factor(fix.eff.markers,levels=fix.eff.markers,ordered=T),traits)
      tmp2 <- tmp2[order(tmp2$Var1),]
      marker.trait <- apply(tmp2,1,paste,collapse=":")
      colnames(tmp) <- marker.trait
    }
    if (n.loc > 1) {
      colnames(tmp) <- sub("loc__","",colnames(tmp),fixed=T)
      tmp2 <- expand.grid(factor(fix.eff.markers,levels=fix.eff.markers,ordered=T),locations)
      tmp2 <- tmp2[order(tmp2$Var1),]
      marker.loc <- apply(tmp2,1,paste,collapse=":")
      colnames(tmp) <- marker.loc
    }
    X <- cbind(X,tmp)
  }
  
  m <- ncol(Z)
  n <- nrow(Z)
  var.u <- crossprod(Gmat$mat)
  
  if (is.null(method)) {
    if (m < n) {
      method <- "MME"
    } else {
      method <- "VINV"
    }
  }
  if (method=="MME") {
    #Construct MME coefficient matrix
    tmp <- mapply(function(a,b){chol(solve(as(a+b,"dpoMatrix")))},omega.list,Rlist)
    Rinv <- bdiag(tmp)
  
    n.fix <- ncol(X)
    RZ <- Rinv %*% Z
    RX <- Rinv %*% X
    Q <- cbind(RX, RZ)
  
    MME <- as(crossprod(Q) + crossprod(cbind(Matrix(0,ncol=n.fix,nrow=m),Gmat$inv)),"symmetricMatrix")
    MME.inv <- as(solve(MME),"symmetricMatrix")
    soln <- MME.inv %*% crossprod(Q,Rinv%*%data$BLUE)
    fixed <- as.numeric(soln[1:n.fix])
    names(fixed) <- colnames(X)
    random.ix <- (n.fix+1):length(soln)
    random <- as.numeric(soln[random.ix])
    var.uhat=var.u - MME.inv[random.ix,random.ix]
    
  } else {
    #invert V
    tmp <- mapply(function(a,b){as(a+b,"dpoMatrix")},omega.list,Rlist)
    Rmat <- bdiag(tmp)
    
    Vinv <- as(solve(as(tcrossprod(Z %*% t(Gmat$mat)) + Rmat,"symmetricMatrix")),"symmetricMatrix")
    chol.Vinv <- chol(Vinv)
    
    W <- as(solve(crossprod(chol.Vinv%*%X)),"symmetricMatrix")
    fixed <- as.numeric(tcrossprod(W,X) %*% Vinv %*% data$BLUE)
    names(fixed) <- colnames(X)
        
    WW <- as(Diagonal(n=n) - chol.Vinv %*% X %*% tcrossprod(W, chol.Vinv %*% X),"symmetricMatrix")
    cholWW <- suppressWarnings(try(chol(WW),silent=TRUE))
    if (class(cholWW)=="try-error") {
      cholWW <- chol(WW + Diagonal(n=n,x=1e-6))
    }
    GZtP <- tcrossprod(var.u,Z) %*% crossprod(cholWW %*% chol.Vinv)
    random <- as.numeric(GZtP %*% data$BLUE)
    var.uhat <- GZtP %*% tcrossprod(Z,var.u)
  }
  
  new(Class="class_prep",id=id,var.u=var.u,var.u.inv=crossprod(Gmat$inv),var.uhat=var.uhat,
      fixed=fixed,random=random,add=vars@add,loc.env=loc.env,trait.env=trait.env,
      fixed.marker=as.character(fix.eff.markers),index.scale=index.scale)
  
}    

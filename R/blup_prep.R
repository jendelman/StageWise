#' Prepare data for BLUP 
#' 
#' Prepare data for BLUP
#' 
#' @param data data frame of BLUEs from Stage 1 
#' @param vcov list of variance-covariance matrices for the BLUEs
#' @param geno object of \code{\link{class_geno}} from \code{\link{read_geno}}
#' @param vars object of \code{\link{class_var}} from \code{\link{Stage2}}
#' @param mask (optional) data frame with possible columns "id","env","trait" 
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
  
  if (!is.na(vars@meanG)) {
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
    stopifnot(!is.na(vars@meanG))
    id <- intersect(data$id,rownames(geno@G))
    data <- data[data$id %in% id,]
    id <- rownames(geno@G)
  } else {
    id <- unique(data$id)
  }
  
  data$env <- factor(data$env,levels=names(vcov))
  
  if (n.trait > 1) {
    tmp <- paste(data$id,data$trait,sep=":")
  } else {
    tmp <- data$id
  }
  if (!is.null(vcov)) {
    omega.list <- mapply(FUN=function(Q,ix){
                          ix2 <- match(ix,rownames(Q))
                          as(Q[ix2,ix2,drop=FALSE],"dpoMatrix")
                        },Q=vcov,ix=split(tmp,data$env))
  } else {
    omega.list <- lapply(split(tmp,data$env),function(rnames){
                                n <- length(rnames)
                                Q <- Matrix(0,nrow=n,ncol=n,dimnames=list(rnames,rnames))
                                return(Q)})
  }
  names(omega.list) <- names(vcov)
  
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
      kron(eigen.I,vars@resid)$full
    })
    Rlist <- mapply(function(Q,rnames){
                    as(Q[rnames,rnames],"dpoMatrix")},
                    Q=tmp2,rnames=lapply(omega.list,rownames))
  }

  tmp <- mapply(function(a,b){solve(as(a+b,"dpoMatrix"))},omega.list,Rlist)
  Rinv <- bdiag(tmp)
  tmp <- mapply(function(a,b){as(a+b,"dpoMatrix")},omega.list,Rlist)
  Rmat <- bdiag(tmp)

  n.id <- length(id)
  data$id <- factor(data$id,levels=id)
  
  if (n.trait > 1)
    data$trait <- factor(data$trait,levels=traits)
  if (n.loc > 1)
    data$loc <- factor(data$loc,levels=locations)
  
  n.mark <- nrow(vars@fixed.marker.var)
  if (n.mark > 0) {
    fix.eff.markers <- rownames(vars@fixed.marker.var)
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
        tmp <- kron(eigen.A=eigen.I, B=vars@g.resid)
        Gmat.half <- tmp$half
        Gmat <- tmp$full
        index.scale <- sqrt(diag(as.matrix(vars@g.resid)))
      } else {
        Gmat1 <- kron(eigen.A=geno@eigen.G, B=vars@add)
        Gmat2 <- kron(eigen.A=eigen.I, B=coerce_dpo(tcrossprod(sqrt(vars@g.resid))))
        Gmat <- as(bdiag(Gmat1$full,Gmat2$full),"symmetricMatrix")
        index.scale <- sqrt(diag(as.matrix(vars@add)))
      } 
    } else {
      index.scale <- numeric(0)
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
      tmp <- kron(eigen.A=eigen.I, B=vars@g.resid)
      Gmat.half <- tmp$half
      Gmat <- tmp$full
      index.scale <- sqrt(diag(as.matrix(vars@g.resid)))
    } else {
      Gmat1 <- kron(eigen.A=geno@eigen.G, B=vars@add)
      Gmat2 <- kron(eigen.A=eigen.I, B=vars@g.resid)
      Gmat <- as(bdiag(Gmat1$full,Gmat2$full),"symmetricMatrix")
      index.scale <- sqrt(diag(as.matrix(vars@add)))
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
  if ((m < n) & (n.loc==1) & (n.trait==1)) {
    cholR <- chol(Rinv)
    GZR <- tcrossprod(GZ,cholR)
    Q <- tcrossprod(GZR) + Gmat
    Qinv <- try(solve(Q),silent=TRUE)
    if (class(Qinv)!="try-error") {
      Vinv <- coerce_dpo(Rinv - crossprod(chol(Qinv)%*%GZR%*%cholR))
    } else {
      success <- FALSE
      warning("switching to direct inversion of V")
    }
  }
     
  if ((m >= n) | !success | (n.loc > 1) | (n.trait > 1)) {
    if (is.null(geno)) {
      V <- tcrossprod(Z%*%Gmat.half) + Rmat
    } else {
      V <- tcrossprod(Z%*%Gmat1$half) + tcrossprod(Z%*%Gmat2$half) + Rmat
    }
    V <- coerce_dpo(V)
    Vinv <- try(solve(V),silent=TRUE)
    if (class(Vinv)=="try-error") {
      stop("V not invertible")
    }  
  }

  W <- solve(crossprod(chol(Vinv)%*%X))
  XV <- crossprod(X,Vinv)
  Pmat <- coerce_dpo(Vinv - crossprod(chol(W)%*%XV))
  fixed <- as.numeric(W%*%XV%*%data$BLUE)
  names(fixed) <- colnames(X)
  random <- as.numeric(GZ%*%Pmat%*%data$BLUE)
  new(Class="class_prep",y=data$BLUE,id=id,Z=Z,var.u=Gmat,Pmat=Pmat,Vinv=Vinv,
      fixed=fixed,random=random,add=vars@add,loc.env=loc.env,trait.env=trait.env,
      fixed.marker=as.character(rownames(vars@fixed.marker.var)),
      index.scale=index.scale)
}    

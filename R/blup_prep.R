#' Prepare data for BLUP 
#' 
#' Prepare data for BLUP
#' 
#' The \code{method} argument can be used to control how the linear system is solved. "MME" leads to inversion of the MME coefficient matrix, while "Vinv" leads to inversion of the overall var-cov matrix for the response vector. If NULL, the software uses whichever method involves inverting the smaller matrix. If the number of random effects (m) is less than the number of BLUEs (n), "MME" is used.
#' 
#' For the multi-location model, if all of the environments for a location are masked, the average of the other locations is used when computed average fixed effects.
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
#' @importFrom MASS ginv
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

  if (vars@model==3L)
    stopifnot(class(geno)=="class_genoD")
  
  n.trait <- 1
  if ("trait" %in% colnames(data)) {
    traits <- rownames(vars@resid)
    n.trait <- length(traits)
    data$trait <- as.character(data$trait)
    stopifnot(n.trait > 1)
  } 
  
  if (length(vars@diagG)>0) {
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
    stopifnot(length(vars@diagG)>0)
    id <- intersect(data$id,rownames(geno@G))
    data <- data[data$id %in% id,]
    id <- rownames(geno@G)
    ploidy <- geno@ploidy
  } else {
    id <- unique(data$id)
    ploidy <- 0L
  }
  
  if (n.trait > 1) {
    data <- data[order(data$env,data$id,data$trait),]
    tmp <- paste(data$id,data$trait,sep=":")
  } else {
    data <- data[order(data$env,data$id),]
    tmp <- data$id
  }
  if (!is.null(vcov)) {
    #stopifnot(nrow(vars@meanOmega) > 0)
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
  envs <- unique(data$env)
  n.env <- length(envs)
  omega.list <- omega.list[envs]
  data$env <- factor(data$env,levels=envs)
  n.obs <- sapply(omega.list,nrow)
  
  n.loc <- 1
  if ("loc" %in% colnames(data)) {
    locations <- rownames(vars@resid)
    n.loc <- length(locations)
    data$loc <- as.character(data$loc)
    stopifnot(n.loc > 1)
    missing.loc <- setdiff(locations,as.character(data$loc))
  } else {
    missing.loc <- character(0)
  }
  
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
      #eigen.I <- list(values=rep(1,n.id),vectors=as(Diagonal(n.id,1),"dgeMatrix"))
      eigen.I <- list(values=rep(1,n.id),vectors=as(as(Diagonal(n.id,1),"generalMatrix"),"unpackedMatrix"))
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
  
  n.mark <- length(vars@fix.eff.marker)
  if (n.mark > 0) {
    dat2 <- data.frame(id=rownames(geno@coeff),as.matrix(geno@coeff[,vars@fix.eff.marker]))
    colnames(dat2) <- c("id",vars@fix.eff.marker)
    data <- merge(data,dat2,by="id")
  }
  
  eigen.I <- list(values=rep(1,n.id),vectors=as(as(Diagonal(n.id,1),"generalMatrix"),"unpackedMatrix"))
  dimnames(eigen.I$vectors) <- list(id,id)
  
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
      
      if (is.null(geno)) {
        Gmat <- kron(eigen.A=eigen.I, B=vars@geno1)
      } else {
        if (vars@model==1L) {
          Gmat <- kron(eigen.A=geno@eigen.G, B=vars@geno1)
        } else {
          Gmat1 <- kron(eigen.A=geno@eigen.G, B=vars@geno1)
          if (vars@model==3L) {
            Gmat2 <- kron(eigen.A=geno@eigen.D, B=vars@geno2)
          
            #add to X matrix
            dom.covariate <- Z %*% kronecker(geno@coeff.D %*% matrix(1,nrow=ncol(geno@coeff.D),ncol=1),
                                                      diag(n.loc))
            colnames(dom.covariate) <- paste("heterosis",locations,sep=":")
            X <- cbind(X,dom.covariate/(geno@scale*(ploidy-1)))
          } else {
            Gmat2 <- kron(eigen.A=eigen.I, B=vars@geno2)
          }
          Gmat <- list(mat=bdiag(Gmat1$mat,Gmat2$mat),inv=bdiag(Gmat1$inv,Gmat2$inv))
        }
      } 
    } else {
      
      Z <- sparse.model.matrix(~id-1,data,sep="__")
      colnames(Z) <- sub("id__","",colnames(Z),fixed=T)

      if (is.null(geno)) {
        Gmat <- kron(eigen.A = eigen.I, B=vars@geno1)
      } else {
        if (vars@model==1L) {
          Gmat <- kron(eigen.A = geno@eigen.G, B=vars@geno1)
        } else {
          Gmat1 <- kron(eigen.A = geno@eigen.G, B=vars@geno1)
          if (vars@model==3L) {
            Gmat2 <- kron(eigen.A=geno@eigen.D, B=vars@geno2)
            dom.covariate <- as.numeric(Z %*% geno@coeff.D %*% matrix(1,nrow=ncol(geno@coeff.D),ncol=1))
            X <- cbind(X,heterosis=dom.covariate/(geno@scale*(ploidy-1)))
          } else {
            Gmat2 <- kron(eigen.A=eigen.I, B=vars@geno2)
          }
          Gmat <- list(mat=bdiag(Gmat1$mat,Gmat2$mat),inv=bdiag(Gmat1$inv,Gmat2$inv))
        }
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

    if (is.null(geno)) {
      Gmat <- kron(eigen.A=eigen.I, B=vars@geno1)
    } else {
      if (vars@model==1L) {
        Gmat <- kron(eigen.A=geno@eigen.G, B=vars@geno1)
      } else {
        Gmat1 <- kron(eigen.A=geno@eigen.G, B=vars@geno1)
        if (vars@model==3L) {
          Gmat2 <- kron(eigen.A=geno@eigen.D, B=vars@geno2)
        
          #add to X matrix
          dom.covariate <- Z %*% kronecker(geno@coeff.D %*% matrix(1,nrow=ncol(geno@coeff.D),ncol=1),
                                         diag(n.trait))
          colnames(dom.covariate) <- paste("heterosis",traits,sep=":")
          X <- cbind(X,dom.covariate/(geno@scale*(ploidy-1)))
        } else {
          Gmat2 <- kron(eigen.A=eigen.I, B=vars@geno2)
        }
        Gmat <- list(mat=bdiag(Gmat1$mat,Gmat2$mat),inv=bdiag(Gmat1$inv,Gmat2$inv))
      }
    } 
  }
  
  if (vars@model > 1L) {
    Z <- cbind(Z,Z)
  }

  if (n.mark > 0) {
    if (n.trait > 1) {
      q <- paste(vars@fix.eff.marker,"trait",sep=":")
    } else {
      if (n.loc > 1) {
        q <- paste(vars@fix.eff.marker,"loc",sep=":")
      } else {
        q <- vars@fix.eff.marker
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
      tmp2 <- expand.grid(factor(vars@fix.eff.marker,levels=vars@fix.eff.marker,ordered=T),traits)
      tmp2 <- tmp2[order(tmp2$Var1),]
      marker.trait <- apply(tmp2,1,paste,collapse=":")
      colnames(tmp) <- marker.trait
    }
    if (n.loc > 1) {
      colnames(tmp) <- sub("loc__","",colnames(tmp),fixed=T)
      tmp2 <- expand.grid(factor(vars@fix.eff.marker,levels=vars@fix.eff.marker,ordered=T),locations)
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
    var.uhat <- var.u - MME.inv[random.ix,random.ix]
    
  } else {
    #invert V
    tmp <- mapply(function(a,b){as(a+b,"dpoMatrix")},omega.list,Rlist)
    Rmat <- bdiag(tmp)
    
    Vinv <- as(solve(as(tcrossprod(Z %*% t(Gmat$mat)) + Rmat,"symmetricMatrix")),"symmetricMatrix")
    chol.Vinv <- chol(Vinv)
    
    tmp <- crossprod(chol.Vinv%*%X)
    tmp2 <- try(solve(tmp),silent=TRUE)
    if (class(tmp2)=="try-error")
      tmp2 <- MASS::ginv(as.matrix(tmp))
    W <- forceSymmetric(tmp2)
    fixed <- as.numeric(tcrossprod(W,X) %*% Vinv %*% data$BLUE)
    names(fixed) <- colnames(X)
    
    tmp <- Diagonal(n=n) - chol.Vinv %*% X %*% tcrossprod(W, chol.Vinv %*% X)
    WW <- forceSymmetric(tmp)
    cholWW <- suppressWarnings(try(chol(WW),silent=TRUE))
    if (class(cholWW)=="try-error") {
      cholWW <- chol(WW + Diagonal(n=n,x=1e-6))
    }
    GZtP <- tcrossprod(var.u,Z) %*% crossprod(cholWW %*% chol.Vinv)
    random <- as.numeric(GZtP %*% data$BLUE)
    var.uhat <- GZtP %*% tcrossprod(Z,var.u)
  }
  
  fixed.marker <- numeric(0)
  heterosis <- numeric(0)
  
  if (n.loc > 1) {
    loc.env <- unique(data[,c("loc","env")])
    loc.env <- loc.env[order(loc.env$loc,loc.env$env),]
    ix <- match(names(fixed)[1:n.env],as.character(loc.env$env))
    avg.env <- tapply(fixed[1:n.env],loc.env$loc[ix],mean)
    tmp <- as.numeric(avg.env)
    names(tmp) <- names(avg.env)
    avg.env <- tmp
  } else {
    if (n.trait > 1) {
      env.trait <- expand.grid(env=envs,trait=traits)
      env.trait$et <- apply(env.trait[,1:2],1,paste,collapse=":")
      ix <- match(env.trait$et,names(fixed))
      avg.env <- tapply(fixed[ix],env.trait$trait,mean)
      tmp <- as.numeric(avg.env)
      names(tmp) <- names(avg.env)
      avg.env <- tmp
      n.env <- n.env*n.trait
    } else {
      avg.env <- mean(fixed[1:n.env])
    }
  }
  
  nlt <- max(n.loc,n.trait)
  if (vars@model==3L) {
    heterosis <- fixed[n.env+1:nlt]
    if (length(missing.loc)>0)
      heterosis[paste("heterosis",missing.loc,sep=":")] <- mean(heterosis,na.rm=T)
    if (n.mark > 0) 
      fixed.marker <- fixed[n.env+nlt+1:(n.mark*nlt)]
  } else {
    if (n.mark > 0) 
      fixed.marker <- fixed[n.env+1:(n.mark*nlt)]
  }
  if (length(missing.loc)>0) {
    avg.env[missing.loc] <- mean(avg.env,na.rm=T)
  }
  
  new(Class="class_prep",id=id,ploidy=ploidy,var.u=var.u,var.uhat=var.uhat,
      avg.env=avg.env,heterosis=heterosis,fixed.marker=fixed.marker,B=vars@B,
      random=random, geno1.var=vars@geno1, geno2.var=vars@geno2, 
      model=vars@model)
  
}    

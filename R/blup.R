#' BLUPs from Stage2
#' 
#' BLUPs from Stage2
#' 
#' Index coefficients in \code{index.coeff} are normalized to sum to 1. When \code{gwas=TRUE}, GWAS -log10(p) scores are computed by standardizing the marker effects. When \code{geno} is not NULL, both breeding value (BV) and (GV) values for the individuals are predicted; otherwise, only GV are returned.
#' 
#' @param data data frame of BLUEs from Stage 1 
#' @param vcov list of variance-covariance matrices for the BLUEs
#' @param geno object of \code{\link{class_geno}} from \code{\link{read_geno}}
#' @param vars object of \code{\link{class_var}} from \code{\link{Stage2}}
#' @param mask (optional) data frame with column "id" and optional columns "env", "trait" 
#' @param reliability TRUE/FALSE whether to compute reliability
#' @param index.coeff index coefficients for the locations or traits (must sum to 1)
#' @param marker.effects TRUE/FALSE
#' @param gwas.ncore Number of cores to use for GWAS (default is 0 = no GWAS)
#' 
#' @return Data frame containing the BLUPs for the individuals, unless marker.effects = TRUE, in which case a list of two data frames: one for the individuals and one for the markers. 
#' 
#' @import methods
#' @import Matrix
#' @importFrom stats formula pnorm
#' @importFrom stats model.matrix
#' @importFrom parallel makeCluster clusterExport parCapply stopCluster
#' @export

blup <- function(data,vcov=NULL,geno=NULL,vars,mask=NULL,reliability=FALSE,
                 index.coeff=NULL,marker.effects=FALSE,gwas.ncore=0) {
  
  stopifnot(inherits(data,"data.frame"))
  stopifnot(inherits(vars,"class_var"))
  data$id <- as.character(data$id)
  data$env <- as.character(data$env)
  
  if (!is.null(mask)) {
    stopifnot(is.element("id",colnames(mask)))
    ix <- which(as.character(data$id) %in% as.character(mask$id))
    if ("env" %in% colnames(mask)) {
      ix2 <- which(as.character(data@data$env) %in% as.character(mask$env))
      ix <- intersect(ix,ix2)
    }
    if ("trait" %in% colnames(mask)) {
      ix2 <- which(as.character(data@data$trait) %in% as.character(mask$trait))
      ix <- intersect(ix,ix2)
    }
    data <- data[-ix,]
  }
  
  if (!is.null(geno)) {
    stopifnot(inherits(geno,"class_geno"))
    stopifnot(!is.na(vars@meanG))
    id <- intersect(data$id,rownames(geno@Ginv))
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
    if (is.null(index.coeff)) {
      index.coeff <- rep(1/n.loc,n.loc)
      names(index.coeff) <- locations
    } else {
      index.coeff <- index.coeff/sum(index.coeff)
      ix <- match(locations,names(index.coeff))
      stopifnot(!is.na(ix))
      index.coeff <- index.coeff[ix]
    }
  } else {
    n.loc <- 1
    index.coeff <- 1
    Rvec <- rep(vars@resid[1,1],n.env)
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
      Imat <- Diagonal(n.id,1)
      dimnames(Imat) <- list(id,id)
      
      if (is.null(geno)) {
        Ginv <- kronecker(Imat,solve(vars@g.resid),make.dimnames = T)
        dimnames(Ginv) <- list(Znames,Znames)
        Gmat <- kronecker(Imat,vars@g.resid,make.dimnames = T)
        dimnames(Gmat) <- list(Znames,Znames)
      } else {
        Ginv1 <- kronecker(geno@Ginv,solve(vars@add),make.dimnames = T)
        Gloc <- cor2cov(vars@g.resid)
        Ginv2 <- kronecker(Imat,solve(Gloc),make.dimnames = TRUE)
        dimnames(Ginv2) <- list(Znames,Znames)
        Ginv <- bdiag(Ginv1,Ginv2)
        Gmat1 <- kronecker(geno@G,vars@add,make.dimnames = T)
        Gloc <- cor2cov(vars@g.resid)
        Gmat2 <- kronecker(Imat,Gloc,make.dimnames = TRUE)
        dimnames(Gmat2) <- list(Znames,Znames)
        Gmat <- bdiag(Gmat1,Gmat2)
      } 
    } else {
      
      Z <- sparse.model.matrix(~id-1,data,sep="__")
      colnames(Z) <- sub("id__","",colnames(Z),fixed=T)
      Z <- Z[,id]
      
      if (is.null(geno)) {
        Ginv <- Diagonal(n=n.id,x=1/as.numeric(vars@g.resid))
        dimnames(Ginv) <- list(id,id)
        Gmat <- Diagonal(n=n.id,x=as.numeric(vars@g.resid))
        dimnames(Ginv) <- list(id,id)
      } else {
        Ginv1 <- geno@Ginv/as.numeric(vars@add)
        Ginv2 <- Diagonal(n=n.id,x=1/as.numeric(vars@g.resid))
        dimnames(Ginv2) <- list(id,id)
        Ginv <- bdiag(Ginv1,Ginv2)
        Gmat1 <- geno@G*as.numeric(vars@add)
        Gmat2 <- Diagonal(n=n.id,x=as.numeric(vars@g.resid))
        dimnames(Gmat2) <- list(id,id)
        Gmat <- bdiag(Gmat1,Gmat2)
      } 
    }
  }

  if (is.null(geno)) {
    Z2 <- Z
  } else {
    Z2 <- cbind(Z,Z)
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
  
  #MME
  XtR <- crossprod(X,Rinv)
  ZtR <- crossprod(Z2,Rinv)
  Q11 <- XtR%*%X
  Q12 <- XtR%*%Z2
  Q22 <- ZtR%*%Z2 + Ginv
  Q <- rbind(cbind(Q11,Q12),cbind(t(Q12),Q22))
  b <- c(as.numeric(XtR%*%data$BLUE),as.numeric(ZtR%*%data$BLUE))
  
  n.fix <- ncol(X)
  if (reliability) {
    Cmat <- solve(Q)
    soln <- Cmat %*% b
    C22 <- Cmat[-(1:n.fix),-(1:n.fix)]
  } else {
    soln <- solve(Q,b)
  }

  fix.value <- mean(soln[1:n.env])
  
  if (n.mark > 0) {
    a <- kronecker(Diagonal(n=n.mark),Matrix(index.coeff,nrow=1)) %*% Matrix(soln[n.env+1:(n.loc*n.mark)],ncol=1)
    fix.value <- fix.value + as.numeric(as.matrix(geno@coeff[,fix.eff.markers]) %*% a)
  } 
    
  random <- soln[-(1:n.fix)]
  M <- kronecker(Diagonal(n=n.id),Matrix(index.coeff,nrow=1))

  if (n.trait==1) {
    if (is.null(geno)) {
      out <- data.frame(id=id,GV=as.numeric(M%*%Matrix(random,ncol=1)) + fix.value)
    } else {
      out <- data.frame(id=id,BV=as.numeric(cbind(M,0*M) %*% Matrix(random,ncol=1)) + fix.value,
                        GV=as.numeric(cbind(M,M) %*% Matrix(random,ncol=1)) + fix.value)
    }
  } else {
    #multi-trait
  }
  
  if (reliability) {
    if (is.null(geno)) {
      numer <- diag(M%*%C22%*%t(M))
      denom <- diag(M%*%Gmat%*%t(M))
      out$GV.r2 <- 1-numer/denom
    } else {
      numer <- diag(cbind(M,0*M)%*%C22%*%t(cbind(M,0*M)))
      denom <- diag(cbind(M,0*M)%*%Gmat%*%t(cbind(M,0*M)))
      out$BV.r2 <- 1-numer/denom
      numer <- diag(cbind(M,M)%*%C22%*%t(cbind(M,M)))
      denom <- diag(cbind(M,M)%*%Gmat%*%t(cbind(M,M)))
      out$GV.r2 <- 1-numer/denom
    }
  }
  
  if (!marker.effects) {
    return(out)
  } else {
    out2 <- vector("list",2)
    out2[[1]] <- out
    V <- as(Z2 %*% tcrossprod(Gmat,Z2) + Rmat,"dpoMatrix")
    y2 <- data$BLUE - X %*% Matrix(soln[1:n.fix],ncol=1)
    m <- ncol(geno@coeff)
    if (n.loc > 1) {
      Ct <- Z %*% kronecker(geno@coeff,Diagonal(n=n.loc,x=1)) %*% 
        kronecker(Matrix(1,nrow=m,ncol=1),Matrix(index.coeff,ncol=1))
    } else {
      Ct <- Z%*%geno@coeff
    }
    out2[[2]] <- data.frame(marker=colnames(geno@coeff),
                              add.effect=as.numeric(crossprod(Ct,solve(V,y2))))
    
    if (gwas.ncore > 0) {
      Vinv <- solve(V)
      Q <- tcrossprod(Vinv%*%X%*%solve(crossprod(X,Vinv%*%X)),Vinv%*%X)
      
      if (gwas.ncore > 1) {
        cl <- makeCluster(gwas.ncore)
        clusterExport(cl=cl,varlist=NULL)
        se <- parCapply(cl=cl,x=Ct,function(z){
          sqrt(crossprod(z,Vinv%*%z) - crossprod(z,Q%*%z))
        })
        std.effect <- out2[[2]]$add.effect/sapply(se,as.numeric)
        stopCluster(cl)
      } else {
        se <- apply(Ct,2,function(z){
          sqrt(crossprod(z,Vinv%*%z) - crossprod(z,Q%*%z))
        })
        std.effect <- out2[[2]]$add.effect/se
      }
      out2[[2]]$score <- -log10(pnorm(q=abs(std.effect),lower.tail=FALSE)*2)
    }
    return(out2)
  }
}

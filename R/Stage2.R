#' Stage 2 analysis of multi-environment trials
#' 
#' Stage 2 analysis of multi-environment trials 
#' 
#' Stage 2 of the two-stage approach described by Damesa et al. 2017, using ASReml-R for variance component estimation (license is required). The variable \code{data} has two mandatory column names: id = individual (genotype identifier), and env = environment at which Stage 1 analysis was performed. The argument \code{traits} is a character vector that must match column names in \code{data}. Missing data are allowed in the multi-trait but not the single-trait analysis. For single-trait analysis, an additional random effect can be included to partition the residual and GxE effects. The variance-covariance matrix of this effect must be named Omega (following notation from Damesa et al. 2017) and defined globally in the workspace, rather than passing it to the function (this is due to limitations with ASReml-R). The function \code{\link{Stage2_prep}} can be used to prepare both \code{data} and Omega. By default, the model includes independent random effects for genotype (id). Additional genetic effects with specific covariance structure (such as the G matrix for genomic breeding values) can be included using the argument \code{kernels}, which is a vector of variable names (for example, "G") defined in the global environment. (Do not use the name "I" for a kernel; it is reserved for the independent genetic effect.) All individuals in \code{data} must be present in the kernel matrices, but the kernels can contain individuals not in \code{data} to make predictions for unphenotyped individuals using \code{\link{predict_MME}}. All kernel matrices must have the same rownames attribute. For numerical stability when inverting the kernel matrices, a small positive number (1e-5) is added to the diagonal elements. By default, the workspace memory for ASReml-R is set at 500mb. If you get an error about insufficient memory, try increasing it. ASReml-R version 4.1.0.148 or later is required. For kernel matrix K, the variance reported in \code{vars} equal the variance component times the mean of the diagonal elements of ZKZ', to facilitate proper calculation of the proportion of variance.
#' 
#' @references Damesa et al. 2017. Agronomy Journal 109: 845-857. doi:10.2134/agronj2016.07.0395
#' 
#' @param data data frame of BLUEs from Stage 1 (see Details)
#' @param traits names of traits, matching columns in \code{data}
#' @param fixed optional names of fixed effects, matching columns in \code{data}
#' @param silent TRUE/FALSE, whether to suppress ASReml-R output
#' @param workspace Memory limit for ASRreml-R variance estimation
#' 
#' @return List containing
#' \describe{
#' \item{aic}{AIC}
#' \item{vars}{variances}
#' \item{fixed}{fixed effects}
#' \item{cov.mat}{genetic variance-covariance matrix for the traits or locations}
#' \item{uniplot}{uniplot of the genetic correlation between locations}
#' \item{MME}{variable of class \code{\link{MME}} for use with \code{\link{predict_MME}}}
#' }
#' 
#' @importFrom stats model.matrix
#' @importFrom methods new
#' @import Matrix
#' @importFrom tidyr pivot_longer
#' @export

Stage2 <- function(data,traits,fixed=NULL,silent=TRUE,workspace="500mb") {
  
  stopifnot(requireNamespace("asreml"))
  stopifnot(length(grep("trait",traits))==0)
  stopifnot(all(c("id","env") %in% colnames(data)))
  
  if ("loc" %in% colnames(data)) {
    data$loc <- factor(as.character(data$loc))
    locations <- levels(data$loc)
    n.loc <- length(locations)
    stopifnot(n.loc > 1)
  } else {
    n.loc <- 1
  }
  
  if (!is.null(fixed)) {
    stopifnot(fixed %in% colnames(data))
  }
  
  traits <- sort(traits)
  n.trait <- length(traits)
  if (n.trait > 1 & n.loc > 1) {
    stop("Genotype x Location effects are not supported for multi-trait analysis")
  }

  if (n.trait > 1) {
    model <- sub("blue",paste(traits,collapse=","),
                 "asreml(data=data,fixed=cbind(blue)~FIX,random=~RANDOM,residual=~id(units):us(trait)",fixed=T)
    if (!is.null(fixed)) {
      model <- sub("FIX",paste(paste0(c("env",fixed),":trait"),collapse="+"),model,fixed=T)
    } else {
      model <- sub("FIX","env:trait",model,fixed=T)
    }
  } else {
    stopifnot(!is.na(data[,traits]))
    model <- sub("blue",traits,"asreml(data=data,fixed=blue~FIX,random=~RANDOM,residual=~idv(units)",fixed=T)
    if (!is.null(fixed)) {
      model <- sub("FIX",paste(c("env",fixed),collapse="+"),model,fixed=T)
    } else {
      model <- sub("FIX","env",model,fixed=T)
    }
    #  model <- sub("blue",traits,"asreml(data=data,fixed=blue~env,random=~RANDOM,residual=~dsum(~id(units)|loc)",fixed=T)
  }
  
  data$env <- factor(as.character(data$env))
  data$id <- factor(as.character(data$id))
  id <- levels(data$id)
  
  if (exists("G")&&(inherits(G,"Matrix")|inherits(G,"matrix"))) {
    cat("G matrix detected\n")
    isG <- TRUE 
    colnames(G) <- rownames(G)
    G <- G + 1e-5*diag(nrow(G))
    stopifnot(is.element(id,rownames(G)))
    id <- sort(rownames(G))
    data$id <- factor(as.character(data$id),levels=id)
    G <- G[id,id]
    meanG <- mean(diag(G)[unique(as.character(data$id))])
  } else {
    isG <- FALSE
  }
  
  if (!isG) {
    if (n.loc > 1) {
      random.effects <- ifelse(n.loc==2,"id:corh(loc)","id:fa(loc,2)")
    } else {
      random.effects <- ifelse(n.trait>1,"id:us(trait)","id")
    }
  } else {
    if (n.trait > 1) {
      random.effects <- "id:idh(trait)+vm(id,source=G,singG='PSD'):us(trait)"
    } else {
      if (n.loc > 1) {
        if (n.loc==2) {
          random.effects <- "id:idh(loc)+vm(id,source=G,singG='PSD'):corh(loc)"
        } else {
          random.effects <- "id:idh(loc)+vm(id,source=G,singG='PSD'):fa(loc,2)"
        }
      } else {
        random.effects <- "id+vm(id,source=G,singG='PSD')"
      }
    }
  }
  
  if (!(exists("Omega")&&(inherits(Omega,"Matrix")|inherits(Omega,"matrix"))) | (n.trait > 1)) {
    skip.omega <- TRUE
  } else {
    cat("Omega matrix detected\n")
    skip.omega <- FALSE
    stopifnot(nrow(Omega)==nrow(data))
    random.effects <- paste0(random.effects,"+vm(units,source=Omega,singG='PSD')")
  }
  
  asreml::asreml.options(workspace=workspace,maxit=30,trace=!silent)
  model <- sub(pattern="RANDOM",replacement=random.effects,model,fixed=T)
  if (!skip.omega) {
    start.table <- eval(parse(text=paste0(model,",start.values = TRUE)")))$vparameters.table
    k <- grep("Omega",start.table$Component,fixed=T)
    start.table$Value[k] <- 1
    start.table$Constraint[k] <- "F"
    ans <- eval(parse(text=paste0(model,",G.param=start.table)")))
  } else {
    ans <- eval(parse(text=paste0(model,")")))
  }
  while (!ans$converge) {
    cat("ASReml-R failed to converge. Do you wish to continue running? y/n \n")
    input <- readLines(n=1)
    if (input=="y") {
      ans <- update.asreml(ans)
    } else {
      return()
    }
  }
  
  sans <- summary(ans,coef=TRUE)
  vc <- sans$varcomp
  if (n.trait==1) {
    remove <- c(grep("Omega",rownames(vc),fixed=T),match("units!R",rownames(vc)))
    vc <- vc[-remove,]
  } else {
    vc <- vc[vc$bound!="F",c("component","std.error")]
  }
  
  if (n.trait > 1) {
    cov.mat <- f.cor(vc[grep("units",rownames(vc),fixed=T),],traits)
    Rmat <- kronecker(Diagonal(nrow(data)),cov.mat)
    id.env <- paste(data$id,data$env,sep=":")
    tmp <- expand.grid(traits,id.env)
    Rnames <- apply(cbind(as.character(tmp$Var2),as.character(tmp$Var1)),1,paste,collapse=":")
    
    vars <- diag(cov.mat)
    data2 <- data[,c("id","env",traits)]
    data2 <- pivot_longer(data=data2,cols=match(traits,colnames(data2)),
                          names_to="trait",values_to="blue")
    data2$trait <- factor(data2$trait,levels=traits)
    ix <- which(!is.na(data2$blue))
    data2 <- as.data.frame(data2[ix,])
    tmp <- apply(data2[,c("id","env","trait")],1,paste,collapse=":")
    ix <- match(tmp,Rnames)
    Rmat <- Rmat[ix,ix]
  
    Z <- Matrix(model.matrix(~id:trait-1,data2))
    colnames(Z) <- sub("trait","",colnames(Z),fixed=T)
  } else {
    vars <- vc[match("units!units",rownames(vc)),1]
    Rmat <- Diagonal(nrow(data))*vars
    if (n.loc > 1) {
      data2 <- data[,c("id","env","loc",traits)]
      colnames(data2) <- c("id","env","loc","blue")
      Z <- Matrix(model.matrix(~id:loc-1,data2))
      colnames(Z) <- sub("loc","",colnames(Z),fixed=T)
    } else {
      data2 <- data[,c("id","env",traits)]
      colnames(data2) <- c("id","env","blue")
      Z <- Matrix(model.matrix(~id-1,data2))
    }
  }
  colnames(Z) <- sub("id","",colnames(Z),fixed=T)
  n.id <- length(id)
  
  if (!skip.omega) {
    vc.out <- rbind('g x env'=vars,residual=mean(diag(Omega)))
    Rmat <- Rmat + Matrix(Omega)
  } else {
    vc.out <- matrix(vars,nrow=1)
    rownames(vc.out) <- "residual"
  }
  colnames(vc.out) <- traits

  K <- vector("list",as.integer(isG)+1)

  if (n.trait==1) {
    if (n.loc==1) {
      vcnames <- rownames(vc)
      K[[2]] <- vc["id",1]*Diagonal(n.id)
      dimnames(K[[2]]) <- list(id,id)
      vars <- matrix(0,ncol=1,nrow=as.integer(isG)+1)
      if (isG) {
        j <- grep("source = G",vcnames,fixed=T)
        K[[1]] <- Matrix(vc[j,1]*G)
        vars[1,1] <- vc[j,1]*meanG
        vars[2,1] <- vc["id",1]
        rownames(vars) <- c("additive","non-additive")
      } else {
        vars[1,1] <- vc["id",1]
        rownames(vars) <- c("genotype")
      }
      vc.out <- rbind(vars,vc.out)
    } else {
      tmp <- expand.grid(locations,id)
      Knames <- apply(cbind(as.character(tmp$Var2),as.character(tmp$Var1)),1,paste,collapse=":")
      iz <- match(colnames(Z),Knames)
      if (!isG) {
        if (n.loc > 2) {
          iz <- grep("id:fa",rownames(vc),fixed=T)
        } else {
          iz <- grep("id:loc",rownames(vc),fixed=T)
        }
        cov.ans <- f.cov.loc(vc=vc[,],locations)
        cov.mat <- cov.ans$cov.mat
        Vg <- mean(cov.mat[upper.tri(cov.mat,diag=F)])
        VgL <- mean(diag(cov.mat)) - Vg
        KK <- kronecker(Diagonal(n.id),cov.mat)
        dimnames(KK) <- list(Knames,Knames)
        K[[1]] <- KK[iz,iz]
        vars <- matrix(c(Vg,VgL),ncol=1,nrow=2)
        rownames(vars) <- c("genotype","g x loc")
        vc.out <- rbind(vars,vc.out)
      } else {
        vars <- f.id(vc,locations,keyword="id:loc")
        KK <- kronecker(Diagonal(n.id),Diagonal(n=n.loc,x=vars))
        dimnames(KK) <- list(Knames,Knames)
        K[[2]] <- KK[iz,iz]
        Vr <- matrix(mean(vars),nrow=1)
        rownames(Vr) <- "non-additive"
        vc.out <- rbind(Vr,vc.out)
      
        cov.ans <- f.cov.loc(vc[grep("source = G",rownames(vc),fixed=T),],locations)
        cov.mat <- cov.ans$cov.mat
        Vg <- mean(cov.mat[upper.tri(cov.mat,diag=F)])
        VgL <- mean(diag(cov.mat)) - Vg
        Vg <- Vg*meanG
        VgL <- VgL*meanG
        KK <- kronecker(G,cov.mat)
        dimnames(KK) <- list(Knames,Knames)
        K[[1]] <- KK[iz,iz]
        vars <- matrix(c(Vg,VgL),ncol=1,nrow=2)
        rownames(vars) <- c("additive","add x loc")
        vc.out <- rbind(vars,vc.out)
      }
    }
  } else {
    #to do
    tmp <- expand.grid(traits,id)
    Knames <- apply(cbind(as.character(tmp$Var2),as.character(tmp$Var1)),1,paste,collapse=":")
    iz <- match(colnames(Z),Knames)
    
    #kernel = I
    if (nK==0) {
      cov.mat <- f.cor(vc[grep("id:trait",rownames(vc),fixed=T),],traits)
      KK <- kronecker(Diagonal(n.id),cov.mat)
      vc.out <- rbind(diag(cov.mat),vc.out)
    } else {
      vars <- f.id(vc,traits,keyword="id:trait")
      KK <- kronecker(Diagonal(n.id),Diagonal(n=n.trait,x=vars))
      vc.out <- rbind(vars,vc.out)
    }
    dimnames(KK) <- list(Knames,Knames)
    K[[1]] <- KK[iz,iz]
    rownames(vc.out) <- replace(rownames(vc.out),1,"kernel=I")
    
    if (nK > 0) {
      cov.mat <- f.cor(vc=vc[grep("):trait",rownames(vc),fixed=T),],traits)
      KK <- eval(parse(text=sub("Q",kernels[1],"Matrix(kronecker(Q,cov.mat))")))
      dimnames(KK) <- list(Knames,Knames)
      K[[2]] <- KK[iz,iz]
      tmp <- eval(parse(text=sub("Q","K[[2]]","diag(Z%*%Q%*%t(Z))")))
      vars <- tapply(tmp,data2$trait,mean)
      vc.out <- rbind(vars,vc.out)
      rownames(vc.out) <- replace(rownames(vc.out),1,paste0("kernel=",kernels[1]))
    }
    
    if (nK > 1) {
      for (i in 2:nK) {
        vars <- f.id(vc,traits,keyword=sub("Q",kernels[i],"source = Q"))
        KK <- eval(parse(text=sub("Q",kernels[i],"kronecker(Q,Diagonal(n=n.trait,x=vars))")))
        dimnames(KK) <- list(Knames,Knames)
        K[[i+1]] <- KK[iz,iz]
        tmp <- eval(parse(text=sub("Q","K[[i+1]]","diag(Z%*%Q%*%t(Z))")))
        vars <- tapply(tmp,data2$trait,mean)
        vc.out <- rbind(vars,vc.out)
        rownames(vc.out) <- replace(rownames(vc.out),1,paste0("kernel=",kernels[i]))
      }
    }
  }
  
  if (!is.null(fixed)) {
    beta <- coefficients(ans)$fixed
    beta <- beta[fixed,1]
    names(beta) <- fixed
    x <- as.matrix(data[,fixed]) %*% diag(beta)
    vars <- matrix(apply(x,2,var),ncol=1)
    rownames(vars) <- fixed
    vc.out <- rbind(vars,vc.out)
  }
  
  out <- new(Class="MME",data=data2,kernels=K,Rmat=Rmat)
  out <- list(aic=round(as.numeric(sans$aic),1),vars=vc.out,MME=out)
  if (!is.null(fixed)) {
    out <- c(out,list(fixed=beta))
  }
  if (n.trait > 1) {
    return(c(out,list(cov.mat=cov.mat,MME=out)))
  }
  if (n.loc > 1) {
    return(c(out,list(cov.mat=cov.ans$cov.mat,uniplot=cov.ans$plot)))
  }
  return(out)
}

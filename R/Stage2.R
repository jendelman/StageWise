#' Stage 2 analysis of multi-environment trials
#' 
#' Stage 2 analysis of multi-environment trials 
#' 
#' Stage 2 of the two-stage approach described by Damesa et al. 2017, using ASReml-R for variance component estimation. The variable \code{data} has three mandatory column: id, env, BLUE. Optionally, \code{data} can have a column labeled "loc", which changes the main effect for genotype into a separable genotype-within-location effect, using a FA2 covariance model for the locations. Optionally, \code{data} can have a column labeled "trait", which introduces an unstructured covariance model for the traits. The multi-location and multi-trait analyses cannot be combined. Missing data are allowed in the multi-trait but not the single-trait analysis. 
#' 
#' The argument \code{vcov} is used to partition the macro- and micro-environmental variation, which are called GxE and residual in the output. \code{vcov} is a named list of variance-covariance matrices for the BLUEs within each environment, with id on the rownames. The order in \code{vcov} and \code{data} should match. Both \code{data} and \code{vcov} can be created using the function \code{\link{Stage1}}. 
#' 
#' Because ASReml-R can only use relationship matrices defined in the global environment, this function creates and then removes global variables when either \code{vcov} or \code{geno} is used.
#' 
#' The argument \code{geno} is used to partition genetic values into additive and non-additive (g.resid) components. Any individuals in \code{data} that are not present in \code{geno} are discarded. For kernel matrix K, the variance reported in \code{vars} equal the variance component times the mean of the diagonal elements of ZKZ', to facilitate proper calculation of the proportion of variance.
#' 
#' By default, the workspace memory for ASReml-R is set at 500mb. If you get an error about insufficient memory, try increasing it. ASReml-R version 4.1.0.148 or later is required. 
#' 
#' @references Damesa et al. 2017. Agronomy Journal 109: 845-857. doi:10.2134/agronj2016.07.0395
#' 
#' @param data data frame of BLUEs from Stage 1 (see Details)
#' @param vcov list of variance-covariance matrices for the BLUEs
#' @param geno output from \code{\link{read_geno}}
#' @param fix.eff.marker markers in \code{geno} to include as additive fixed effect covariates
#' @param silent TRUE/FALSE, whether to suppress ASReml-R output
#' @param workspace Memory limit for ASRreml-R variance estimation
#' 
#' @return List containing
#' \describe{
#' \item{aic}{AIC}
#' \item{vars}{variances, as variable of class \code{\link{class_var}}}
#' \item{fixed}{Fixed effect estimates for env and markers}
#' \item{random}{Random effect predictions}
#' \item{uniplot}{uniplot of the genetic correlation between locations}
#' }
#' 
#' @importFrom stats model.matrix var
#' @importFrom methods new
#' @importFrom rlang .data
#' @import Matrix
#' @import ggplot2
#' @import ggrepel
#' @export

Stage2 <- function(data,vcov=NULL,geno=NULL,fix.eff.marker=NULL,silent=TRUE,workspace="500mb") {
  
  stopifnot(inherits(data,"data.frame"))
  stopifnot(c("id","env","BLUE") %in% colnames(data))
  data$id <- as.character(data$id)
  data$env <- factor(as.character(data$env))
  
  if (!is.null(geno)) {
    stopifnot(inherits(geno,"class_geno"))
    id <- sort(intersect(data$id,rownames(geno@G)))
    .GlobalEnv$asremlG <- geno@G[id,id] + Diagonal(length(id))*1e-6
    meanG <- mean(diag(.GlobalEnv$asremlG))
    data <- data[data$id %in% id,]
  } else {
    meanG <- as.numeric(NA)
  }
  
  if (!is.null(fix.eff.marker)) {
    stopifnot(!is.null(geno))
    n.mark <- length(fix.eff.marker)
    stopifnot(fix.eff.marker %in% colnames(geno@coeff))
    dat2 <- data.frame(id=rownames(geno@coeff),as.matrix(geno@coeff[,fix.eff.marker]))
    colnames(dat2) <- c("id",fix.eff.marker)
    data <- merge(data,dat2,by="id")
  }
  
  id <- sort(unique(data$id))
  data$id <- factor(data$id,levels=id)
  stopifnot(table(data$env)>0)
  
  if ("loc" %in% colnames(data)) {
    data$loc <- factor(as.character(data$loc))
    data <- data[order(data$loc),]
    locations <- levels(data$loc)
    n.loc <- length(locations)
    stopifnot(n.loc > 1)
  } else {
    n.loc <- 1
  }
  
  if ("trait" %in% colnames(data)) {
    stopifnot(n.loc==1)
    data$Trait <- factor(data$trait)
    traits <- levels(data$Trait)
    n.trait <- length(traits)
    stopifnot(n.trait > 1)
  } else {
    n.trait <- 1
  }

  if (n.trait==1) {
    if (n.loc==1) {
      model <- "asreml(data=data,fixed=BLUE~FIX-1,random=~RANDOM,residual=~idv(units)"
    } else {
      model <- "asreml(data=data,fixed=BLUE~FIX-1,random=~RANDOM,residual=~dsum(~idv(units)|loc)"
    }
    if (!is.null(fix.eff.marker)) {
      if (n.loc==1) {
        model <- sub("FIX",paste(c("env",fix.eff.marker),collapse="+"),model,fixed=T)
      } else {
        model <- sub("FIX",paste(c("env",paste0(fix.eff.marker,":loc")),collapse="+"),model,fixed=T)
      }
    } else {
      model <- sub("FIX","env",model,fixed=T)
    }
  } else {
    model <- "asreml(data=data,fixed=BLUE~FIX-1,random=~RANDOM,residual=~id(units):us(trait)"
    if (!is.null(fix.eff.marker)) {
      model <- sub("FIX",paste(paste0(c("env",fix.eff.marker),":Trait"),collapse="+"),model,fixed=T)
    } else {
      model <- sub("FIX","env:Trait",model,fixed=T)
    }
  }
  
  if (is.null(geno)) {
    if (n.loc > 1) {
      random.effects <- ifelse(n.loc==2,"id:corh(loc)","id:fa(loc,2)")
    } else {
      random.effects <- ifelse(n.trait>1,"id:us(Trait)","id")
    }
  } else {
    if (n.trait > 1) {
      random.effects <- "id:us(Trait)+vm(id,source=asremlG,singG='PSD'):us(Trait)"
    } else {
      if (n.loc > 1) {
        if (n.loc==2) {
          random.effects <- "id:corh(loc)+vm(id,source=asremlG,singG='PSD'):corh(loc)"
        } else {
          random.effects <- "id:corh(loc)+vm(id,source=asremlG,singG='PSD'):fa(loc,2)"
        }
      } else {
        random.effects <- "id+vm(id,source=asremlG,singG='PSD')"
      }
    }
  }
  
  if (!is.null(vcov)) {
    stopifnot(is.list(vcov))
    stopifnot(sort(levels(data$env))==sort(names(vcov)))
    id.ix <- lapply(vcov,function(y){which(rownames(y) %in% id)})
    omega.list<- mapply(FUN=function(Q,ix){as(Q[ix,ix],"dpoMatrix")},Q=vcov,ix=id.ix)
    diagOmega <- unlist(mapply(FUN=function(Q,ix){diag(Q)[ix]},Q=vcov,ix=id.ix))
    meanOmega <- mean(tapply(diagOmega,data$id,mean))
    .GlobalEnv$asremlOmega <- direct_sum(lapply(omega.list,solve))
    attr(.GlobalEnv$asremlOmega,"INVERSE") <- TRUE
    random.effects <- paste0(random.effects,"+vm(units,source=asremlOmega)")
  } else {
    meanOmega <- as.numeric(NA)
  }
  
  asreml.options(workspace=workspace,maxit=30,trace=!silent)
  model <- sub(pattern="RANDOM",replacement=random.effects,model,fixed=T)
  if (!is.null(vcov)) {
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
  if (n.loc > 1) {
    out <- vector("list",5)
    names(out) <- c("aic","vars","fixed","random","uniplot")
  } else {
    out <- vector("list",4)
    names(out) <- c("aic","vars","fixed","random")
  }
  out$aic <- as.numeric(sans$aic)
  
  #fixed effects
  if (n.trait==1) {
    beta <- sans$coef.fixed
    rownames(beta) <- gsub("env_","",rownames(beta))
    ix <- match(levels(data$env),rownames(beta))
    out$fixed$env <- data.frame(env=levels(data$env),effect=as.numeric(beta[ix,1]))
  
    beta.stat <- matrix(NA,nrow=0,ncol=0)
    marker.cov <- matrix(NA,nrow=0,ncol=0)
    if (!is.null(fix.eff.marker)) {
      var.x <- var(as.matrix(geno@coeff[,fix.eff.marker]))
      if (n.loc==1) {
        ix <- match(fix.eff.marker,rownames(beta))
        out$fixed$marker <- data.frame(marker=fix.eff.marker,effect=as.numeric(beta[ix,1]))
        beta.stat <- matrix(diag(var.x) * beta[ix,1]^2,ncol=1)
        rownames(beta.stat) <- fix.eff.marker
      } else {
        rownames(beta) <- gsub("loc_","",rownames(beta))
        ix <- unlist(lapply(fix.eff.marker,grep,rownames(beta),fixed=T))
        loc.marker <- expand.grid(loc=locations,marker=fix.eff.marker,stringsAsFactors = FALSE)
        out$fixed$marker <- data.frame(marker=loc.marker$marker,
                                       loc=loc.marker$loc,
                                       effect=as.numeric(beta[ix,1]))

        tmp <- expand.grid(1:n.loc,1:n.loc)
        tmp <- tmp[tmp$Var2 >= tmp$Var1,]
        tmp$value <- numeric(nrow(tmp))
        x <- out$fixed$marker
        for (i in 1:nrow(tmp)) {
          lv <- matrix(x$effect[x$loc==locations[tmp$Var1[i]]],nrow=1)
          rv <- matrix(x$effect[x$loc==locations[tmp$Var2[i]]],ncol=1)
          tmp$value[i] <- lv %*% var.x %*% rv
        }
        marker.cov <- as.matrix(sparseMatrix(i=tmp$Var1,j=tmp$Var2,x=tmp$value,dims=c(n.loc,n.loc),
                                             dimnames=list(locations,locations),symmetric=T))
        tmp <- split(out$fixed$marker$effect,out$fixed$marker$marker)
        tmp2 <- lapply(tmp,function(x){crossprod(matrix(x,nrow=1))})
        beta.stat <- t(mapply(function(x,y){partition(x*y)},x=tmp2,y=as.list(diag(var.x))))
        colnames(beta.stat) <- c("marker","marker x loc")
      }
    }
  } else {
    #multi-trait
  }
    
  #variances 
  vc <- sans$varcomp
  vc <- vc[-which(vc$bound=="F" & vc$component==1),c("component","std.error")]
  vc.names <- rownames(vc)
  
  if (n.trait==1) {
    resid.vc <- Matrix(vc[grep("units",vc.names,fixed=T),1],ncol=1)
    if (n.loc > 1) {
      rownames(resid.vc) <- locations
      if (is.null(geno)) {
        if (n.loc > 2) {
          iz <- grep("id:fa",vc.names,fixed=T)
        } else {
          iz <- grep("id:loc",vc.names,fixed=T)
        }
        cov.ans <- f.cov.loc(vc=vc[iz,],locations)
        g.resid.vc <- cov.ans$cov.mat
        add.vc <- Matrix(NA,nrow=0,ncol=0)
      } else {
        g.resid.vc <- Matrix(vc[grep("id:loc",vc.names,fixed=T),1],ncol=1)
        rownames(g.resid.vc) <- c("cor",locations)
        cov.ans <- f.cov.loc(vc[grep("source = asremlG",vc.names,fixed=T),],locations)
        add.vc <- cov.ans$cov.mat
      }
      out$uniplot <- cov.ans$plot
    } else {
      g.resid.vc <- Matrix(vc[match("id",vc.names),1],ncol=1)
      if (!is.null(geno)) {
        add.vc <- Matrix(vc[grep("source = asremlG",vc.names,fixed=T),1],ncol=1)
      } else {
        add.vc <- Matrix(NA,nrow=0,ncol=0)
      }
    }
  } else {
    #multi-trait

  }
  
  out$vars <- new(Class="class_var",add=add.vc,g.resid=g.resid.vc,resid=resid.vc,
                  meanG=meanG,meanOmega=meanOmega,
                  fixed.marker.var=beta.stat,fixed.marker.cov=marker.cov)
  
  #random effects
  if (n.trait==1) {
    u <- sans$coef.random 
    if (n.loc > 1) {
      id.loc.names <- expand.grid(loc=locations,id=id,stringsAsFactors = F)[,c(2,1)]
      if (!is.null(geno)) {
        id.loc <- expand.grid(paste0("loc_",locations),paste0("id_",id),stringsAsFactors = F)[,c(2,1)]
        ix2 <- match(apply(id.loc,1,paste,collapse=":"),rownames(u))
        if (n.loc==2) {
          id.loc <- expand.grid(paste0("loc_",locations),paste0("vm(id, source = asremlG, singG = \"PSD\")_",id),stringsAsFactors = F)[,c(2,1)]
        } else {
          id.loc <- expand.grid(paste0("fa(loc, 2)_",locations),paste0("vm(id, source = asremlG, singG = \"PSD\")_",id),stringsAsFactors = F)[,c(2,1)]
        }
        ix1 <- match(apply(id.loc,1,paste,collapse=":"),rownames(u))
        out$random <- data.frame(id.loc.names,add=as.numeric(u[ix1,1]),g.resid=as.numeric(u[ix2,1]))
      } else {
        if (n.loc==2) {
          id.loc <- expand.grid(paste0("loc_",locations),paste0("id_",id),stringsAsFactors = F)[,c(2,1)]
        } else {
          id.loc <- expand.grid(paste0("fa(loc, 2)_",locations),paste0("id_",id),stringsAsFactors = F)[,c(2,1)]
        }
        ix <- match(apply(id.loc,1,paste,collapse=":"),rownames(u))
        out$random <- data.frame(id.loc.names,value=as.numeric(u[ix,1]))
      }
    } else {
      if (!is.null(geno)) {
        ix1 <- match(paste("vm(id, source = asremlG, singG = \"PSD\")",id,sep="_"),rownames(u))
        ix2 <- match(paste("id",id,sep="_"),rownames(u))
        out$random <- data.frame(id=id,add=as.numeric(u[ix1,1]),g.resid=as.numeric(u[ix2,1]))
      } else {
        ix <- match(paste("id",id,sep="_"),rownames(u))
        out$random <- data.frame(id=id,value=as.numeric(u[ix,1]))
      }
    }
  } else {
    #multi-trait
  }
  
  if (!is.null(geno))
    rm("asremlG",envir = .GlobalEnv)
  if (!is.null(vcov))
    rm("asremlOmega",envir = .GlobalEnv)
  
  return(out)
}

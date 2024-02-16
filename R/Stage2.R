#' Stage 2 analysis of multi-environment trials
#' 
#' Stage 2 analysis of multi-environment trials 
#' 
#' Stage 2 of the two-stage approach described by Damesa et al. 2017, using ASReml-R for variance component estimation. The variable \code{data} has three mandatory column: id, env, BLUE. Optionally, \code{data} can have a column labeled "loc", which changes the main effect for genotype into a separable genotype-within-location effect, using a FA2 covariance model for the locations. Optionally, \code{data} can have a column labeled "trait", which uses an unstructured covariance model. The multi-location and multi-trait analyses cannot be combined. Missing data are allowed in the multi-trait but not the single-trait analysis. The argument \code{geno} is used to partition genetic values into additive and non-additive components. Any individuals in \code{data} that are not present in \code{geno} are discarded. 
#' 
#' The argument \code{vcov} is used to partition the macro- and micro-environmental variation, which are called GxE and residual in the output. \code{vcov} is a named list of variance-covariance matrices for the BLUEs within each environment, with id for rownames (single trait) or id:trait. The order in \code{vcov} and \code{data} should match. Both \code{data} and \code{vcov} can be created using the function \code{\link{Stage1}}. 
#' 
#' Because ASReml-R can only use relationship matrices defined in the global environment, this function creates and then removes global variables when either \code{vcov} or \code{geno} is used. By default, the workspace memory for ASReml-R is set at 500mb. If you get an error about insufficient memory, try increasing it. ASReml-R version 4.1.0.148 or later is required. 
#' 
#' The \code{covariates} option is only available for single trait/loc analysis.
#' 
#' Argument \code{pairwise} was added in package version 1.04, which specifies that multi-trait analysis is performed as multiple bivariate analyses, which often converges better. The returned object is a list of the results from the bivariate analyses, as well as "vars" for all traits, which is needed for \code{\link{blup_prep}}.
#' 
#' @references Damesa et al. 2017. Agronomy Journal 109: 845-857. doi:10.2134/agronj2016.07.0395
#' 
#' @param data data frame of BLUEs from Stage 1 (see Details)
#' @param vcov named list of variance-covariance matrices for the BLUEs
#' @param geno output from \code{\link{read_geno}}
#' @param fix.eff.marker markers in \code{geno} to include as additive fixed effect covariates
#' @param silent TRUE/FALSE, whether to suppress ASReml-R output
#' @param workspace Memory limit for ASRreml-R variance estimation
#' @param non.add one of the following: "none","g.resid","dom"
#' @param max.iter maximum number of iterations for asreml
#' @param covariates names of other covariates in \code{data}
#' @param pairwise TRUE/FALSE should multi-trait analysis proceed pairwise
#' 
#' @return List containing
#' \describe{
#' \item{aic}{AIC}
#' \item{vars}{variance components for \code{\link{blup_prep}}, as variable of class \code{\link{class_var}}}
#' \item{params}{Estimates and SE for fixed effects and variance components}
#' \item{random}{Random effect predictions}
#' \item{loadings}{scaled loadings for the FA2 multi-loc model}
#' }
#' 
#' @importFrom stats model.matrix var
#' @importFrom methods new
#' @importFrom rlang .data
#' @import CVXR
#' @import Matrix
#' @import ggplot2
#' @import ggrepel
#' 
#' @export

Stage2 <- function(data,vcov=NULL,geno=NULL,fix.eff.marker=NULL,
                   silent=TRUE,workspace="500mb",non.add="g.resid",max.iter=20,
                   covariates=NULL,pairwise=FALSE) {
  
  stopifnot(is(data,"data.frame"))
  stopifnot(requireNamespace("asreml"))
  stopifnot(non.add %in% c("none","g.resid","dom"))
  if (non.add=="dom")
    stopifnot(is(geno,"class_genoD"))
  library(asreml)
  stopifnot(c("id","env","BLUE") %in% colnames(data))
  
  if ("trait" %in% colnames(data) & pairwise) {
    data$trait <- as.character(data$trait)
    traits <- unique(data$trait)
    n.trait <- length(traits)
    
    if (n.trait > 2) {
      geno1 <- matrix(0,n.trait,n.trait)
      dimnames(geno1) <- list(traits,traits)
      resid.vc <- B <- geno1
      
      if (non.add=="none") {
        geno2 <- matrix(0,0,0)
      } else {
        geno2 <- geno1
      }
      pairs <- t(combn(traits,2))
      npair <- nrow(pairs)
      result <- vector("list",npair)
      for (i in 1:npair) {
        if (!silent)
          print(paste(c("Pairwise:",pairs[i,]),collapse=" "))
        
        data2 <- data[data$trait %in% pairs[i,],]
        if (!is.null(vcov)) {
          n.env <- length(vcov)
          vcov2 <- vector("list",n.env)
          names(vcov2) <- names(vcov)
          for (j in 1:n.env) {
            dname <- strsplit(rownames(vcov[[j]]),split=":",fixed=T)
            traits <- sapply(dname,"[[",2)
            ix <- which(traits %in% pairs[i,])
            vcov2[[j]] <- vcov[[j]][ix,ix]
          }
        } else {
          vcov2 <- NULL
        }
        result[[i]] <- Stage2(data2,vcov2,geno,fix.eff.marker,
            silent,workspace,non.add,max.iter,covariates,pairwise=FALSE)
        trait.pair <- rownames(result[[i]]$vars@geno1)
        geno1[trait.pair,trait.pair] <- geno1[trait.pair,trait.pair] + as.matrix(result[[i]]$vars@geno1)
        B[trait.pair,trait.pair] <- B[trait.pair,trait.pair] + as.matrix(result[[i]]$vars@B)
        resid.vc[trait.pair,trait.pair] <- resid.vc[trait.pair,trait.pair] + as.matrix(result[[i]]$vars@resid)
        if (non.add!="none")
          geno2[trait.pair,trait.pair] <- geno2[trait.pair,trait.pair] + as.matrix(result[[i]]$vars@geno2)
      }
      diag(geno1) <- diag(geno1)/(n.trait-1)
      geno1 <- fu(geno1)
      diag(resid.vc) <- diag(resid.vc)/(n.trait-1)
      resid.vc <- fu(resid.vc)
      diag(B) <- diag(B)/(n.trait-1)
      B <- fu(B)
      
      if (non.add!="none") {
        diag(geno2) <- diag(geno2)/(n.trait-1)
        geno2 <- fu(geno2)
      }
      
      return(list(pairs=result, 
                  vars=new(Class="class_var",
                           geno1=geno1, geno2=geno2, resid=resid.vc, B=as.matrix(B),
                           model=result[[1]]$vars@model,
                           diagG=mean(sapply(result,function(x){x$vars@diagG})),
                           diagD=mean(sapply(result,function(x){x$vars@diagD})),
                           vars=array(NA,dim=c(0,0,0)),
                           fix.eff.marker=result[[1]]$vars@fix.eff.marker)))
    }  
  }
  
  data$id <- as.character(data$id)
  data$env <- as.character(data$env)
  data$env.id <- apply(data[,c("env","id")],1,paste,collapse=":")
  
  missing <- which(is.na(data$BLUE))
  if (length(missing) > 0) {
    data <- data[!(data$env.id %in% data$env.id[missing]),]
  }
  
  diagG <- diagD <- numeric(0)
  dom <- NULL
  if (!is.null(geno)) {
    stopifnot(is(geno,"class_geno"))
    id <- sort(intersect(data$id,rownames(geno@G)))
    n <- length(id)
    data <- data[data$id %in% id,]
    id.weights <- table(factor(data$id,levels=id))
    .GlobalEnv$asremlG <- geno@G[id,id]
    meanG <- as.numeric(pvar(V=.GlobalEnv$asremlG,weights=id.weights))
    diagG <- mean(diag(.GlobalEnv$asremlG))

    if (non.add=="dom") {
      .GlobalEnv$asremlD <- geno@D[id,id]
      meanD <- as.numeric(pvar(V=.GlobalEnv$asremlD,weights=id.weights))

      dom.covariate <- as.numeric(geno@coeff.D[id,] %*% matrix(1,nrow=ncol(geno@coeff.D),ncol=1))
      dom.covariate <- dom.covariate/(geno@scale*(geno@ploidy-1))
      names(dom.covariate) <- id
      data$dom <- dom.covariate[data$id]
  
      diagD <- mean(diag(.GlobalEnv$asremlD))
      dom <- "dom"
    } 
  } else {
    id <- sort(unique(data$id))
    n <- length(id)
  }
  
  if (!is.null(fix.eff.marker)) {
    stopifnot(!is.null(geno))
    n.mark <- length(fix.eff.marker)
    stopifnot(fix.eff.marker %in% colnames(geno@coeff))
    dat2 <- data.frame(id=rownames(geno@coeff),as.matrix(geno@coeff[,fix.eff.marker]))
    colnames(dat2) <- c("id",fix.eff.marker)
    data <- merge(data,dat2,by="id")
  } else {
    n.mark <- 0
  }
  
  envs <- unique(data$env)
  n.env <- length(envs)
  if (n.env==1 & (non.add=="g.resid")) {
    stop("Need more than one environment for g.resid model")
  }
  data$env <- factor(data$env,levels=envs)
  data$id <- factor(data$id,levels=id)
  env.id <- sort(unique(data$env.id))
  data$env.id <- factor(data$env.id,levels=env.id)
  
  out <- vector("list",4)
  names(out) <- c("aic","vars","fixed","random")
  n.loc <- 1
  n.trait <- 1
  
  if ("trait" %in% colnames(data)) {
    data$trait <- as.character(data$trait)
    traits <- sort(unique(data$trait))
    n.trait <- length(traits)
    stopifnot(n.trait > 1)
    
    data$env.id.trait <- apply(data[,c("env.id","trait")],1,paste,collapse=":")
    env.id.trait <- unique(data$env.id.trait)
    data$env.id.trait <- factor(data$env.id.trait,levels=env.id.trait)
    
    if (!is.null(vcov)) {
      dname <- strsplit(rownames(vcov[[1]])[1:n.trait],split=":",fixed=T)
      traits <- sapply(dname,"[[",2)
    } 
    data$Trait <- factor(data$trait,levels=traits)
    data <- data[order(data$env.id,data$Trait),]
    
  } else {
    
    traits <- ""  #NULL
    if ("loc" %in% colnames(data)) {
      data$loc <- factor(as.character(data$loc))
      data <- data[order(data$loc),]
      locations <- levels(data$loc)
      n.loc <- length(locations)
      stopifnot(n.loc > 1)
      loc.weights <- table(data$loc)
      data$gL <- paste(as.character(data$id),as.character(data$loc),sep=":")
      tmp <- expand.grid(loc=locations,id=id)
      gL <- apply(tmp[,c("id","loc")],1,paste,collapse=":")
      data$gL <- factor(data$gL,levels=gL)
      gL.weights <- table(data$gL)
    } 
  }
  
  if (!is.null(covariates)) {
    stopifnot(covariates %in% colnames(data))
    stopifnot(sapply(data[,covariates],class) %in% c("numeric","integer"))
    if (n.trait > 1 | n.loc > 1)
      stop("Additional covariates only supported for single trait/location.")
  }
  
  if (is.null(geno)) {
    tmp <- c("env","covariates","fixed.marker","additive","add x loc","dominance","heterosis","genotype","g x loc","g x env","Stage1.error","residual")
  } else {
    tmp <- c("env","covariates","fixed.marker","additive","add x loc","dominance","heterosis","g.resid","g x loc","g x env","Stage1.error","residual")
  }
  vars <- array(data=numeric(0),dim=c(n.trait,n.trait,12),dimnames=list(traits,traits,tmp))

  if (n.trait==1) {
    if (n.loc==1) {
      model <- "asreml(data=data,fixed=BLUE~FIX-1,random=~RANDOM,residual=~idv(units)"
      model <- sub("FIX",paste(c("env",covariates,fix.eff.marker,dom),collapse="+"),model,fixed=T)
    } else {
      model <- "asreml(data=data,fixed=BLUE~FIX-1,random=~RANDOM,residual=~dsum(~idv(units)|loc)"
      tmp <- gsub("+",":loc+",paste(c("env",fix.eff.marker,dom),collapse="+"),fixed=T)
      model <- sub("FIX",paste(tmp,"loc",sep=":"),model,fixed=T)
    }
    #if (n.mark > 0 | (non.add=="dom") | !is.null(fix.eff)) {
    #if (n.loc==1) {
    #  model <- sub("FIX",paste(c("env",fix.eff,fix.eff.marker,dom),collapse="+"),model,fixed=T)
    #} else {
    #  model <- sub("FIX",paste(c("env",paste0(c(fix.eff.marker,dom),":loc")),collapse="+"),model,fixed=T)
    #}
    #} else {
    #  model <- sub("FIX","env",model,fixed=T)
    #}
  } else {
    model <- "asreml(data=data,fixed=BLUE~FIX-1,random=~RANDOM,residual=~id(env.id):us(Trait)"
    tmp <- gsub("+",":Trait+",paste(c("env",fix.eff.marker,dom),collapse="+"),fixed=T)
    model <- sub("FIX",paste(tmp,"Trait",sep=":"),model,fixed=T)
    
    # if (n.mark > 0 | (non.add=="dom")) {
    #   model <- sub("FIX",paste(paste0(c("env",fix.eff.marker,dom),":Trait"),collapse="+"),model,fixed=T)
    # } else {
    #   model <- sub("FIX","env:Trait",model,fixed=T)
    # }
  }
  
  if (is.null(geno)) {
    if (n.loc > 1) {
      random.effects <- ifelse(n.loc==2,"id:corh(loc)","id:fa(loc,2)")
    } else {
      if (n.trait > 1) {
        random.effects <- ifelse(n.trait==2,"id:corh(Trait)","id:us(Trait)")
      } else {
        random.effects <- "id"
      }
    }
  } else {
    if (n.trait > 1) {
      if (n.trait==2) {
        random.effects <- "vm(id,source=asremlG,singG='PSD'):corh(Trait)"
      } else {
        random.effects <- "vm(id,source=asremlG,singG='PSD'):us(Trait)"
      }
      if (non.add=="g.resid") {
        if (n.trait==2) {
          random.effects <- paste(random.effects,"id:corh(Trait)",sep="+")
        } else {
          random.effects <- paste(random.effects,"id:us(Trait)",sep="+")
        }
      } 
      if (non.add=="dom") {
        if (n.trait==2) {
          random.effects <- paste(random.effects,"vm(id,source=asremlD,singG='PSD'):corh(Trait)",sep="+")
        } else {
          random.effects <- paste(random.effects,"vm(id,source=asremlD,singG='PSD'):us(Trait)",sep="+")
        }
      }
    } else {
      if (n.loc > 1) {
        if (n.loc==2) {
          random.effects <- "vm(id,source=asremlG,singG='PSD'):corh(loc)"
        } else {
          random.effects <- "vm(id,source=asremlG,singG='PSD'):fa(loc,2)"
        }
      } else {
        random.effects <- "vm(id,source=asremlG,singG='PSD')"
      }
      
      if (non.add=="g.resid") {
        if (n.loc > 1) {
          random.effects <- paste(random.effects,"id:corh(loc)",sep="+")
        } else {
          random.effects <- paste(random.effects,"id",sep="+")
        }
      }
      if (non.add=="dom") {
        if (n.loc > 1) {
          random.effects <- paste(random.effects,"vm(id,source=asremlD,singG='PSD'):corh(loc)",sep="+")
        } else {
          random.effects <- paste(random.effects,"vm(id,source=asremlD,singG='PSD')",sep="+")
        }
      }
    }
  }
  
  n.gE <- nrow(data)/n.trait
  if (!is.null(vcov)) {
    stopifnot(is.list(vcov))
    stopifnot(envs %in% names(vcov))
    vcov <- vcov[envs]
    vcov <- mapply(FUN=function(x,y){
      tmp <- paste(y,rownames(x),sep=":")
      dimnames(x) <- list(tmp,tmp)
      return(x)},x=vcov,y=as.list(names(vcov)))
    if (n.trait==1) {
      ix <- lapply(vcov,function(y){which(rownames(y) %in% env.id)})
      random.effects <- paste0(random.effects,"+vm(env.id,source=asremlOmega)")
    } else {
      ix <- lapply(vcov,function(y){which(rownames(y) %in% env.id.trait)})
      random.effects <- paste0(random.effects,"+vm(env.id.trait,source=asremlOmega)")
    }
    
    omega.list <- vector("list",n.env)
    for (j in 1:n.env) {
      omega.list[[j]] <- as(vcov[[j]][ix[[j]],ix[[j]]],"dpoMatrix")
    }
    .GlobalEnv$asremlOmega <- direct_sum(lapply(omega.list,solve))
    dname <- lapply(omega.list,rownames)
    dimnames(.GlobalEnv$asremlOmega) <- list(unlist(dname),unlist(dname))
    attr(.GlobalEnv$asremlOmega,"INVERSE") <- TRUE
    
    Omega <- bdiag(omega.list)
    ix <- split(1:(n.gE*n.trait),rep(1:n.trait,times=n.gE))
    for (i in 1:n.trait) 
      for (j in i:n.trait) {
        Omega_ij <- Omega[ix[[i]],ix[[j]]]
        vars[i,j,"Stage1.error"] <- mean(diag(Omega_ij)) - mean(Omega_ij)
      }
  }

  asreml::asreml.options(workspace=workspace,maxit=max.iter,trace=!silent)
  model <- sub(pattern="RANDOM",replacement=random.effects,model,fixed=T)

  start.table <- eval(parse(text=paste0(model,",start.values = TRUE)")))$vparameters.table
  if (!is.null(vcov)) {
    k <- grep("Omega",start.table$Component,fixed=T)
    start.table$Value[k] <- 1
    start.table$Constraint[k] <- "F"
  }
  if (!is.null(geno) & (n.loc > 1)) {
    k <- grep(":loc!loc!cor",start.table$Component,fixed=T)
    k2 <- grep("asremlG",start.table$Component,fixed=T)
    k <- setdiff(k,k2)
    start.table$Value[k] <- 0.9999
    start.table$Constraint[k] <- "F"
  }
  ans <- eval(parse(text=paste0(model,",G.param=start.table)")))
  if (!ans$converge) {
    stop("ASReml-R did not converge. Try increasing max.iter")
  }
  # while (!ans$converge) {
  #   cat("ASReml-R failed to converge. Do you wish to continue running? y/n \n")
  #   input <- readLines(n=1)
  #   if (input=="y") {
  #     ans <- asreml::update.asreml(ans)
  #   } else {
  #     return()
  #   }
  # }
  
  sans <- summary(ans,coef=TRUE)
  out$aic <- as.numeric(sans$aic)
  
  B <- matrix(0,nrow=n.trait,ncol=n.trait)
  dimnames(B) <- list(traits,traits)
  
  #fixed effects
  if (n.trait==1) {
    beta <- sans$coef.fixed
    beta.names <- rownames(beta)
    beta.names <- gsub("env_","",beta.names)
    ix <- match(levels(data$env),beta.names)
    if (any(is.na(ix))) {
      beta.names <- sapply(strsplit(beta.names,split=":",fixed=T),"[[",1)
      ix <- match(levels(data$env),beta.names)
    }
    out$params$env <- data.frame(env=levels(data$env),
                                 estimate=as.numeric(beta[ix,1]),
                                 SE=as.numeric(beta[ix,2]))
    vars[1,1,"env"] <- pvar(mu=out$params$env$estimate,weights=as.numeric(table(data$env)))
  
    if (non.add=="dom") {
      ix <- grep("dom",beta.names,fixed=T)
      if (n.loc==1) {
        vars[1,1,"heterosis"] <- pvar(mu=dom.covariate*beta[ix,1],weights=id.weights)
        out$params$heterosis <- data.frame(estimate=as.numeric(beta[ix,1]),
                                           SE=as.numeric(beta[ix,2]))
      } else {
        beta.names[ix] <- gsub("loc_","",beta.names[ix])
        beta.names[ix] <- gsub("dom","",beta.names[ix])
        beta.names[ix] <- gsub(":","",beta.names[ix])
        ix <- ix[match(locations,beta.names[ix])]
        mu <- as.numeric(kronecker(matrix(dom.covariate,ncol=1),beta[ix,1]))
        vars[1,1,"heterosis"] <- pvar(mu=mu,weights=gL.weights)
        out$params$heterosis <- data.frame(loc=locations,
                                           estimate=as.numeric(beta[ix,1]),
                                           SE=as.numeric(beta[ix,2]))
      }
    }
    if (!is.null(covariates)) {
      ix <- match(covariates,beta.names)
      out$params$covariates <- data.frame(name=covariates,
                                      estimate=as.numeric(beta[ix,1]),
                                      SE=as.numeric(beta[ix,2]))
      mu <- as.numeric(as.matrix(data[,covariates,drop=FALSE])%*%beta[ix,1])
      vars[1,1,"covariates"] <- pvar(mu=mu)
    }
    
    if (n.mark > 0) {
      if (n.loc==1) {
        ix <- match(fix.eff.marker,beta.names)
        out$params$marker <- data.frame(marker=fix.eff.marker,
                                        estimate=as.numeric(beta[ix,1]),
                                        SE=as.numeric(beta[ix,2]))
        mu <- as.numeric(as.matrix(geno@coeff[id,fix.eff.marker])%*%beta[ix,1])
        vars[1,1,"fixed.marker"] <- pvar(mu=mu,weights=id.weights)
      } else {
        ix <- grep("loc_",beta.names)
        tmp <- strsplit(beta.names[ix],split=":",fixed=T)
        beta.names[ix] <- sapply(tmp,function(z){
          tmp2 <- apply(array(z),1,grep,pattern="loc_")
          paste(z[order(sapply(tmp2,length))],collapse=":")
        })
        beta.names[ix] <- gsub("loc_","",beta.names[ix])
        loc.marker <- expand.grid(loc=locations,marker=fix.eff.marker,
                                  stringsAsFactors = FALSE)
        loc.marker <- loc.marker[,c(2,1)]
        ix <- match(apply(loc.marker,1,paste,collapse=":"),beta.names)
        x <- out$params$marker <- data.frame(marker=loc.marker$marker,
                                       loc=loc.marker$loc,
                                       estimate=as.numeric(beta[ix,1]),
                                       SE=as.numeric(beta[ix,2]))
        mu_gL <- as.numeric(kronecker(as.matrix(geno@coeff[id,fix.eff.marker]),diag(n.loc))%*%beta[ix,1])
        vars[1,1,"fixed.marker"] <- pvar(mu=mu_gL,weights=gL.weights) 
      }
    }
  } else {
    #multi-trait 
    
    beta <- sans$coef.fixed
    ix <- grep("env_",rownames(beta),fixed=T)
    rownames(beta) <- gsub("env_","",rownames(beta))
    rownames(beta) <- gsub("Trait_","",rownames(beta))
    tmp <- strsplit(rownames(beta)[ix],split=":",fixed=T)
    out$params$env <- data.frame(env=sapply(tmp,"[[",1),
                                 trait=sapply(tmp,"[[",2),
                                 estimate=as.numeric(beta[ix,1]),
                                 SE=as.numeric(beta[ix,2]))
    
    weights <- as.numeric(table(data$env[data$trait==traits[1]]))
    weights <- weights/sum(weights)
    for (i in 1:n.trait) {
       for (j in i:n.trait) {
        mu1 <- out$params$env$estimate[out$params$env$trait==traits[i]]
        mu2 <- out$params$env$estimate[out$params$env$trait==traits[j]]
        vars[i,j,"env"] <- sum(mu1*mu2*weights) - sum(mu1*weights)*sum(mu2*weights)
      }
    }

    if (non.add=="dom") {
      gamma <- (geno@ploidy/2 - 1)/(geno@ploidy - 1)
      ix <- grep("dom",rownames(beta),fixed=T)
      out$params$heterosis <- data.frame(trait=traits,
                                         estimate=as.numeric(beta[ix,1]),
                                         SE=as.numeric(beta[ix,2]))
      
      weights <- id.weights/sum(id.weights)
      for (i in 1:n.trait) 
        for (j in i:n.trait) {
          mu1d <- beta[ix[i],1]*dom.covariate
          mu2d <- beta[ix[j],1]*dom.covariate
          vars[i,j,"heterosis"] <- sum(mu1d*mu2d*weights) - sum(mu1d*weights)*sum(mu2d*weights)
          B[j,i] <- B[i,j] <- gamma^2*(mean(mu1d*mu2d) - mean(mu1d)*mean(mu2d))
        }
    } else {
      gamma <- 0
    }

    if (n.mark > 0) {
      trait.marker <- expand.grid(trait=traits,marker=fix.eff.marker,stringsAsFactors = FALSE)
      ix <- match(apply(trait.marker,1,paste,collapse=":"),rownames(beta))
      x <- out$params$marker <- data.frame(marker=trait.marker$marker,
                                       trait=trait.marker$trait,
                                       estimate=as.numeric(beta[ix,1]),
                                       SE=as.numeric(beta[ix,2]))
      weights <- id.weights/sum(id.weights)
      for (i in 1:n.trait) {
        for (j in i:n.trait) {
          b1 <- matrix(x[x$trait==traits[i],"estimate"],ncol=1)
          mu1 <- as.numeric(as.matrix(geno@coeff[id,fix.eff.marker])%*%b1)
          b2 <- matrix(x[x$trait==traits[j],"estimate"],ncol=1)
          mu2 <- as.numeric(as.matrix(geno@coeff[id,fix.eff.marker])%*%b2)
          vars[i,j,"fixed.marker"] <- sum(mu1*mu2*weights) - sum(mu1*weights)*sum(mu2*weights)
          if (non.add=="dom") {
            mu1 <- mu1+gamma*mu1d
            mu2 <- mu2+gamma*mu2d
          }
          B[j,i] <- B[i,j] <- mean(mu1*mu2) - mean(mu1)*mean(mu2)
        }
      }
    }
  }
    
  #variances 
  vc <- sans$varcomp
  vc <- vc[-which(vc$bound=="F" & round(vc$component)==1L),]
  vc.names <- rownames(vc)
  name2 <- vc.names
  name2 <- gsub("vm(id, source = asremlG, singG = \"PSD\")","additive",name2,fixed=T)
  name2 <- gsub("vm(id, source = asremlD, singG = \"PSD\")","dominance",name2,fixed=T)
  name2 <- gsub("units!units","residual",name2,fixed=T)
  name2 <- gsub("units","residual",name2,fixed=T)
  name2 <- gsub("env.id","residual",name2,fixed=T)
  name2 <- gsub("Trait!Trait!","",name2,fixed=T)
  name2 <- gsub("Trait!Trait_","",name2,fixed=T)
  out$params$vc <- data.frame(name=name2, 
                              estimate=vc$component, SE=vc$std.error)
  
  Imat <- Diagonal(n=n,x=1)
  dimnames(Imat) <- list(id,id)
  
  if (n.trait==1) {
    resid.vc <- Matrix(vc[grep("units",vc.names,fixed=T),1],ncol=1)
    
    if (n.loc > 1) {
      rownames(resid.vc) <- locations
      if (!is.null(vcov)) {
        vars[1,1,"g x env"] <- pvar(V=diag(as.numeric(resid.vc)),weights=loc.weights)
      } else {
        vars[1,1,"residual"] <- pvar(V=diag(as.numeric(resid.vc)),weights=loc.weights)
      }
      
      if (is.null(geno)) {
        if (n.loc > 2) {
          iz <- grep("id:fa",vc.names,fixed=T)
        } else {
          iz <- grep("id:loc",vc.names,fixed=T)
        }
        cov.ans <- f.cov.loc(vc=vc[iz,],locations)
        geno1.vc <- cov.ans$cov.mat
        geno2.vc <- Matrix(NA,nrow=0,ncol=0)
        model <- 0L
        vars[1,1,"genotype"] <- mean(geno1.vc[upper.tri(geno1.vc,diag=FALSE)])*(1 - 1/n)
        
        K <- kronecker(Imat,geno1.vc,make.dimnames = T)
        vars[1,1,"g x loc"] <- pvar(V=K,weights=gL.weights) - vars[1,1,"genotype"]
        
      } else {
        cov.ans <- f.cov.loc(vc[grep("source = asremlG",vc.names,fixed=T),],locations)
        geno1.vc <- cov.ans$cov.mat
        geno2.vc <- Matrix(NA,nrow=0,ncol=0)
        vars[1,1,"additive"] <- mean(geno1.vc[upper.tri(geno1.vc,diag=FALSE)]) * meanG
        K <- kronecker(.GlobalEnv$asremlG,geno1.vc,make.dimnames = T)
        vars[1,1,"add x loc"] <- pvar(V=K,weights=gL.weights) - vars[1,1,"additive"]
        model <- 1L
        
        if (non.add=="g.resid") {
          tmp <- Matrix(vc[grep("id:loc",vc.names,fixed=T),1],ncol=1)
          rownames(tmp) <- locations
          geno2.vc <- coerce_dpo(tcrossprod(sqrt(tmp)))
          K <- kronecker(Imat,geno2.vc,make.dimnames = T)
          vars[1,1,"g.resid"] <- pvar(V=K,weights=gL.weights)
          model <- 2L
        }
        if (non.add=="dom") {
          #tmp <- Matrix(vc[grep("source = asremlD):loc",vc.names,fixed=T),1],ncol=1)
          tmp <- Matrix(vc[grep("source = asremlD",vc.names,fixed=T),1],ncol=1)
          rownames(tmp) <- locations
          geno2.vc <- coerce_dpo(tcrossprod(sqrt(tmp)))
          K <- kronecker(.GlobalEnv$asremlD,geno2.vc,make.dimnames = T)
          vars[1,1,"dominance"] <- pvar(V=K,weights=gL.weights)
          model <- 3L
        }
        
      }
      out <- c(out,loadings=list(cov.ans$loadings))
      
    } else {
      
      #one location
      if (!is.null(vcov)) {
        vars[1,1,"g x env"] <- as.numeric(resid.vc)*(1 - 1/n.gE)
      } else {
        vars[1,1,"residual"] <- as.numeric(resid.vc)*(1 - 1/n.gE)
      }
      
      if (is.null(geno)) {
        model <- 0L
        geno1.vc <- Matrix(vc["id",1],ncol=1)
        geno2.vc <- Matrix(NA,nrow=0,ncol=0)
        vars[1,1,"genotype"] <- as.numeric(geno1.vc)*(1 - 1/n)
      } else {
        model <- 1L
        geno1.vc <- Matrix(vc[grep("source = asremlG",vc.names,fixed=T),1],ncol=1)
        geno2.vc <- Matrix(NA,nrow=0,ncol=0)
        vars[1,1,"additive"] <- as.numeric(geno1.vc)*meanG
        if (non.add=="g.resid") {
          geno2.vc <- Matrix(vc["id",1],ncol=1)
          model <- 2L
          vars[1,1,"g.resid"] <- as.numeric(geno2.vc)*(1 - 1/n)
        }
        if (non.add=="dom") {
          geno2.vc <- Matrix(vc[grep("source = asremlD",vc.names,fixed=T),1],ncol=1)
          model <- 3L
          vars[1,1,"dominance"] <- as.numeric(geno2.vc)*meanD
        }
      }
    }
  } else {
    #multi-trait
    
    iv <- grep("env.id:Trait!",vc.names,fixed=T)
    resid.vc <- f.cov.trait(vc[iv,],traits,us=TRUE)
    if (!is.null(vcov)) {
      tmp <- "g x env"
    } else {
      tmp <- "residual"
    }
    
    for (i in 1:n.trait)
      for (j in i:n.trait) 
        vars[i,j,tmp] <- resid.vc[i,j]*(1 - 1/n.gE)
    
    vc <- vc[-iv,]
    if (is.null(geno)) {
      model <- 0L
      geno2.vc <- Matrix(NA,nrow=0,ncol=0)
      geno1.vc <- f.cov.trait(vc[grep("id:Trait!",rownames(vc),fixed=T),],
                              traits,us=(n.trait>2))
      for (i in 1:n.trait)
        for (j in i:n.trait)
          vars[i,j,"genotype"] <- geno1.vc[i,j]*(1 - 1/n)
      
    } else {
      model <- 1L
      geno2.vc <- Matrix(NA,nrow=0,ncol=0)
      geno1.vc <- f.cov.trait(vc[grep("source = asremlG",rownames(vc),fixed=T),],
                            traits,us=(n.trait>2))
      for (i in 1:n.trait) 
        for (j in i:n.trait) 
          vars[i,j,"additive"] <- geno1.vc[i,j]*meanG
      
      if (non.add=="g.resid") {
        geno2.vc <- f.cov.trait(vc[grep("id:Trait!",rownames(vc),fixed=T),],
                              traits,us=(n.trait>2))
        model <- 2L
        for (i in 1:n.trait) 
          for (j in i:n.trait) 
            vars[i,j,"g.resid"] <- geno2.vc[i,j]*(1 - 1/n)
        
      } 
      if (non.add=="dom") {
        geno2.vc <- f.cov.trait(vc[grep("source = asremlD",rownames(vc),fixed=T),],
                            traits,us=(n.trait>2))
        model <- 3L
        for (i in 1:n.trait) 
          for (j in i:n.trait) 
            vars[i,j,"dominance"] <- geno2.vc[i,j]*meanD
        
      }
    }
  }
  
  if (n.mark==0)
    fix.eff.marker <- character(0)
  
  if (!is.null(covariates))
    attributes(fix.eff.marker) <- list(covariates=covariates)
  
  if (n.trait > 1) {
    for (i in 1:12) 
      vars[,,i][lower.tri(vars[,,i])] <- vars[,,i][upper.tri(vars[,,i])]
  }
  
  out$vars <- new(Class="class_var",geno1=geno1.vc,geno2=geno2.vc,model=model,
                  resid=resid.vc,diagG=diagG,diagD=diagD,
                  vars=vars,B=B,fix.eff.marker=fix.eff.marker)
  
  #random effects
  if (n.trait==1) {
    u <- sans$coef.random 
    if (n.loc > 1) {
      id.loc.names <- expand.grid(loc=locations,id=id,stringsAsFactors = F)[,c(2,1)]
      if (!is.null(geno)) {
        if (n.loc==2) {
          id.loc <- expand.grid(paste0("loc_",locations),
                                paste0("vm(id, source = asremlG, singG = \"PSD\")_",id),
                                stringsAsFactors = F)[,c(2,1)]
        } else {
          id.loc <- expand.grid(paste0("fa(loc, 2)_",locations),
                                paste0("vm(id, source = asremlG, singG = \"PSD\")_",id),
                                stringsAsFactors = F)[,c(2,1)]
        }
        ix1 <- match(apply(id.loc,1,paste,collapse=":"),rownames(u))
        out$random <- data.frame(id.loc.names,add=as.numeric(u[ix1,1]))
        
        if (non.add=="g.resid") {
          id.loc <- expand.grid(paste0("loc_",locations),
                                paste0("id_",id),stringsAsFactors = F)[,c(2,1)]
          ix2 <- match(apply(id.loc,1,paste,collapse=":"),rownames(u))
          out$random$g.iid=as.numeric(u[ix2,1])
        } 
        if (non.add=="dom") {
          id.loc <- expand.grid(paste0("loc_",locations),
                                paste0("vm(id, source = asremlD)_",id),stringsAsFactors = F)[,c(2,1)]
          ix2 <- match(apply(id.loc,1,paste,collapse=":"),rownames(u))
          out$random$dom=as.numeric(u[ix2,1])
        }
      } else {
        if (n.loc==2) {
          id.loc <- expand.grid(paste0("loc_",locations),paste0("id_",id),stringsAsFactors = F)[,c(2,1)]
        } else {
          id.loc <- expand.grid(paste0("fa(loc, 2)_",locations),paste0("id_",id),stringsAsFactors = F)[,c(2,1)]
        }
        ix <- match(apply(id.loc,1,paste,collapse=":"),rownames(u))
        out$random <- data.frame(id.loc.names,g.iid=as.numeric(u[ix,1]))
      }
    } else {
      if (!is.null(geno)) {
        ix1 <- match(paste("vm(id, source = asremlG, singG = \"PSD\")",id,sep="_"),rownames(u))
        out$random <- data.frame(id=id,add=as.numeric(u[ix1,1]))
        if (non.add=="g.resid") {
          ix2 <- match(paste("id",id,sep="_"),rownames(u))
          out$random$g.iid=as.numeric(u[ix2,1])
        }
        if (non.add=="dom") {
          ix2 <- match(paste("vm(id, source = asremlD)",id,sep="_"),rownames(u))
          out$random$dom=as.numeric(u[ix2,1])
        }
      } else {
        ix <- match(paste("id",id,sep="_"),rownames(u))
        out$random <- data.frame(id=id,g.iid=as.numeric(u[ix,1]))
      }
    }
  } else {
    #multi-trait
    u <- sans$coef.random 
    id.trait.names <- expand.grid(trait=traits,id=id,stringsAsFactors = F)[,c(2,1)]
    
    if (is.null(geno)) {
      id.trait <- expand.grid(paste0("Trait_",traits),paste0("id_",id),stringsAsFactors = F)[,c(2,1)]
      ix <- match(apply(id.trait,1,paste,collapse=":"),rownames(u))
      out$random <- data.frame(id.trait.names,g.iid=as.numeric(u[ix,1]))
    } else {
      id.trait <- expand.grid(paste0("Trait_",traits),paste0("vm(id, source = asremlG, singG = \"PSD\")_",id),
                              stringsAsFactors = F)[,c(2,1)]
      ix1 <- match(apply(id.trait,1,paste,collapse=":"),rownames(u))
      out$random <- data.frame(id.trait.names,add=as.numeric(u[ix1,1]))
      if (non.add=="g.resid") {
        id.trait <- expand.grid(paste0("Trait_",traits),paste0("id_",id),stringsAsFactors = F)[,c(2,1)]
        ix2 <- match(apply(id.trait,1,paste,collapse=":"),rownames(u))
        out$random$g.iid=as.numeric(u[ix2,1])
      }
      if (non.add=="dom") {
        id.trait <- expand.grid(paste0("Trait_",traits),paste0("vm(id, source = asremlD)_",id),
                                stringsAsFactors = F)[,c(2,1)]
        ix2 <- match(apply(id.trait,1,paste,collapse=":"),rownames(u))
        out$random$dom=as.numeric(u[ix2,1])
      }
    }
  }
  
  if (!is.null(geno))
    rm("asremlG",envir = .GlobalEnv)
  if (!is.null(vcov))
    rm("asremlOmega",envir = .GlobalEnv)
  if (non.add=="dom")
    rm("asremlD",envir = .GlobalEnv)
  
  return(out)
}

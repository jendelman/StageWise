#' Stage 1 analysis of multi-environment trials
#' 
#' Computes genotype BLUEs for each experiment
#' 
#' The input file must have one column labeled "id" for the individuals and one labeled "env" for the environments. The data for each environment are analyzed independently with a linear mixed model. Although not used in Stage1, to include a genotype x location effect in \code{\link{Stage2}}, a column labeled "loc" should be present in the input file. 
#' 
#' Argument \code{effects} is used to specify other i.i.d. effects besides genotype and has three columns: name, fixed, factor. The "name" column is a string that must match a column in the input file. The fixed column is a logical variable to indicate whether the effect is fixed (TRUE) or random (FALSE). The factor column is a logical variable to indicate whether the effect is a factor (TRUE) or numeric (FALSE). 
#' 
#' Argument \code{solver} specifies which software to use for REML. Current options are "asreml" and "spats". For "spats", the argument \code{spline} must be a vector of length two, with the names of the x and y variables (respectively) for the 2D spline.
#' 
#' The heritability and residuals in the output are based on a random effects model for id.
#' 
#' Missing response values are omitted for single-trait analysis but retained for multi-trait analysis (unless both traits are missing), to allow for prediction in Stage 2. 
#' 
#' Argument \code{workspace} is a vector of length two containing the workspace and pworkspace limits for ASReml-R, with default values of 500mb. If you get an error about insufficient memory, try increasing the appropriate value (workspace for variance estimation and pworkspace for BLUE computation).
#' 
#' For multiple traits, only "asreml" is supported, and only the BLUE model is run, so the returned object does not contain H2. 
#' 
#' If the input file has a column "expt", this allows for the use of separate spatial models for multiple experiments within an environment (only for single trait): each experiment is first analyzed separately, and then the BLUEs from all experiments per env are jointly analyzed to compute a single BLUE per env. The estimation errors from each experiment are propagated into the multi-expt model using ASReml-R. The situation is different with multi-trait analysis, as all experiments are analyzed jointly per env, with a fixed effect for expt but a common residual model. Any additional cofactors (e.g., block) that are nested within expt need to be explicitly nested! 
#' 
#' @param filename Name of CSV file
#' @param traits trait names (see Details)
#' @param effects data frame specifying other effects in the model (see Details)
#' @param solver one of the following: "asreml","spats"
#' @param spline vector of variable names for 2D spline with SpATS
#' @param silent TRUE/FALSE, whether to suppress REML output
#' @param workspace memory limits for ASRreml-R
#' 
#' @return List containing
#' \describe{
#' \item{blues}{data frame of BLUEs}
#' \item{vcov}{list of variance-covariance matrices for the BLUEs, one per experiment (env)}
#' \item{fit}{data frame with broad-sense H2 (plot basis) and/or AIC}
#' \item{resid}{For single trait, list of diagnostic plots and data frame of residuals. For multi-trait, list of resid var-cov matrices.}
#' }
#' 
#' @importFrom utils combn read.csv
#' @import ggplot2
#' @import SpATS
#' @importFrom ggpubr ggarrange
#' @importFrom rlang .data
#' @importFrom stats resid
#' @importFrom spam as.dgCMatrix.spam
#' @import Matrix
#' @export

Stage1 <- function(filename,traits,effects=NULL,solver="asreml",
                   spline=NULL,silent=TRUE,workspace=c("500mb","500mb")) {
  
  data <- read.csv(file=filename,check.names=F)
  solver <- toupper(solver)
  stopifnot(solver %in% c("ASREML","SPATS"))
  stopifnot(traits %in% colnames(data))
  stopifnot(effects$name %in% colnames(data))
  n.trait <- length(traits)
  
  stopifnot(requireNamespace("asreml"))
  library(asreml)
  asreml::asreml.options(maxit=30,workspace=workspace[1],pworkspace=workspace[2],trace=FALSE)
  
  if (solver=="SPATS") {
    if (n.trait > 1) {
      stop("Use asreml for multiple traits")
    }
    stopifnot(requireNamespace("SpATS"))
    stopifnot(spline %in% colnames(data))
  }
  
  if (!is.null(effects)) {
    stopifnot(effects$names %in% colnames(data))
  }
  stopifnot(c("id","env") %in% colnames(data))
  
  if (!is.element("expt",colnames(data))) {
    data$expt <- data$env
    no.expt <- TRUE
  } else {
    no.expt <- FALSE
    tmp <- split(data$env,data$expt)
    tmp2 <- sapply(tmp,function(z){length(unique(z))})
    if (any(tmp2 > 1)) {
      stop("expt names must be unique within env")
    }
  }
  
  data$env <- as.character(data$env)
  data$expt <- as.character(data$expt)
  data$id <- as.character(data$id)
  expt.og <- unique(data$expt)
  
  iz <- apply(as.matrix(data[,traits]),1,function(z){!all(is.na(z))})
  data <- data[iz,]
  
  tmp <- split(data$id,data$expt)
  replicated <- sapply(tmp,function(x){
    y <- table(table(x))
    any(as.integer(names(y)) > 1)
  })
  expts <- names(which(replicated))
  data <- data[data$expt %in% expts,]
  
  if (nrow(data)==0) {
    stop("No experiments with replication")
  }
  envs <- unique(data$env)
  expt.missing <- setdiff(expt.og,expts)
  
  if (length(expt.missing)>1) {
    cat("Some experiments removed due to missing data or lack of replication:\n")
    cat(paste0(paste(expt.missing,collapse="\n"),"\n"))
  }
  n.expt <- length(expts)
  n.env <- length(envs)
  
  effect.table <- matrix("",nrow=2,ncol=2)
  rownames(effect.table) <- c("blue","blup")
  colnames(effect.table) <- c("fixed","random")
  if (n.trait == 1) {
    if (solver=="ASREML") {
      model <- sub("traits",traits,
                 "asreml(data=data1,na.action=na.method(x='omit'),fixed=traits~FIX,random=~RANDOM,residual=~idv(units))",fixed=T)
      effect.table[1,1] <- "id"
      effect.table[2,2] <- "id"
    }
    if (solver=="SPATS") 
      model <- sub("traits",traits,
                   "SpATS(data=data1,response='traits',genotype='id',fixed=~FIX,random=~RANDOM,spatial=~SAP(spline.x,spline.y),genotype.as.random=GARgar",fixed=T)
  } else {
    model <- sub("response",paste(traits,collapse=","),
                 "asreml(data=data1,na.action=na.method(x='omit'),fixed=cbind(response)~FIX,random=~RANDOM,residual=~id(units):us(trait))",fixed=T)
    effect.table[1,1] <- "id:trait"
  }

  factor.vars <- "id"

  if (!is.null(effects)) {
    if (n.trait > 1) {
      effects$name2 <- apply(array(effects$name),1,paste,"trait",sep=":")
    } else {
      effects$name2 <- effects$name
    }
    for (i in 1:nrow(effects)) {
      if (effects$fixed[i]) {
        effect.table[,1] <- paste(effect.table[,1],effects$name2[i],sep="+")
      } else {
        effect.table[,2] <- paste(effect.table[,2],effects$name2[i],sep="+")
      }
    }
    factor.vars <- c(factor.vars,effects$name[which(effects$factor)])
    numeric.vars <- effects$name[which(!effects$factor)]
  } else {
    numeric.vars <- character(0)
  }
  n.numeric <- length(numeric.vars)
  
  #eliminate leading "+"
  effect.table <- apply(effect.table,c(1,2),function(z){if(substr(z,1,1)=="+"){substr(z,2,nchar(z))}else{z}})
  
  #BLUE model
  if (effect.table[1,2]=="") {
    blue.model <- sub("random=~RANDOM,","",model,fixed=T)
  } else {
    blue.model <- sub("RANDOM",effect.table[1,2],model,fixed=T)
  }
  if (effect.table[1,1]=="") {
    blue.model <- sub("fixed=~FIX,","",blue.model,fixed=T)
  } else {
    blue.model <- sub("FIX",effect.table[1,1],blue.model,fixed=T)
  }
  
  #BLUP model
  if (solver=="ASREML" & effect.table[2,1]=="") {
    effect.table[2,1] <- "1"
  }
  if (effect.table[2,2]=="") {
    blup.model <- sub("random=~RANDOM,","",model,fixed=T)
  } else {
    blup.model <- sub("RANDOM",effect.table[2,2],model,fixed=T)
  }
  if (effect.table[2,1]=="") {
    blup.model <- sub("fixed=~FIX,","",blup.model,fixed=T)
  } else {
    blup.model <- sub("FIX",effect.table[2,1],blup.model,fixed=T)
  }
  if (solver=="SPATS") {
    blup.model <- sub("GARgar","TRUE",blup.model,fixed=T)
    blue.model <- sub("GARgar","FALSE",blue.model,fixed=T)
    blup.model <- sub("spline.x",spline[1],blup.model,fixed=T)
    blup.model <- sub("spline.y",spline[2],blup.model,fixed=T)
    blue.model <- sub("spline.x",spline[1],blue.model,fixed=T)
    blue.model <- sub("spline.y",spline[2],blue.model,fixed=T)
    if (silent) {
      blup.model <- paste0(blup.model,",control=list(monitoring=0))")
      blue.model <- paste0(blue.model,",control=list(monitoring=0))")
    } else {
      blup.model <- paste0(blup.model,")")
      blue.model <- paste0(blue.model,")")
    }
  }
  
  resid.blup <- NULL
  
  if ("loc" %in% colnames(data)) {
    fit <- data[!duplicated(data$expt),c("loc","env","expt")]
  } else {
    fit <- data[!duplicated(data$expt),c("env","expt")]
  }
  
  fit <- fit[match(expts,fit$expt),]
  if (n.trait==1) {
    fit$H2 <- as.numeric(NA)
    if (solver=="ASREML") 
      fit$AIC <- as.numeric(NA)
  } else {
    fit$AIC <- as.numeric(NA)
  }
  if (no.expt)
    fit <- fit[,-match("expt",colnames(fit))]

  blue.out <- NULL
  blup.resid <- NULL
  resid.vc <- vector("list",n.expt)
  names(resid.vc) <- expts
  if (solver=="SPATS") {
    spatial.plot <- vector("list",n.expt)
    names(spatial.plot) <- expts
  }
  if (!silent)
    cat(sub("X",paste(traits,collapse=" "),"Traits: X\n"))
  
  if (n.trait==1) {
    vcov <- vector("list",n.expt)
    names(vcov) <- expts
    
    for (j in 1:n.expt) {
      if (!silent)
        cat(sub("X",expts[j],"expt: X\n"))
    
      ix <- which(data$expt==expts[j])
      data1 <- data[ix,]
    
      for (q in 1:length(factor.vars)) {
        eval(parse(text="data1[,factor.vars[q]] <- factor(as.character(data1[,factor.vars[q]]))"))
      }
      if (n.numeric > 0) {
        for (q in 1:n.numeric) {
          eval(parse(text="data1[,numeric.vars[q]] <- as.numeric(data1[,numeric.vars[q]])"))
        }
      }
    
      #BLUP model
      ans <- try(eval(parse(text=blup.model)),silent=TRUE)
      if ((is(ans,"try-error")) || ((solver=="ASREML")&&(!ans$converge))) {
        cat("BLUP model failed to converge.\n")
        next
      }
      residuals <- resid(ans)
      blup.resid <- rbind(blup.resid,
                          data.frame(id=as.character(data1$id),expt=expts[j],resid=residuals))
      if (solver=="ASREML") {
        vc <- summary(ans)$varcomp
        Vg <- vc[match("id",rownames(vc)),1]
        Ve <- vc[match("units!units",rownames(vc)),1]
        fit$H2[j] <- round(Vg/(Vg+Ve),2)
      }
      if (solver=="SPATS") {
        fit$H2[j] <- round(as.numeric(getHeritability(ans)),2)
        x.coord <- ans$data[,ans$terms$spatial$terms.formula$x.coord]
        y.coord <- ans$data[,ans$terms$spatial$terms.formula$y.coord]
        fit.spatial.trend <- obtain.spatialtrend(ans)
        p1.data <- data.frame(x=x.coord,y=y.coord,z=residuals)
        p1 <- ggplot(p1.data,aes(x=.data$x,y=.data$y,fill=.data$z)) + geom_tile() + scale_fill_viridis_c(name="") + xlab(spline[1]) + ylab(spline[2]) + ggtitle("Residuals")
        p2.data <-data.frame(expand.grid(x=fit.spatial.trend$col.p, y=fit.spatial.trend$row.p), z=as.numeric(t(fit.spatial.trend$fit)))
        p2 <- ggplot(p2.data,aes(x=.data$x,y=.data$y,fill=.data$z)) + geom_tile() + scale_fill_viridis_c(name="") + xlab(spline[1]) + theme(axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank()) + ggtitle("Spatial Trend")
        spatial.plot[[j]] <- ggarrange(p1,p2,common.legend = TRUE,legend = "right")
      }
  
      #BLUE model
      asreml::asreml.options(trace=!silent)
      ans <- try(eval(parse(text=blue.model)),silent=TRUE)
      if ((is(ans,"try-error")) || ((solver=="ASREML")&&(!ans$converge))) {
        cat("BLUE model failed to converge.\n")
      } else {
        asreml::asreml.options(trace=FALSE)
      
        if (solver=="ASREML") {
          predans <- asreml::predict.asreml(ans,classify="id",vcov = TRUE)
          tmp <- predans$pvals[,c("id","predicted.value")]
          colnames(tmp) <- c("id","BLUE")
          vcov[[j]] <- predans$vcov
          fit$AIC[j] <- round(as.numeric(summary(ans)$aic),1)
        }
        if (solver=="SPATS") {
          predans <- predict.SpATS(object=ans,which="id",predFixed="marginal",
                                   return.vcov.matrix = TRUE)
          tmp <- predans[,c("id","predicted.values")]
          colnames(tmp) <- c("id","BLUE")
          tmp2 <- Matrix(spam::as.dgCMatrix.spam(attr(predans,"vcov")))
          vcov[[j]] <- as(forceSymmetric(tmp2),"packedMatrix")
        }
        tmp$id <- as.character(tmp$id)
        dimnames(vcov[[j]]) <- list(tmp$id,tmp$id)
        blue.out <- rbind(blue.out,data.frame(expt=expts[j],tmp))
      }
    }
  
    ik <- which(!sapply(vcov,is.null))
    vcov <- vcov[ik]
    tmp <- data[data$expt %in% names(ik),]
    expt.in.env <- lapply(split(tmp$expt,factor(tmp$env,levels=envs)),unique)
    blue.out$env <- data$env[match(blue.out$expt,data$expt)]
  
    nee <- sapply(expt.in.env,length)
    iu <- which(nee==1)
    if (length(iu) > 0) {
      tmp <- names(vcov)[iu]
      names(vcov)[iu] <- blue.out$env[match(tmp,blue.out$expt)]
    }
  
    for (j in which(nee > 1)) {
      omega.list <- vcov[expt.in.env[[j]]]
      vcov2 <- mapply(FUN=function(x,y){
        tmp <- paste(y,rownames(x),sep=":")
        dimnames(x) <- list(tmp,tmp)
        return(x)},x=omega.list,y=as.list(expt.in.env[[j]]))
    
      .GlobalEnv$asremlOmega <- direct_sum(lapply(vcov2,solve))
      dname <- lapply(vcov2,rownames)
      dimnames(.GlobalEnv$asremlOmega) <- list(unlist(dname),unlist(dname))
      attr(.GlobalEnv$asremlOmega,"INVERSE") <- TRUE
    
      ix <- which(blue.out$env==envs[j])
      data2 <- blue.out[ix,]
      data2$expt.id <- factor(paste(data2$expt,data2$id,sep=":"))
      data2$id <- factor(data2$id)
      data2$expt <- factor(data2$expt)
      start.table <- asreml::asreml(data=data2,fixed=BLUE~expt-1+id,
                            random=~vm(expt.id,source=asremlOmega),
                            residual=~idv(units),start.values=TRUE)$vparameters.table
      k <- grep("Omega",start.table$Component,fixed=T)
      start.table$Value[k] <- 1
      start.table$Constraint[k] <- "F"
      ans3 <- asreml::asreml(data=data2,fixed=BLUE~expt-1+id,
                     random=~vm(expt.id,source=asremlOmega),
                     residual=~idv(units),G.param=start.table)
      
      predans3 <- asreml::predict.asreml(ans3,classify="id",vcov = TRUE)
      blue3 <- data.frame(expt=NA,predans3$pvals[,c("id","predicted.value")],env=envs[j])
      colnames(blue3) <- c("expt","id","BLUE","env")
      vcov3 <- predans3$vcov
      dimnames(vcov3) <- list(blue3$id,blue3$id)
      vcov2 <- vcov[-match(expt.in.env[[j]],names(vcov))]
      vcov <- c(vcov2,vcov3)
      names(vcov) <- c(names(vcov2),envs[j])
      blue.out <- rbind(blue.out[-ix,],blue3)
      rm("asremlOmega",envir = .GlobalEnv)
    }
  }
  
  if (n.trait > 1) {
    vcov <- vector("list",n.env)
    names(vcov) <- envs
    
    expt.in.env <- lapply(split(data$expt,factor(data$env,levels=envs)),unique)
    nexpt.per.env <- sapply(expt.in.env,length)
    #analysis by environment, not experiment
    for (j in 1:n.env) {   
      if (!silent)
        cat(sub("X",envs[j],"env: X\n"))
      
      ix <- which(data$env==envs[j])
      data1 <- data[ix,]
      
      for (q in 1:length(factor.vars)) {
        eval(parse(text="data1[,factor.vars[q]] <- factor(as.character(data1[,factor.vars[q]]))"))
      }
      if (n.numeric > 0) {
        for (q in 1:n.numeric) {
          eval(parse(text="data1[,numeric.vars[q]] <- as.numeric(data1[,numeric.vars[q]])"))
        }
      }
      
      if (nexpt.per.env[j] > 1) {
        blue.model2 <- sub("id:trait","id:trait+expt:trait",blue.model)
        data1$expt <- factor(data1$expt)
        #blue.model2 <- sub("id(units):us(trait)","dsum(~units:us(trait)|expt)",blue.model2,fixed=T)
      } else {
        blue.model2 <- blue.model
      }
      asreml::asreml.options(trace=!silent)
      ans <- try(eval(parse(text=blue.model2)),silent=TRUE)
      if ((is(ans,"try-error")) || ((solver=="ASREML")&&(!ans$converge))) {
        cat("BLUE model failed to converge.\n")
      } else {
        asreml::asreml.options(trace=FALSE)
        vc <- summary(ans)$varcomp
        vc <- vc[-which(vc$bound=="F" & round(vc$component)==1L),c("component","std.error")]
        vc.names <- rownames(vc)
        iv <- grep("units:trait!",vc.names,fixed=T)
        resid.vc[[j]] <- f.cov.trait(vc[iv,],traits,us=TRUE)
    
        fit$AIC[j] <- round(as.numeric(summary(ans)$aic),1)
        predans <- asreml::predict.asreml(ans,classify="id:trait",vcov = TRUE)
        tmp <- predans$pvals[,c("id","trait","predicted.value")]
        colnames(tmp) <- c("id","trait","BLUE")
        vcov[[j]] <- predans$vcov
        tmp$id <- as.character(tmp$id)
        id.trait <- apply(tmp[,1:2],1,paste,collapse=":")
        dimnames(vcov[[j]]) <- list(id.trait,id.trait)
        blue.out <- rbind(blue.out,data.frame(env=envs[j],tmp))
      }
    }
  }
  
  if (n.trait==1) {
    blue.out <- blue.out[,c("env","id","BLUE")]
  } else {
    blue.out <- blue.out[,c("env","id","trait","BLUE")]
  }
  blue.out <- blue.out[order(blue.out$env,blue.out$id),]
  if ("loc" %in% colnames(data)) {
    blue.out$loc <- as.character(data$loc[match(blue.out$env,data$env)])
  }
  
  if (n.trait==1) {
    p1 <- ggplot(data=blup.resid,aes(y=.data$resid,x=.data$expt)) + ylab("Residual") + xlab("") +
      stat_boxplot(outlier.color="red") + theme_bw() + theme(axis.text.x=element_text(angle=90,vjust=0.5))
    p2 <- ggplot(data=blup.resid,aes(sample=.data$resid)) + stat_qq() + stat_qq_line() + facet_wrap(~expt) + theme_bw() + xlab("Expected") + ylab("Observed")
    if (solver=="ASREML")
      return(list(blues=blue.out,vcov=vcov,fit=fit,
                  resid=list(boxplot=p1,qqplot=p2,table=blup.resid)))
    if (solver=="SPATS")
      return(list(blues=blue.out,vcov=vcov,fit=fit,
                  resid=list(boxplot=p1,qqplot=p2,spatial=spatial.plot,table=blup.resid)))
  } else {
    #Multi-trait
    return(list(blues=blue.out,vcov=vcov,fit=fit,resid=resid.vc))
  }
}

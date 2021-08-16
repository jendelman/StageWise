#' Stage 1 analysis of multi-environment trials
#' 
#' Computes genotype BLUEs within each environment
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
#' \item{blue}{data frame of BLUEs for all environments}
#' \item{vcov}{list of variance-covariance matrices for the BLUEs, one per env}
#' \item{H2}{broad-sense H2 on a plot basis}
#' \item{resid}{list containing diagnostic plots and data frame of residuals}
#' }
#' 
#' @importFrom utils combn read.csv
#' @import ggplot2
#' @import asreml
#' @import SpATS
#' @importFrom ggpubr ggarrange
#' @importFrom rlang .data
#' @importFrom stats resid
#' @importFrom spam as.dgCMatrix.spam
#' @export

Stage1 <- function(filename,traits,effects=NULL,solver="asreml",
                   spline=NULL,silent=TRUE,workspace=c("500mb","500mb")) {
  
  data <- read.csv(file=filename,check.names=F)
  solver <- toupper(solver)
  stopifnot(solver %in% c("ASREML","SPATS"))
  if (solver=="ASREML") {
    stopifnot(requireNamespace("asreml"))
    asreml::asreml.options(maxit=30,workspace=workspace[1],pworkspace=workspace[2],trace=!silent)
  }
  if (solver=="SPATS") {
    stopifnot(requireNamespace("SpATS"))
    stopifnot(spline %in% colnames(data))
  }
  stopifnot(traits %in% colnames(data))
  n.trait <- length(traits)
  if (n.trait > 1) {
    stop("Multiple traits not supported yet. Check back soon.")
  }
  if (!is.null(effects)) {
    stopifnot(effects$names %in% colnames(data))
  }
  stopifnot(c("id","env") %in% colnames(data))
  
  data$env <- as.character(data$env)
  data$id <- as.character(data$id)
  envs <- unique(data$env)
  n.env <- length(envs)
  
  iz <- apply(as.matrix(data[,traits]),1,function(z){!all(is.na(z))})
  data <- data[iz,]
  
  effect.table <- matrix("",nrow=2,ncol=2)
  rownames(effect.table) <- c("blue","blup")
  colnames(effect.table) <- c("fixed","random")
  if (n.trait == 1) {
    if (solver=="ASREML") {
      model <- sub("traits",traits,
                 "asreml(data=data1,na.action=na.method(y='omit',x='omit'),fixed=traits~FIX,random=~RANDOM,residual=~idv(units))",fixed=T)
      effect.table[1,1] <- "id"
      effect.table[2,2] <- "id"
    }
    if (solver=="SPATS") 
      model <- sub("traits",traits,
                   "SpATS(data=data1,response='traits',genotype='id',fixed=~FIX,random=~RANDOM,spatial=~SAP(spline.x,spline.y),genotype.as.random=GARgar",fixed=T)
  } else {
    model <- sub("response",paste(traits,collapse=","),
                 "asreml(data=data1,fixed=cbind(response)~FIX,random=~RANDOM,residual=~id(units):corh(trait))",fixed=T)
    effect.table[1,1] <- "id:trait"
    effect.table[2,] <- NA
  }

  factor.vars <- "id"
  if (!is.null(effects)) {
    for (i in 1:nrow(effects)) {
      if (effects$fixed[i]) {
        effect.table[,1] <- paste(effect.table[,1],effects$name[i],sep="+")
      } else {
        effect.table[,2] <- paste(effect.table[,2],effects$name[i],sep="+")
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
  vcov <- vector("list",n.env)
  names(vcov) <- envs
  if (n.trait==1) {
    H2 <- data.frame(env=envs,H2=as.numeric(NA))
  }
  if ("loc" %in% colnames(data)) {
    H2$loc <- data$loc[match(H2$env,data$env)]
  }

  blue.out <- NULL
  blup.resid <- NULL
  if (solver=="SPATS") {
    spatial.plot <- vector("list",n.env)
    names(spatial.plot) <- envs
  }
  cat(sub("X",paste(traits,collapse=" "),"Traits: X\n"))
  for (j in 1:n.env) {
    cat(sub("X",envs[j],"Env: X\n"))
    ix <- which(data$env==envs[j])
    if (length(ix)==0) {
      cat("Warning: No data present\n")
    } else {
      data1 <- data[ix,]
      for (q in 1:length(factor.vars)) {
        eval(parse(text="data1[,factor.vars[q]] <- factor(as.character(data1[,factor.vars[q]]))"))
      }
      if (n.numeric > 0) {
        for (q in 1:n.numeric) {
          eval(parse(text="data1[,numeric.vars[q]] <- as.numeric(data1[,numeric.vars[q]])"))
        }
      }
    }
    
    ans <- try(eval(parse(text=blue.model)),silent=TRUE)
    if (class(ans)=="try-error") {
      stop("BLUE model failed to converge.")
    }
    if ((solver=="ASREML")&&(!ans$converge)) {
      stop("BLUE model failed to converge.")
    }
    if (n.trait==1) {
      if (solver=="ASREML") {
        predans <- predict.asreml(ans,classify="id",vcov = TRUE)
        tmp <- predans$pvals[,c("id","predicted.value")]
        colnames(tmp) <- c("id","BLUE")
        vcov[[j]] <- predans$vcov
      }
      if (solver=="SPATS") {
        predans <- predict.SpATS(object=ans,which="id",predFixed="marginal",
                                   return.vcov.matrix = TRUE)
        tmp <- predans[,c("id","predicted.values")]
        colnames(tmp) <- c("id","BLUE")
        tmp2 <- Matrix(spam::as.dgCMatrix.spam(attr(predans,"vcov")))
        vcov[[j]] <- as(forceSymmetric(tmp2),"dspMatrix")
      }
      tmp$id <- as.character(tmp$id)
      dimnames(vcov[[j]]) <- list(tmp$id,tmp$id)
      blue.out <- rbind(blue.out,data.frame(env=envs[j],tmp))
      
    } else {
      predans <- predict.asreml(ans,classify="id:trait",vcov = TRUE)
      tmp <- predans$pvals[,c("id","trait","predicted.value")]
      colnames(tmp) <- c("id","trait","BLUE")
    }
    
    ans <- try(eval(parse(text=blup.model)),silent=TRUE)
    if (class(ans)=="try-error") {
      stop("BLUP model failed to converge.")
    }
    if ((solver=="ASREML")&&(!ans$converge)) {
      stop("BLUP model failed to converge.")
    }
    if (n.trait==1) {
      residuals <- resid(ans)
      blup.resid <- rbind(blup.resid,
                          data.frame(id=as.character(data1$id),env=envs[j],resid=residuals))
      if (solver=="ASREML") {
        vc <- summary(ans)$varcomp
        Vg <- vc[match("id",rownames(vc)),1]
        Ve <- vc[match("units!units",rownames(vc)),1]
        H2$H2[j] <- round(Vg/(Vg+Ve),2)
      }
      if (solver=="SPATS") {
        H2$H2[j] <- round(as.numeric(getHeritability(ans)),2)
        x.coord <- ans$data[,ans$terms$spatial$terms.formula$x.coord]
        y.coord <- ans$data[,ans$terms$spatial$terms.formula$y.coord]
        fit.spatial.trend <- obtain.spatialtrend(ans)
        p1.data <- data.frame(x=x.coord,y=y.coord,z=residuals)
        p1 <- ggplot(p1.data,aes(x=.data$x,y=.data$y,fill=.data$z)) + geom_tile() + scale_fill_viridis_c(name="") + xlab(spline[1]) + ylab(spline[2]) + ggtitle("Residuals")
        p2.data <-data.frame(expand.grid(x=fit.spatial.trend$col.p, y=fit.spatial.trend$row.p), z=as.numeric(t(fit.spatial.trend$fit)))
        p2 <- ggplot(p2.data,aes(x=.data$x,y=.data$y,fill=.data$z)) + geom_tile() + scale_fill_viridis_c(name="") + xlab(spline[1]) + theme(axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank()) + ggtitle("Spatial Trend")
        spatial.plot[[j]] <- ggarrange(p1,p2,common.legend = TRUE,legend = "right")
      }
    } else {
      #to do
    } 
  }
  
  if ("loc" %in% colnames(data)) {
    blue.out$loc <- as.character(data$loc[match(blue.out$env,data$env)])
  }
  if (n.trait==1) {
    p1 <- ggplot(data=blup.resid,aes(y=.data$resid,x=.data$env)) + ylab("Residual") + xlab("") +
      stat_boxplot(outlier.color="red") + theme_bw() + theme(axis.text.x=element_text(angle=90,vjust=0.5))
    p2 <- ggplot(data=blup.resid,aes(sample=.data$resid)) + stat_qq() + stat_qq_line() + facet_wrap(~env) + theme_bw() + xlab("Expected") + ylab("Observed")
    if (solver=="ASREML")
      return(list(blue=blue.out,vcov=vcov,H2=H2,
                  resid=list(boxplot=p1,qqplot=p2,table=blup.resid)))
    if (solver=="SPATS")
      return(list(blue=blue.out,vcov=vcov,H2=H2,
                  resid=list(boxplot=p1,qqplot=p2,spatial=spatial.plot,table=blup.resid)))
  } else {
    # to do
  }
}

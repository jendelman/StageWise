#' Stage 1 analysis of multi-environment trials
#' 
#' Computes genotype BLUEs within each environment using ASReml-R
#' 
#' The variable \code{data} must have one column labeled "id" for the individuals and one labeled "env" for the environments. The data for each environment are analyzed independently with a linear mixed model. Although not used in Stage1, to include a genotype x location effect in \code{\link{Stage2}}, a column labeled "loc" should be included in \code{data}. 
#' 
#' Including multiple traits in \code{trait} triggers a multivariate analysis, but for computational reasons, only 2 traits are analyzed at a time. With more than 2 traits, the software analyzes all pairs of traits. For single trait analysis, broad-sense H2 on a plot basis is computed from the variance components for each env, with genotype as a random effect. The residuals from this analysis are also returned as a table and plots.
#' 
#' Argument \code{effects} is used to specify other i.i.d. effects besides genotype and has three columns: name, fixed, factor. The "name" column is a string that must match a column in \code{data}. The fixed column is a logical variable to indicate whether the effect is fixed (TRUE) or random (FALSE). The factor column is a logical variable to indicate whether the effect is a factor (TRUE) or numeric (FALSE). 
#' 
#' Missing response values are omitted for single-trait analysis but retained for multi-trait analysis (unless both traits are missing), to allow for prediction in Stage 2. By default, the workspace and pworkspace limits for ASReml-R are set at 500mb. If you get an error about insufficient memory, try increasing the appropriate value (workspace for variance estimation and pworkspace for BLUE computation).
#' 
#' @references Damesa et al. 2017. Agronomy Journal 109: 845-857. doi:10.2134/agronj2016.07.0395
#' 
#' @param data input data frame (see Details)
#' @param traits trait names (see Details)
#' @param effects data frame specifying other effects in the model (see Details)
#' @param workspace memory limit for ASRreml-R variance estimation
#' @param pworkspace memory limit for ASRreml-R BLUE computation
#' @param silent TRUE/FALSE, whether to suppress ASReml-R output
#' 
#' @return List containing
#' \describe{
#' \item{blue}{data frame of BLUEs for all environments}
#' \item{vcov}{list of variance-covariance matrices for the BLUEs, one per env}
#' \item{H2}{broad-sense H2 on a plot basis (only for single trait)}
#' \item{resid}{list containing boxplot, qqplot, and table of residuals (only for single trait)}
#' }
#' 
#' @importFrom stats complete.cases
#' @importFrom utils combn
#' @import asreml
#' @export

Stage1 <- function(data,traits,effects=NULL,silent=TRUE,workspace="500mb",pworkspace="500mb") {
  
  stopifnot(requireNamespace("asreml"))
  stopifnot(traits %in% colnames(data))
  n.trait <- length(traits)
  if (n.trait > 2) {
    trait2 <- combn(traits,2)
    ans <- apply(trait2,2,function(traits){Stage1(data,traits,effects,silent,workspace,pworkspace)})
    return(ans)
  } else {
  if (!is.null(effects)) {
    stopifnot(effects$names %in% colnames(data))
  }
  stopifnot(c("id","env") %in% colnames(data))
  
  data$env <- as.character(data$env)
  data$id <- as.character(data$id)
  if (n.trait > 1) {
    iz <- apply(data[,traits],1,function(z){!all(is.na(z))})
    data <- data[iz,]
  }
  envs <- unique(data$env)
  n.env <- length(envs)
  
  #prepare asreml command
  asreml.options(maxit=30,workspace=workspace,pworkspace=pworkspace,trace=!silent)
  effect.table <- matrix("",nrow=2,ncol=2)
  rownames(effect.table) <- c("blue","blup")
  colnames(effect.table) <- c("fixed","random")
  if (n.trait > 1) {
    model <- sub("response",paste(traits,collapse=","),
                  "asreml(data=data1,na.action=na.method(y='include',x='omit'),fixed=cbind(response)~FIX,random=~RANDOM,residual=~id(units):corh(trait))",fixed=T)
    effect.table[1,1] <- "id:trait"
    effect.table[2,] <- NA
  } else {
    model <- sub("response",traits,
                 "asreml(data=data1,na.action=na.method(y='omit',x='omit'),fixed=response~FIX,random=~RANDOM,residual=~idv(units))",fixed=T)
    effect.table[1,1] <- "id"
    effect.table[2,2] <- "id"
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
  
  if (effect.table[1,2]=="") {
    blue.model <- sub("random=~RANDOM,","",model,fixed=T)
  } else {
    blue.model <- sub("RANDOM",effect.table[1,2],model,fixed=T)
  }
  blue.model <- sub("FIX",effect.table[1,1],blue.model,fixed=T)
  
  tmp <- effect.table[2,1]
  if (tmp=="") {
    blup.model <- sub("FIX","1",model,fixed=T)
  } else {
    if (substr(tmp,1,1)=="+") {
      tmp <- substr(tmp,2,nchar(tmp))
    } 
    blup.model <- sub("FIX",tmp,model,fixed=T)
  }
  blup.model <- sub("RANDOM",effect.table[2,2],blup.model,fixed=T)
  
  resid.blup <- NULL
  vcov <- vector("list",n.env)
  names(vcov) <- envs
  if (n.trait==1) {
    H2 <- data.frame(env=envs,H2=as.numeric(NA))
  }

  blue.out <- NULL
  blup.resid <- NULL
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
    
    ans <- eval(parse(text=blue.model))
    while (!ans$converge) {
      cat("Fixed effects model failed to converge. Do you wish to continue running? y/n \n")
      input <- readLines(n=1)
      if (input=="y") {
        ans <- update.asreml(ans)
      } else {
        return()
      }
    }
    if (n.trait > 1) {
      predans <- predict.asreml(ans,classify="id:trait",vcov = TRUE)
      tmp <- predans$pvals[,c("id","trait","predicted.value")]
      colnames(tmp) <- c("id","trait","BLUE")
    } else {
      predans <- predict.asreml(ans,classify="id",vcov = TRUE)
      tmp <- predans$pvals[,c("id","predicted.value")]
      colnames(tmp) <- c("id","BLUE")
    }
    tmp$id <- as.character(tmp$id)
    
    blue.out <- rbind(blue.out,data.frame(env=envs[j],tmp))
    vcov[[j]] <- predans$vcov
    rownames(vcov[[j]]) <- colnames(vcov[[j]]) <- paste(tmp$id,tmp$env,sep=":")
    
    if (n.trait==1) {
      ans <- eval(parse(text=blup.model))
      while (!ans$converge) {
        cat("Random effects model failed to converge. Do you wish to continue running? y/n \n")
        input <- readLines(n=1)
        if (input=="y") {
          ans <- update.asreml(ans)
        } else {
          return()
        }
      }
      vc <- summary(ans)$varcomp
      Vg <- vc[match("id",rownames(vc)),1]
      Ve <- vc[match("units!units",rownames(vc)),1]
      H2$H2[j] <- round(Vg/(Vg+Ve),2)
      blup.resid <- rbind(blup.resid,data.frame(id=as.character(data1$id),env=envs[j],
                          resid=resid(ans)))
    } 
  }
  if ("loc" %in% colnames(data)) {
    blue.out$loc <- as.character(data$loc[match(data$env,data$loc)])
  }
  if (n.trait==1) {
    p1 <- ggplot(data=blup.resid,aes(y=resid,x=env)) + ylab("Residual") + xlab("") +
      stat_boxplot(outlier.color="red") + theme_bw() + theme(axis.text.x=element_text(angle=90,vjust=0.5))
    p2 <- ggplot(data=blup.resid,aes(sample=resid)) + stat_qq() + stat_qq_line() + facet_wrap(~env) + theme_bw() + xlab("Expected") + ylab("Observed")
    
    return(list(blue=blue.out,vcov=vcov,resid=list(boxplot=p1,qqplot=p2,table=blup.resid),H2=H2))
  } else {
    return(list(blue=blue.out,vcov=vcov))
  }
  }
}

#' Displays variances and correlations
#'
#' Displays variances and correlations
#' 
#' For a single trait, the 'var' output is a data frame with two columns of information for the various effects: the first is the variance and the second is the proportion of variance explained (PVE), excluding the environment effect. For multiple locations or traits, the 'cor' output is the correlation matrix for additive effects (does not include fixed effect markers). 
#' 
#' For multiple traits, the variance and PVE results are returned as separate data frames, unless \code{index.coeff} is used to create an index. To create a similar output with multiple locations, use \code{separate.loc=TRUE}. 
#' 
#' The \code{index.coeff} are the coefficients of the standardized true values (see also \code{\link{blup}}). The argument \code{gamma} is the relative weight for non-additive to additive genetic merit (see also \code{\link{gain}}). 
#' 
#' @param object object of \code{\link{class_var}}
#' @param digits number of digits for rounding
#' @param index.coeff merit index coefficients
#' @param gamma contribution of non-additive values for genetic merit
#' @param separate.loc TRUE/FALSE controls output for multi-location
#'
#' @return List output that varies depending on the situation (see Details)
#'
#' @include class_var.R
#' @name summary.var
#' @exportMethod summary
NULL

setGeneric("summary")
setMethod("summary",c(object="class_var"),
          definition=function(object,digits=3,index.coeff=NULL,
                              gamma=0,separate.loc=F){
            
            n.trait <- ncol(object@resid)
            n.loc <- ifelse(n.trait > 1,1,nrow(object@resid))
            
            if (n.loc==1 & separate.loc)
              stop("Only one location")
            
            if (dim(object@vars)[1]==1 & n.loc > 1) {
              #older style
              if (separate.loc)
                stop("Update package and rerun Stage2 to utilize this option")
            }
            
            if (n.trait > 1 & !is.null(index.coeff)) {
              traits <- rownames(object@geno1)
              stopifnot(sort(names(index.coeff))==sort(traits))
              index.coeff <- index.coeff[match(traits,names(index.coeff))]
              trait.scale <- sqrt(diag(object@geno1)+gamma^2*diag(object@geno2))
              index.coeff <- index.coeff/trait.scale
            }
            #if (dim(object@vars)[1]==0) {
            #  pvar <- FALSE
            #} else {
            #  pvar <- TRUE
            #}
            pvar <- TRUE
            if (pvar) {
              if (is.null(index.coeff)) {
                if (n.loc==1 | (n.loc > 1 & separate.loc)) {
                  vars <- apply(object@vars,3,diag)
                } else {
                  if (dim(object@vars)[1]==1 & n.loc > 1) {
                    #older style
                    vars <- object@vars[1,1,]
                  } else {
                    vars <- object@vars[1,2,]
                  }
                }
              } else {
                vars <- apply(object@vars,3,function(V){
                  as.numeric(crossprod(index.coeff,V%*%index.coeff))
                })
              }
            }
              
            if ((n.trait > 1 & is.null(index.coeff)) | separate.loc) {
              cor.mat <- cov_to_cor(object@geno1)
              if (pvar) {
                vars <- t(vars)
                prop.var <- t(t(vars[-1,])/apply(vars[-1,],2,sum,na.rm=T))
                prop.var <- round(prop.var[!is.na(prop.var[,1]),],3)
                
                vars <- apply(vars[!is.na(vars[,1]),],2,sigdig,digits=digits)
                return(list(var=data.frame(vars),
                            PVE=data.frame(prop.var),
                            cor.mat=round(cor.mat,3)))
              } else {
                return(round(cor.mat,3))
              }
            } else {
              if (n.loc > 1) {
                cor.mat <- cov_to_cor(object@geno1)
                x <- data.frame(Variance=sigdig(vars,digits),
                                PVE=round(c(NA,vars[-1]/sum(vars[-1],na.rm=T)),3))
                return(list(var=x[!is.na(x[,1]),],cor=round(cor.mat,3)))
              } else {
                x <- data.frame(Variance=sigdig(vars,digits),
                                PVE=round(c(NA,NA,vars[-(1:2)]/sum(vars[-(1:2)],na.rm=T)),3))
                return(x[!is.na(x[,1]),])
              }
            }
          })

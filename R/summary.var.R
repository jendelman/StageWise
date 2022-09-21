#' Displays variances and correlations
#'
#' Displays variances and correlations
#' 
#' For a single trait, the 'var' output is a data frame with two columns of information for the various effects: the first is the variance and the second is the proportion of variance explained (PVE), excluding the environment effect. For multiple locations or traits, the 'cor' output is the correlation matrix for additive effects (does not include fixed effect markers). For multiple traits, the variance and PVE results are returned as separate data frames.
#' 
#' @param object object of \code{\link{class_var}}
#' @param digits number of digits for rounding
#'
#' @return List output that varies depending on the situation (see Details)
#'
#' @include class_var.R
#' @name summary.var
#' @exportMethod summary
NULL

setGeneric("summary")
setMethod("summary",c(object="class_var"),
          definition=function(object,digits=3){
            
            n.trait <- ncol(object@resid)
            vars <- apply(object@vars,3,diag)

            if (n.trait > 1) {
              cor.mat <- cov_to_cor(object@geno1)
              vars <- t(vars)
              prop.var <- t(t(vars[-1,])/apply(vars[-1,],2,sum,na.rm=T))
              prop.var <- round(prop.var[!is.na(prop.var[,1]),],3)
              
              vars <- apply(vars[!is.na(vars[,1]),],2,sigdig,digits=digits)
              return(list(var=data.frame(vars),
                          PVE=data.frame(prop.var),
                          cor.mat=round(cor.mat,3)))
            } else {
              if (nrow(object@geno1) > 1) {
                #multiple locations
                cor.mat <- cov_to_cor(object@geno1)
                x <- data.frame(Variance=sigdig(vars,digits),
                           PVE=round(c(NA,vars[-1]/sum(vars[-1],na.rm=T)),3))
                return(list(var=x[!is.na(x[,1]),],cor=round(cor.mat,3)))
              } else {
                x <- data.frame(Variance=sigdig(vars,digits),
                           PVE=round(c(NA,vars[-1]/sum(vars[-1],na.rm=T)),3))
                return(x[!is.na(x[,1]),])
                }
            }
          })

#' Summarize variances and correlations
#'
#' Summarize variances and correlations
#' 
#' When reporting the partitioning of variance, the variance component for the additive effect is multiplied by the mean diagonal of the G matrix. 
#' 
#' @param object object of \code{\link{class_var}}
#'
#' @return For univariate analysis, a matrix with the proportion of variance. For multi-location or multi-trait analysis, a correlation matrix is also returned. If a G matrix was included in the analysis, the correlation is between additive values; otherwise, it is between genotypic values. If fixed effect markers are present, the additive correlation is above the diagonal, and below the diagonal is the proportion of additive covariance due to the fixed effect markers.
#'
#' @include class_var.R
#' @name summary
#' @exportMethod summary
NULL

setGeneric("summary")
setMethod("summary",c(object="class_var"),
          definition=function(object){
            cor.mat <- NULL
            n.mark <- nrow(object@fixed.marker.var)
            if (ncol(object@resid)>1) {
              #multiple traits
              traits <- rownames(object@resid)
              if (is.na(object@meanG)) {
                vtab <- data.frame(matrix(c(diag(object@g.resid),diag(object@resid)),nrow=2,byrow = T))
                colnames(vtab) <- traits
                if (any(!is.na(object@meanOmega))) {
                  vtab <- rbind(vtab,object@meanOmega)
                  rownames(vtab) <- c("genotype","g x env","residual")
                } else {
                  rownames(vtab) <- c("genotype","residual")
                }
                cor.mat <- round(cov_to_cor(object@g.resid),3)
              } else {
                vtab <- data.frame(matrix(c(diag(object@add)*object@meanG,diag(object@g.resid),
                                            diag(object@resid)),nrow=3,byrow = T))
                colnames(vtab) <- traits
                if (any(!is.na(object@meanOmega))) {
                  vtab <- rbind(vtab,object@meanOmega)
                  rownames(vtab) <- c("additive","g.resid","g x env","residual")
                } else {
                  rownames(vtab) <- c("additive","g.resid","residual")
                }
                
                if (n.mark > 0) {
                  vtab <- rbind(object@fixed.marker.var,vtab)
                  cor.mat <- round(cov_to_cor(object@add+object@fixed.marker.cov),3)
                  tmp <- object@fixed.marker.cov/(object@add+object@fixed.marker.cov)
                  cor.mat[lower.tri(cor.mat)] <- tmp[lower.tri(tmp)]
                  cor.mat <- round(cor.mat,3)    
                } else {
                  cor.mat <- round(cov_to_cor(object@add),3)    
                }
                
              }
              vtab <- as.matrix(vtab)
              pvar <- round(vtab %*% diag(1/apply(vtab,2,sum)),3)
              colnames(pvar) <- traits
            } else {
              
              #single trait
              if (!is.na(object@meanG)) {
                if (nrow(object@add) > 1) {
                  #multiple locations
                  tmp <- partition(object@add)
                  Va <- tmp[1]*object@meanG
                  VaL <- tmp[2]*object@meanG
                  
                  Vr <- mean(object@g.resid[,1])
                  V1 <- matrix(c(Va,VaL,Vr),ncol=1)
                  rownames(V1) <- c("additive","add x loc","g.resid")
                  
                  if (n.mark > 0) {
                    fix.eff.markers <- rownames(object@fixed.marker.var)
                    tmp <- matrix(t(object@fixed.marker.var),ncol=1)
                    rownames(tmp) <- as.vector(matrix(rbind(fix.eff.markers,paste(fix.eff.markers,"x loc")),ncol=1))
                    V1 <- rbind(tmp,V1)
                    
                    cor.mat <- round(cov_to_cor(object@add+object@fixed.marker.cov),3)
                    tmp <- object@fixed.marker.cov/(object@add+object@fixed.marker.cov)
                    cor.mat[lower.tri(cor.mat)] <- tmp[lower.tri(tmp)]
                    cor.mat <- round(cor.mat,3)    
                  } else {
                    cor.mat <- round(cov_to_cor(object@add),3)    
                  }
                } else {
                  #one location
                  Va <- as.numeric(object@add)*object@meanG
                  Vr <- as.numeric(object@g.resid)
                  
                  V1 <- matrix(c(Va,Vr),ncol=1)
                  rownames(V1) <- c("additive","g.resid")
                  
                  if (n.mark > 0) 
                    V1 <- rbind(object@fixed.marker.var,V1)
                }
              } else {
                
                #no markers
                if (nrow(object@g.resid) > 1) {
                  #multiple locations
                  Vg <- mean(object@g.resid[upper.tri(object@g.resid,diag=F)])
                  VgL <- mean(diag(object@g.resid)) - Vg
                  
                  V1 <- matrix(c(Vg,VgL),ncol=1)
                  rownames(V1) <- c("genotype","g x loc")
                  
                  cor.mat <- round(cov_to_cor(object@g.resid),3)
                } else {
                  Vg <- as.numeric(object@g.resid)
                  V1 <- matrix(c(Vg),ncol=1)
                  rownames(V1) <- c("genotype")
                }
              }
            
              if (!is.na(object@meanOmega)) {
                V2 <- matrix(c(mean(object@resid),object@meanOmega),ncol=1)
                rownames(V2) <- c("g x env","residual")
              } else {
                V2 <- matrix(mean(object@resid),ncol=1)
                rownames(V2) <- "residual"
              }
              vtab <- rbind(V1,V2)
              pvar <- round(vtab / sum(vtab),3)
              colnames(pvar) <- "Var. Proportion"
              
            } #single trait
            
            if (!is.null(cor.mat)) {
              return(list(prop.var=pvar,cor=cor.mat))
            } else {
              return(pvar) 
            }
          })
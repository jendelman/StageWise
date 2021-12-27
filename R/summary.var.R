#' Summarize variances and correlations
#'
#' Summarize variances and correlations
#' 
#' When reporting the partitioning of variance, the variance component for the additive effect is multiplied by the mean diagonal of the G matrix. 
#' 
#' @param object object of \code{\link{class_var}}
#'
#' @return A matrix with the proportion of variance and covariance (for multi-trait). For multi-location or multi-trait analysis, a correlation matrix is also returned. If a G matrix was included, the additive correlation is above the diagonal, and genotypic correlation below. 
#'
#' @include class_var.R
#' @name summary
#' @exportMethod summary
NULL

setGeneric("summary")
setMethod("summary",c(object="class_var"),
          definition=function(object){
            cor.mat <- NULL
            n.mark <- dim(object@fixed.marker.cov)[3]
            if (n.mark > 0)
              fix.eff.markers <- dimnames(object@fixed.marker.cov)[[3]]
            meanG <- (length(object@meanG) > 0)
            meanOmega <- sum(dim(object@meanOmega)) > 0
            if (ncol(object@resid)>1) {
              #multiple traits
              traits <- rownames(object@resid)
              n.trait <- length(traits)
              tmp <- expand.grid(1:n.trait,1:n.trait)
              tmp <- tmp[tmp$Var1 <= tmp$Var2,]
              tmp$band <- tmp$Var2-tmp$Var1
              tmp <- tmp[order(tmp$band,tmp$Var1),]
              
              tnames <- apply(tmp[,1:2],1,function(z){
                 if(diff(z)==0) {
                   return(traits[z[1]])
                 } else {
                   return(paste(traits[z],collapse=":"))
                 }
               })
              
              ix <- cbind(tmp$Var1,tmp$Var2)
              vtab <- rbind(object@g.resid[ix], object@resid[ix])
              colnames(vtab) <- tnames
              if (!meanG) {
                if (meanOmega) {
                  vtab <- rbind(vtab,object@meanOmega[ix])
                  rownames(vtab) <- c("genotype","g x env","residual")
                } else {
                  rownames(vtab) <- c("genotype","residual")
                }
                cor.mat <- cov_to_cor(object@g.resid)
              } else {
                vtab <- rbind(object@add[ix]*object@meanG,vtab)
                #colnames(vtab) <- traits
                if (meanOmega) {
                  vtab <- rbind(vtab,object@meanOmega[ix])
                  rownames(vtab) <- c("additive","g.resid","g x env","residual")
                } else {
                  rownames(vtab) <- c("additive","g.resid","residual")
                }
                
                cov.mat <- object@add 
                if (n.mark > 0) {
                  vtab2 <- NULL
                  for (i in 1:n.mark) {
                    vtab2 <- rbind(vtab2,object@fixed.marker.cov[cbind(ix,i)])
                  }
                  rownames(vtab2) <- dimnames(object@fixed.marker.cov)[[3]]
                  colnames(vtab2) <- colnames(vtab)
                  vtab <- rbind(vtab2,vtab)
        
                  cov.mat <- cov.mat + apply(object@fixed.marker.cov,c(1,2),sum)
                } 
                cor.mat <- cov_to_cor(cov.mat)
                tmp <- cov_to_cor(cov.mat+object@g.resid)
                cor.mat[lower.tri(cor.mat)] <- tmp[lower.tri(tmp)]
                cor.mat <- round(cor.mat,3)
              }
              
              #vtab <- as.matrix(vtab)
              pvar <- round(vtab %*% diag(1/apply(vtab,2,sum)),3)
              colnames(pvar) <- colnames(vtab)
            } else {
              
              #single trait
              if (meanG) {
                if (nrow(object@add) > 1) {
                  #multiple locations
                  tmp <- partition(object@add)
                  Va <- tmp[1]*object@meanG
                  VaL <- tmp[2]*object@meanG
                  
                  Vr <- mean(object@g.resid[,1])
                  V1 <- matrix(c(Va,VaL,Vr),ncol=1)
                  rownames(V1) <- c("additive","add x loc","g.resid")
                  
                  cov.mat <- object@add 
                  if (n.mark > 0) {
                    tmp <- matrix(NA,nrow=2*n.mark,ncol=1)
                    rownames(tmp) <- as.vector(matrix(rbind(fix.eff.markers,paste(fix.eff.markers,"x loc")),ncol=1))
                    for (i in 1:n.mark) {
                      tmp[c(2*i-1, 2*i),1] <- partition(object@fixed.marker.cov[,,i])
                    }
                    V1 <- rbind(tmp,V1)
                    cov.mat <- cov.mat + apply(object@fixed.marker.cov,c(1,2),sum)
                  } 
                  cor.mat <- cov_to_cor(cov.mat)
                  tmp <- cov_to_cor(cov.mat+tcrossprod(sqrt(object@g.resid)))
                  cor.mat[lower.tri(cor.mat)] <- tmp[lower.tri(tmp)]
                  cor.mat <- round(cor.mat,3)
                  
                } else {
                  #one location
                  Va <- as.numeric(object@add)*object@meanG
                  Vr <- as.numeric(object@g.resid)
                  
                  V1 <- matrix(c(Va,Vr),ncol=1)
                  rownames(V1) <- c("additive","g.resid")
                  
                  if (n.mark > 0) {
                    V1 <- rbind(matrix(object@fixed.marker.cov[1,1,],ncol=1),V1)
                    rownames(V1) <- replace(rownames(V1),1:n.mark,fix.eff.markers)
                  }
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
            
              if (meanOmega) {
                V2 <- matrix(c(mean(object@resid),object@meanOmega),ncol=1)
                rownames(V2) <- c("g x env","residual")
              } else {
                V2 <- matrix(mean(object@resid),ncol=1)
                rownames(V2) <- "residual"
              }
              vtab <- rbind(V1,V2)
              pvar <- round(vtab / sum(vtab),3)
              colnames(pvar) <- "R2"
              
            } #single trait
            
            if (!is.null(cor.mat)) {
              return(list(R2=pvar,cor=cor.mat))
            } else {
              return(pvar) 
            }
          })